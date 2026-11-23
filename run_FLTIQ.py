# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Identify full-length transcripts, merge annotations and quantify paired-end RNA-seq."""

import argparse
import csv
import logging
import math
from pathlib import Path
import sys

from FLTIQ.flair import FLAIR
from FLTIQ.hisat2 import Hisat2, count_expression, verify_featurecounts
from FLTIQ.kallisto import historical_transcriptome, quantify
from FLTIQ.merge import merge_annotations
from FLTIQ.preprocess import prepare_long_reads, prepare_short_reads, restore_read_order
from FLTIQ.sqanti import SQANTI


def input_file(value):
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise ValueError(f'Input file does not exist: {path}')
    return path


def valid_name(value):
    if not value or value in {'.', '..'} or Path(value).name != value:
        raise ValueError(f'Expected a single project or sample name, got {value!r}')
    return value


def read_samples(args):
    if args.samples:
        if args.r1 or args.r2:
            raise ValueError('Choose either --samples or the --r1/--r2 pair')
        sheet = input_file(args.samples)
        with sheet.open() as handle:
            rows = csv.DictReader(handle, delimiter='\t')
            if not {'sample', 'r1', 'r2'}.issubset(rows.fieldnames or []):
                raise ValueError('Sample sheet requires sample, r1 and r2 columns')
            samples = []
            for row in rows:
                if not all(row.get(field) for field in ('sample', 'r1', 'r2')):
                    raise ValueError(f'Incomplete sample sheet row {rows.line_num}')
                samples.append((valid_name(row['sample']), input_file(sheet.parent / row['r1']),
                                input_file(sheet.parent / row['r2'])))
        if not samples or (not args.historical_raw and len({s[0] for s in samples}) != len(samples)):
            raise ValueError('Supply at least one sample and use unique sample names')
        if args.historical_raw and len({path for _, r1, r2 in samples for path in (r1, r2)}) != 2 * len(samples):
            raise ValueError('Each raw mate file must occur exactly once in the sample sheet')
        return samples
    if bool(args.r1) != bool(args.r2):
        raise ValueError('Supply both --r1 and --r2 for paired-end quantification')
    return [(args.proj_name, input_file(args.r1), input_file(args.r2))] if args.r1 else []


def boolean(value):
    if value.lower() in {'true', '1', 'yes'}:
        return True
    if value.lower() in {'false', '0', 'no'}:
        return False
    raise argparse.ArgumentTypeError('Expected True or False')


def raw_options(parser):
    group = parser.add_argument_group('Historical workflow from raw cDNA FASTQ')
    group.add_argument('--historical-raw', action='store_true', help='Reproduce the raw cDNA and paired-end workflow')
    group.add_argument('--pychopper-threads', type=int, help='Pychopper thread count (default: -p)')
    group.add_argument('--pychopper-cutoffs', type=float, nargs='+', metavar='Q',
                       help='One recorded Pychopper cutoff per library in input order (default: autotune)')
    group.add_argument('--long-read-order', help='Optional ordered list of all filtered read IDs, one per line')
    group.add_argument('--cage-peaks', help='Reference CAGE peaks, BED')
    group.add_argument('--polya-motifs', help='Ranked poly(A) motif list')
    group.add_argument('--adapters', help='Trimmomatic paired-end adapter FASTA')
    group.add_argument('--transcript-count-gtf', help='Independent GTF for historical transcript counts')
    group.add_argument('--gene-featurecounts', default='featureCounts', help='Subread 2.0.6 featureCounts')
    group.add_argument('--transcript-featurecounts', help='Subread 2.0.1 featureCounts')
    group.add_argument('--reference-transcripts', help='Reference transcript FASTA for historical Kallisto')


def validate_raw(args, libraries, samples, count, kallisto):
    required = ['cage_peaks', 'polya_motifs']
    if samples:
        required.append('adapters')
    if count:
        required.extend(['transcript_count_gtf', 'transcript_featurecounts'])
    if kallisto:
        required.append('reference_transcripts')
    for key in required:
        if not getattr(args, key):
            raise ValueError(f"The selected historical workflow requires --{key.replace('_', '-')}")
    if args.pychopper_cutoffs and (len(args.pychopper_cutoffs) != len(libraries) or
                                  not all(math.isfinite(q) and q > 0 for q in args.pychopper_cutoffs)):
        raise ValueError('Provide one positive Pychopper cutoff per raw long-read library')
    if args.pychopper_threads is not None and args.pychopper_threads < 1:
        raise ValueError('Pychopper thread count must be positive')
    if args.is_native_rna is False:
        raise ValueError('Historical FLAIR uses --nvrna for alignment of the Pychopper cDNA reads')
    if count:
        verify_featurecounts(args.gene_featurecounts, '2.0.6')
        verify_featurecounts(args.transcript_featurecounts, '2.0.1')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-r', '--reference', required=True, help='Reference genome FASTA')
    parser.add_argument('-a', '--annotation', required=True, help='Reference GTF on the same assembly')
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument('-f', '--fl_reads', nargs='+', action='extend', help='Long reads, one file per library')
    inputs.add_argument('--raw-library', nargs='+', action='append', metavar='FASTQ',
                        help='Ordered FASTQ chunks for one raw library. Repeat for additional libraries.')
    raw_options(parser)
    parser.add_argument('-s', '--sqanti_dir', required=True, help='Directory containing sqanti3_qc.py')
    parser.add_argument('--sqanti-report', choices=['html', 'pdf', 'both', 'skip'], default='html',
                        help='SQANTI3 report format (default: html)')
    parser.add_argument('--sqanti_python', default=sys.executable, help='Python in the SQANTI3 environment')
    parser.add_argument('-o', '--output_dir', required=True)
    parser.add_argument('-n', '--proj_name', default='myProject')
    parser.add_argument('-p', '--n_cpu', type=int, default=4)
    parser.add_argument('--is_native_rna', type=boolean, default=True)
    parser.add_argument('-i', '--hisat2_idx', help='HISAT2 genome index prefix for quantification')
    parser.add_argument('--quant_method', choices=['featurecounts', 'kallisto', 'both'],
                        help='Short-read quantification method (default with samples: featurecounts)')
    parser.add_argument('--bootstraps', type=int, default=100, help='Number of Kallisto bootstrap samples')
    parser.add_argument('-1', '--r1')
    parser.add_argument('-2', '--r2')
    parser.add_argument('--samples', help='TSV with sample, r1 and r2 columns')
    parser.add_argument('-d', '--data_dir', help=argparse.SUPPRESS)  # accepted by the original shell wrapper
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    try:
        valid_name(args.proj_name)
        reference = input_file(args.reference)
        annotation = input_file(args.annotation)
        input_file(Path(args.sqanti_dir) / 'sqanti3_qc.py')
        reads = []
        for group in args.fl_reads or []:
            parts = [group] if Path(group).is_file() else group.replace(',', ' ').split()
            reads.extend(input_file(p) for p in parts)
        libraries = [[input_file(p) for p in chunks] for chunks in args.raw_library or []] or [[p] for p in reads]
        if not libraries or args.n_cpu < 1:
            raise ValueError('Provide long reads and a positive thread count')
        samples = read_samples(args)
        if args.quant_method and not samples:
            raise ValueError('--quant_method requires paired-end reads or a sample sheet')
        method = args.quant_method or 'featurecounts'
        count = bool(samples) and method in {'featurecounts', 'both'}
        kallisto = bool(samples) and method in {'kallisto', 'both'}
        if count and not args.hisat2_idx:
            raise ValueError('HISAT2 / featureCounts requires --hisat2_idx')
        if args.bootstraps < 0:
            raise ValueError('--bootstraps must be nonnegative')
        if args.historical_raw:
            validate_raw(args, libraries, samples, count, kallisto)
        elif any([args.raw_library, args.pychopper_threads is not None, args.pychopper_cutoffs, args.long_read_order, args.adapters, args.reference_transcripts,
                  args.transcript_count_gtf, args.transcript_featurecounts, args.gene_featurecounts != 'featureCounts']):
            raise ValueError('Raw preprocessing and historical quantification options require --historical-raw')
        for key in ('long_read_order', 'cage_peaks', 'polya_motifs', 'adapters', 'transcript_count_gtf', 'reference_transcripts'):
            if getattr(args, key):
                setattr(args, key, input_file(getattr(args, key)))
    except (ValueError, FileNotFoundError) as error:
        parser.error(str(error))

    run_dir = Path(args.output_dir).resolve() / args.proj_name
    run_dir.parent.mkdir(parents=True, exist_ok=True)
    try:
        run_dir.mkdir()
    except FileExistsError:
        parser.error(f'Project directory already exists: {run_dir}. Choose a new project or output directory.')
    if args.historical_raw:
        threads = args.pychopper_threads or args.n_cpu
        pooled = prepare_long_reads(libraries, args.pychopper_cutoffs, run_dir / 'long_qc', threads)
        if args.long_read_order:
            restore_read_order(pooled, args.long_read_order)
        reads = [pooled]
    flair = FLAIR(reference=reference, annotation=annotation, fl_reads=reads, n_cpu=args.n_cpu,
                  is_native_rna=args.is_native_rna, output_dir=run_dir / 'flair', proj_name=args.proj_name)
    flair.align()
    flair.correct(historical=args.historical_raw)
    flair.collapse()
    sqanti = SQANTI(sqanti_dir=args.sqanti_dir, isoforms=flair.isoforms, annotation=annotation,
                    reference=reference, n_cpu=args.n_cpu, output_dir=run_dir / 'sqanti',
                    proj_name=args.proj_name, python=args.sqanti_python)
    sqanti.qc(cage_peaks=args.cage_peaks, polya_motifs=args.polya_motifs, force_id_ignore=not args.historical_raw,
              report=args.sqanti_report)
    merged = merge_annotations(sqanti.classification, sqanti.corrected_gtf, annotation, run_dir / 'merge',
                               historical=args.historical_raw)
    kallisto_samples = samples
    if args.historical_raw and samples:
        samples, kallisto_samples = prepare_short_reads(samples, args.adapters, run_dir / 'short_qc', args.n_cpu)
    if count:
        aligner = Hisat2(hisat2_idx=args.hisat2_idx, n_cpu=args.n_cpu, output_dir=run_dir / 'hisat2')
        bams = [aligner.sorted_bam(r1, r2, name) for name, r1, r2 in samples]
        count_expression(merged, bams, run_dir / 'counts', args.n_cpu, historical=args.historical_raw,
                         gene_tool=args.gene_featurecounts, transcript_tool=args.transcript_featurecounts or 'featureCounts',
                         transcript_annotation=args.transcript_count_gtf)
    if kallisto:
        directory = run_dir / 'kallisto'
        transcripts = None
        if args.historical_raw:
            transcripts = historical_transcriptome(sqanti.corrected_fasta, sqanti.classification,
                                                   args.reference_transcripts, directory)
        quantify(reference, merged, kallisto_samples, directory, args.n_cpu, args.bootstraps, transcripts=transcripts)
    print(merged)


if __name__ == '__main__':
    main()
