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
from pathlib import Path
import sys

from FLTIQ.flair import FLAIR
from FLTIQ.hisat2 import Hisat2, count_expression
from FLTIQ.kallisto import quantify
from FLTIQ.merge import merge_annotations
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
            samples = [(valid_name(row['sample']), input_file(sheet.parent / row['r1']),
                        input_file(sheet.parent / row['r2'])) for row in rows]
        if not samples or len({s[0] for s in samples}) != len(samples):
            raise ValueError('Supply at least one sample and use unique sample names')
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-r', '--reference', required=True, help='Reference genome FASTA')
    parser.add_argument('-a', '--annotation', required=True, help='Reference GTF on the same assembly')
    parser.add_argument('-f', '--fl_reads', required=True, nargs='+', action='extend', help='Long-read FASTA/FASTQ files')
    parser.add_argument('-s', '--sqanti_dir', required=True, help='Directory containing sqanti3_qc.py')
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
        for group in args.fl_reads:
            parts = [group] if Path(group).is_file() else group.replace(',', ' ').split()
            reads.extend(input_file(p) for p in parts)
        if not reads or args.n_cpu < 1:
            raise ValueError('Provide long reads and a positive thread count')
        samples = read_samples(args)
        if args.quant_method and not samples:
            raise ValueError('--quant_method requires paired-end reads or a sample sheet')
        method = args.quant_method or 'featurecounts'
        if samples and method in {'featurecounts', 'both'} and not args.hisat2_idx:
            raise ValueError('HISAT2 / featureCounts requires --hisat2_idx')
        if args.bootstraps < 0:
            raise ValueError('--bootstraps must be nonnegative')
    except ValueError as error:
        parser.error(str(error))

    run_dir = Path(args.output_dir).resolve() / args.proj_name
    run_dir.parent.mkdir(parents=True, exist_ok=True)
    try:
        run_dir.mkdir()
    except FileExistsError:
        parser.error(f'Project directory already exists: {run_dir}. Choose a new project or output directory.')
    flair = FLAIR(reference=reference, annotation=annotation, fl_reads=reads, n_cpu=args.n_cpu,
                  is_native_rna=args.is_native_rna, output_dir=run_dir / 'flair', proj_name=args.proj_name)
    flair.align()
    flair.correct()
    flair.collapse()
    sqanti = SQANTI(sqanti_dir=args.sqanti_dir, isoforms=flair.isoforms, annotation=annotation,
                    reference=reference, n_cpu=args.n_cpu, output_dir=run_dir / 'sqanti',
                    proj_name=args.proj_name, python=args.sqanti_python)
    sqanti.qc()
    merged = merge_annotations(sqanti.classification, sqanti.corrected_gtf, annotation, run_dir / 'merge')
    if samples and method in {'featurecounts', 'both'}:
        aligner = Hisat2(hisat2_idx=args.hisat2_idx, n_cpu=args.n_cpu, output_dir=run_dir / 'hisat2')
        bams = [aligner.sorted_bam(r1, r2, name) for name, r1, r2 in samples]
        count_expression(merged, bams, run_dir / 'counts', args.n_cpu)
    if samples and method in {'kallisto', 'both'}:
        quantify(reference, merged, samples, run_dir / 'kallisto', args.n_cpu, args.bootstraps)
    print(merged)


if __name__ == '__main__':
    main()
