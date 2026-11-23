"""Quantify paired-end RNA-seq against transcripts in the merged annotation."""

import csv

from FLTIQ.preprocess import concatenate
from FLTIQ.utils import run_tool


def historical_transcriptome(corrected_fasta, classification, reference_transcripts, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    with classification.open() as handle:
        categories = {row['isoform']: row['structural_category']
                      for row in csv.DictReader(handle, delimiter='\t')}
    novel = output_dir / 'non_fsm.fa'
    with corrected_fasta.open() as source, novel.open('w') as dest:
        keep = False
        for line in source:
            if line.startswith('>'):
                if keep:
                    dest.write('\n')
                name = line[1:].split()[0]
                keep = categories[name] != 'full-splice_match'
                if keep:
                    dest.write(f'>{name}\n')
            elif keep:
                dest.write(line.strip())
        if keep:
            dest.write('\n')
    return concatenate([novel, reference_transcripts], output_dir / 'transcripts.fa')


def quantify(reference, annotation, samples, output_dir, n_cpu, bootstraps, *, transcripts=None):
    output_dir.mkdir(parents=True, exist_ok=True)
    if transcripts is None:
        transcripts = output_dir / 'transcripts.fa'
        run_tool(['gffread', annotation, '-g', reference, '-w', transcripts], output=transcripts)
    index = output_dir / 'transcripts.idx'
    run_tool(['kallisto', 'index', '-i', index, transcripts], output=index)
    for name, r1, r2 in samples:
        sample_dir = output_dir / 'samples' / name
        sample_dir.mkdir(parents=True, exist_ok=True)
        run_tool(['kallisto', 'quant', '-t', n_cpu, '-b', bootstraps, '-i', index,
                  '-o', sample_dir, r1, r2], output=sample_dir / 'abundance.tsv')
