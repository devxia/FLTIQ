"""Quantify paired-end RNA-seq against transcripts in the merged annotation."""

from FLTIQ.utils import run_tool


def quantify(reference, annotation, samples, output_dir, n_cpu, bootstraps):
    output_dir.mkdir(parents=True, exist_ok=True)
    transcripts = output_dir / 'transcripts.fa'
    index = output_dir / 'transcripts.idx'
    run_tool(['gffread', annotation, '-g', reference, '-w', transcripts], output=transcripts)
    run_tool(['kallisto', 'index', '-i', index, transcripts], output=index)
    for name, r1, r2 in samples:
        sample_dir = output_dir / 'samples' / name
        sample_dir.mkdir(parents=True, exist_ok=True)
        run_tool(['kallisto', 'quant', '-t', n_cpu, '-b', bootstraps, '-i', index,
                  '-o', sample_dir, r1, r2], output=sample_dir / 'abundance.tsv')
