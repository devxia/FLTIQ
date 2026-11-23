"""Prepare basecalled cDNA and paired-end FASTQ for the historical workflow."""

from collections import defaultdict
import gzip
from pathlib import Path
import sqlite3
import tempfile

from FLTIQ.utils import run_tool


def concatenate(inputs, output):
    """Stream ordered FASTQ/FASTA files into one uncompressed file."""
    with output.open('wb') as dest:
        needs_newline = False
        for path in inputs:
            opener = gzip.open if str(path).endswith('.gz') else open
            with opener(path, 'rb') as source:
                block = source.read(1024 * 1024)
                if block and needs_newline:
                    dest.write(b'\n')
                while block:
                    dest.write(block)
                    needs_newline = not block.endswith(b'\n')
                    block = source.read(1024 * 1024)
    if not output.stat().st_size:
        raise ValueError(f'No records in {output}')
    return output


def restore_read_order(pooled, order):
    """Reorder four-line Chopper FASTQ by external read IDs, checking the complete set."""
    with tempfile.TemporaryDirectory(dir=pooled.parent) as tmp:
        ordered = Path(tmp) / 'ordered.fq'
        db = sqlite3.connect(Path(tmp) / 'offsets.db')
        try:
            db.execute('CREATE TABLE reads (id BLOB PRIMARY KEY, offset INTEGER) WITHOUT ROWID')
            with pooled.open('rb') as source:
                while True:
                    offset = source.tell()
                    header = source.readline()
                    if not header:
                        break
                    seq, plus, qual = source.readline(), source.readline(), source.readline()
                    if not header.startswith(b'@') or not plus.startswith(b'+') or not qual or not seq:
                        raise ValueError(f'Invalid Chopper FASTQ record at byte {offset}')
                    read_id = header[1:].split()[0]
                    try:
                        db.execute('INSERT INTO reads VALUES (?, ?)', (read_id, offset))
                    except sqlite3.IntegrityError as error:
                        raise ValueError(f'Duplicate FASTQ read ID: {read_id!r}') from error
                with order.open('rb') as ids, ordered.open('wb') as dest:
                    for line in ids:
                        read_id = line.strip()
                        row = db.execute('SELECT offset FROM reads WHERE id=?', (read_id,)).fetchone()
                        if row is None:
                            raise ValueError(f'Unknown or repeated read ID in order file: {read_id!r}')
                        source.seek(row[0])
                        for _ in range(4):
                            dest.write(source.readline())
                        db.execute('DELETE FROM reads WHERE id=?', (read_id,))
                if db.execute('SELECT 1 FROM reads LIMIT 1').fetchone():
                    raise ValueError('Order file is missing reads from the filtered FASTQ')
        finally:
            db.close()
        ordered.replace(pooled)


def prepare_long_reads(libraries, cutoffs, output_dir, threads):
    filtered = []
    for number, chunks in enumerate(libraries, 1):
        folder = output_dir / str(number)
        folder.mkdir(parents=True)
        source = chunks[0] if len(chunks) == 1 else concatenate(chunks, folder / 'input.fq')
        full_length = folder / 'full_length.fq'
        cutoff = ['-q', cutoffs[number - 1]] if cutoffs else []
        run_tool(['pychopper', '-k', 'PCS109', '-t', threads, *cutoff, '-r', folder / 'report.pdf',
                  '-A', folder / 'hits.bed', '-S', folder / 'statistics.tsv',
                  '-u', folder / 'unclassified.fq', '-w', folder / 'rescued.fq', source, full_length])
        passed = folder / 'filtered.fq'
        with passed.open('wb') as handle:
            run_tool(['chopper', '-q', 7, '-l', 500, '-i', full_length], stdout=handle)
        filtered.append(passed)
    return concatenate(filtered, output_dir / 'pooled.fq')


def prepare_short_reads(samples, adapters, output_dir, threads):
    """Keep paired survivors, pool lanes per sample and retain the first lane for Kallisto."""
    lanes = defaultdict(list)
    for sample, r1, r2 in samples:
        folder = output_dir / sample / str(len(lanes[sample]) + 1)
        folder.mkdir(parents=True)
        paired = (folder / 'R1.fq', folder / 'R2.fq')
        run_tool(['trimmomatic', 'PE', '-threads', threads, r1, r2, paired[0],
                  folder / 'unpaired_R1.fq', paired[1], folder / 'unpaired_R2.fq',
                  f'ILLUMINACLIP:{adapters}:2:30:10:2:true', 'SLIDINGWINDOW:10:20', 'MINLEN:150'])
        lanes[sample].append(paired)
    pooled = []
    first = []
    for sample in sorted(lanes):
        pairs = lanes[sample]
        first.append((sample, *pairs[0]))
        mates = pairs[0] if len(pairs) == 1 else tuple(
            concatenate(paths, output_dir / sample / f'pooled_R{mate}.fq')
            for mate, paths in enumerate(zip(*pairs), 1))
        if not all(path.stat().st_size for path in mates):
            raise ValueError(f'No paired reads passed trimming for {sample}')
        pooled.append((sample, *mates))
    return pooled, first
