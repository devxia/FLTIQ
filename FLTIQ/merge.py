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

"""Merge non-FSM SQANTI isoforms with the reference annotation using AGAT and gffread."""

import csv
from pathlib import Path
import re
import shutil
import tempfile

from FLTIQ.utils import run_tool


def gtf_rows(path):
    with Path(path).open() as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) != 9:
                raise ValueError(f'{path}, line {number} requires 9 GTF columns')
            attrs = dict(re.findall(r'(?:^|;)\s*([^\s;]+)\s+"([^"]*)"', fields[8]))
            yield fields, attrs


def write_row(handle, fields, updates):
    # Preserve the original attributes, including repeated GENCODE tags and exon metadata.
    for key, value in updates.items():
        pattern = rf'((?:^|;)\s*{key}\s+)"[^"]*"'
        fields[8], substitutions = re.subn(pattern, lambda m: f'{m[1]}"{value}"', fields[8])
        if not substitutions:
            fields[8] = fields[8].rstrip('; ') + f'; {key} "{value}";'
    handle.write('\t'.join(fields) + '\n')


class Cleaner:
    def __init__(self, sqanti_classification, sqanti_gtf, known_annotation, output_dir='result_SQANTI', *, historical=False):
        self.historical = historical
        self.sqanti_gtf = Path(sqanti_gtf)
        self.known_annotation = Path(known_annotation).resolve()
        self.output_dir = Path(output_dir).resolve()
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.sample_name = self.sqanti_gtf.stem
        self.noFSM_gtf = self.output_dir / f'{self.sample_name}_noFSM.gtf'
        self.merged_gff = self.output_dir / 'merged_init.gff'
        self.merged_gtf = self.output_dir / 'merged_init.gtf'
        self.final_gtf = self.output_dir / f'{self.sample_name}.final.gtf'
        self.fsm = set()
        self.novel_genes = {}
        with Path(sqanti_classification).open() as handle:
            rows = csv.DictReader(handle, delimiter='\t')
            required = {'isoform', 'structural_category', 'associated_gene'}
            if not required.issubset(rows.fieldnames or []):
                raise ValueError(f'{sqanti_classification} requires columns {sorted(required)}')
            for row in rows:
                if row['structural_category'] == 'full-splice_match':
                    self.fsm.add(row['isoform'])
                else:
                    self.novel_genes[row['isoform']] = row['associated_gene']
        self.gene_names = {}
        self.reference_transcripts = set()
        for _, attrs in gtf_rows(self.known_annotation):
            if 'gene_id' in attrs and 'gene_name' in attrs:
                self.gene_names[attrs['gene_id']] = attrs['gene_name']
            if 'transcript_id' in attrs:
                self.reference_transcripts.add(attrs['transcript_id'])

    def del_FSM(self):
        with self.noFSM_gtf.open('w') as handle:
            for fields, attrs in gtf_rows(self.sqanti_gtf):
                tid = attrs.get('transcript_id')
                if not tid or tid in self.fsm:
                    continue
                if tid not in self.novel_genes:
                    raise ValueError(f'Transcript {tid} is absent from the SQANTI classification')
                if self.historical:
                    handle.write('\t'.join(fields) + '\n')
                    continue
                updates = {}
                if 'ENS' in tid or tid in self.reference_transcripts:
                    updates['transcript_id'] = f'novel_{tid}'
                if 'gene_id' not in attrs or ':' in attrs['gene_id']:
                    updates['gene_id'] = self.novel_genes[tid]
                write_row(handle, fields, updates)

    def merge(self):
        if self.noFSM_gtf.stat().st_size:
            return run_tool(['agat_sp_merge_annotations.pl', '--gff', self.known_annotation,
                             '--gff', self.noFSM_gtf, '--out', self.merged_gff], output=self.merged_gff,
                            cwd=self.output_dir)
        return []

    def convert(self):
        if not self.noFSM_gtf.stat().st_size:
            shutil.copyfile(self.known_annotation, self.merged_gtf)
            return []
        return run_tool(['gffread', self.merged_gff, '-T', '-o', self.merged_gtf], output=self.merged_gtf)

    def rm(self):
        self.merged_gff.unlink(missing_ok=True)

    def clean(self):
        with self.final_gtf.open('w') as handle:
            for fields, attrs in gtf_rows(self.merged_gtf):
                marker_source = fields[8] if self.historical else attrs.get('transcript_id', '')
                if 'IDmodified-transcript' in marker_source:
                    continue
                if self.historical:
                    tid = attrs.get('transcript_id')
                    if not tid or ('|' in tid and tid in self.fsm):
                        continue
                    gid = attrs['gene_id']
                    name = attrs.get('gene_name')
                    if not name and ':' in tid and ':' in gid:
                        gid = self.novel_genes.get(tid, gid)
                    if not name or (gid.startswith('ENS') and '_' not in gid):
                        name = self.gene_names.get(gid) or name or gid
                    fields[8] = f'transcript_id "{tid}"; gene_id "{gid}"; gene_name "{name}";'
                    handle.write('\t'.join(fields) + '\n')
                    continue
                gid = attrs['gene_id']
                name = self.gene_names.get(gid) or attrs.get('gene_name') or gid
                write_row(handle, fields, {'gene_name': name})
        if not self.final_gtf.stat().st_size:
            raise ValueError(f'Empty merged annotation at {self.final_gtf}')


def merge_annotations(classification, corrected_gtf, annotation, output_dir, *, historical=False):
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.merge-', dir=output_dir) as work:
        cleaner = Cleaner(classification, corrected_gtf, annotation, work, historical=historical)
        cleaner.del_FSM()
        cleaner.merge()
        cleaner.convert()
        cleaner.rm()
        cleaner.clean()
        for product in [cleaner.noFSM_gtf, cleaner.merged_gtf, cleaner.final_gtf]:
            product.replace(output_dir / product.name)
        return output_dir / cleaner.final_gtf.name
