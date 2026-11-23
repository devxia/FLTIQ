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

"""Map paired-end RNA-seq reads and count exons grouped by gene or transcript."""

from pathlib import Path
import re
import subprocess

from FLTIQ.utils import run_tool


class Hisat2:
    def __init__(self, *, hisat2_idx, n_cpu=4, output_dir='result_hisat2'):
        self.hisat2_idx = hisat2_idx
        self.n_cpu = n_cpu
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def run(self, R1_path, R2_path, sample_name=None):
        name = sample_name or Path(R1_path).name.split('.')[0]
        sam = self.output_dir / f'{name}.sam'
        result = run_tool(['hisat2', '-p', self.n_cpu, '-x', self.hisat2_idx,
                           '-1', R1_path, '-2', R2_path, '-S', sam], output=sam)
        result[0]['samfile_path'] = str(sam)
        return result

    def sorted_bam(self, r1, r2, sample_name):
        sam = Path(self.run(r1, r2, sample_name)[0]['samfile_path'])
        bam = sam.with_suffix('.bam')
        run_tool(['samtools', 'sort', '-@', self.n_cpu, '-o', bam, sam], output=bam)
        sam.unlink()
        return bam


def verify_featurecounts(executable, expected_version):
    result = subprocess.run([executable, '-v'], check=True, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if re.search(r'featureCounts v' + re.escape(expected_version) + r'(?:\s|$)', result.stdout) is None:
        raise ValueError(f'Expected featureCounts {expected_version} at {executable}, got {result.stdout.strip()}')


def count_expression(annotation, bams, output_dir, n_cpu, *, gene_tool='featureCounts',
                     transcript_tool='featureCounts', transcript_annotation=None, historical=False):
    output_dir.mkdir(parents=True, exist_ok=True)
    targets = [(gene_tool, 'gene_id', annotation, 'geneCounts.tsv'),
               (transcript_tool, 'transcript_id', transcript_annotation or annotation, 'transcriptCounts.tsv')]
    for executable, attribute, gtf, filename in targets:
        output = output_dir / filename
        # Subread 2.0.1 counts pairs with -p alone. Newer versions require --countReadPairs.
        pair_flags = ['--countReadPairs'] if attribute == 'transcript_id' and not historical else []
        run_tool([executable, '-T', n_cpu, '-p', *pair_flags, '-t', 'exon',
                  '-g', attribute, '-a', gtf, '-o', output, *bams], output=output)
