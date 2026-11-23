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

"""Long-read alignment, splice junction correction and isoform discovery with FLAIR 2."""

from pathlib import Path

from FLTIQ.utils import run_tool


class FLAIR:
    def __init__(self, *, reference, annotation, fl_reads, n_cpu=4, is_native_rna=True,
                 output_dir='result_flair', proj_name='flair'):
        self.reference = reference
        self.annotation = annotation
        self.fl_reads = [fl_reads] if isinstance(fl_reads, (str, Path)) else list(fl_reads)
        self.n_cpu = n_cpu
        self.native_args = ['--nvrna'] if is_native_rna else []
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        self.prefix = Path(output_dir) / proj_name
        self.bed = f'{self.prefix}.bed'
        self.corrected_bed = f'{self.prefix}_all_corrected.bed'
        self.isoforms = f'{self.prefix}.isoforms.gtf'

    def align(self):
        return run_tool(['flair', 'align', '-g', self.reference, '-r', *self.fl_reads,
                         '-t', self.n_cpu, '-o', self.prefix, *self.native_args], output=self.bed)

    def correct(self, *, historical=False):
        return run_tool(['flair', 'correct', '-g', self.reference, '-q', self.bed,
                         '-f', self.annotation, '-t', self.n_cpu, '-o', self.prefix,
                         *([] if historical else self.native_args)], output=self.corrected_bed)

    def collapse(self):
        return run_tool(['flair', 'collapse', '-g', self.reference, '-q', self.corrected_bed,
                         '-r', *self.fl_reads, '-f', self.annotation, '-t', self.n_cpu,
                         '-o', self.prefix], output=self.isoforms)
