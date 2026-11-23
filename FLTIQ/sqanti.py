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

"""Classify FLAIR isoforms with SQANTI3 5.2."""

import os
from pathlib import Path
import shutil
import sys

from FLTIQ.utils import run_tool


class SQANTI:
    def __init__(self, *, sqanti_dir, isoforms, annotation, reference, n_cpu=4,
                 output_dir='result_SQANTI', proj_name='flair', python=sys.executable):
        self.script = Path(sqanti_dir) / 'sqanti3_qc.py'
        self.python = shutil.which(str(python)) or str(python)
        self.isoforms = isoforms
        self.annotation = annotation
        self.reference = reference
        self.n_cpu = n_cpu
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.name = proj_name
        self.classification = self.output_dir / f'{proj_name}_classification.txt'
        self.corrected_gtf = self.output_dir / f'{proj_name}_corrected.gtf'

    def qc(self):
        env = os.environ.copy()
        env['PATH'] = str(Path(self.python).absolute().parent) + os.pathsep + env['PATH']
        result = run_tool([self.python, self.script, self.isoforms, self.annotation, self.reference,
                           '-t', self.n_cpu, '-d', self.output_dir, '-o', self.name,
                           '--force_id_ignore'], output=self.classification, env=env)
        if not self.corrected_gtf.is_file() or not self.corrected_gtf.stat().st_size:
            raise RuntimeError(f'SQANTI3 produced no corrected annotation at {self.corrected_gtf}')
        return result
