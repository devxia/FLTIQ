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

"""Library to run SQANTI from Python."""

import glob
import os
import subprocess
from typing import Any, List, Mapping, Optional, Sequence

from absl import logging

from FLTIQ import utils


class SQANTI:
    """Python wrapper of the SQANTI binary."""

    def __init__(self,
                 *,
                 sqanti_dir: str,
                 isoforms: str,
                 annotation: str,
                 reference: str,
                 n_cpu: int = 4,
                 output_dir: str = "result_SQANTI"):
        """Initializes the Python SQANTI wrapper.

    Args:
        sqanti_dir: The path to the SQANTI directory.
        isoforms: Isoforms (FASTA/FASTQ) or GTF format.
        reference: Reference genome (Fasta format).
        annotation: GTF annotation file.
        n_cpu: The number of CPUs to give SQANTI.

    Raises:
        RuntimeError: If SQANTI binary not found within the path.
    """
        self.sqanti_dir = sqanti_dir
        self.isoforms = isoforms
        self.annotation = annotation
        self.reference = reference

        for reads_path in self.reads:
            if not glob.glob(reads_path + '_*'):
                logging.error('Could not find SQANTI reads %s', reads_path)
                raise ValueError(f'Could not find SQANTI reads {reads_path}')
        self.n_cpu = n_cpu
        self.output_dir = output_dir

    def qc(self, ) -> List[Mapping[str, Any]]:
        """SQANTI Quality Control"""

        cmd = [
            "python",
            self.sqanti_dir.rstrip("/")+"/sqanti3_qc.py",
            '-t', str(self.n_cpu),
            '-d', "self.output_dir",
            self.isoforms,
            self.annotation,
            self.reference]

        logging.info('Launching subprocess "%s"', ' '.join(cmd))
        process = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        with utils.timing('SQANTI'):
            stdout, stderr = process.communicate()
            retcode = process.wait()

        if retcode:
            # Logs have a 15k character limit, so log SQANTI error line by line.
            logging.error('SQANTI failed. SQANTI stderr begin:')
            for error_line in stderr.decode('utf-8').splitlines():
                if error_line.strip():
                    logging.error(error_line.strip())
            logging.error('SQANTI stderr end')
            raise RuntimeError('SQANTI failed\nstdout:\n%s\n\nstderr:\n%s\n' % (
                stdout.decode('utf-8'), stderr[:500_000].decode('utf-8')))

        raw_output = dict(
            output=stdout,
            stderr=stderr)
        return [raw_output]
