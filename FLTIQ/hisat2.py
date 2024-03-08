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

"""Library to run hisat2 from Python."""

import glob
import os
import subprocess
from typing import Any, List, Mapping, Optional, Sequence

from absl import logging

from FLTIQ import utils


class Hisat2:
    """Python wrapper of the Hisat2 binary."""

    def __init__(self,
                 *,
                 hisat2_idx: str,
                 n_cpu: int = 4,
                 output_dir: str = "result_hisat2"):
        """Initializes the Python hisat2 wrapper.

    Args:
        hisat2_idx: The path to the hisat2 index.
        n_cpu: The number of CPUs to give hisat2.

    Raises:
        RuntimeError: If hisat2 binary not found within the path.
    """
        self.hisat2_idx = hisat2_idx

        for reads_path in self.reads:
            if not glob.glob(reads_path + '_*'):
                logging.error('Could not find hisat2 reads %s', reads_path)
                raise ValueError(f'Could not find hisat2 reads {reads_path}')
        self.n_cpu = n_cpu
        self.output_dir = output_dir

    def run(self, R1_path, R2_path,) -> List[Mapping[str, Any]]:
        """hisat2 maping"""
        sample_name = R1_path.strip().split("/")[-1].split(".")[-2]
        samfile_path = self.output_dir.rstrip("/")+"/"+f"{sample_name}.sam"

        cmd = [
            "hisat2",
            '-p', str(self.n_cpu),
            '-x', self.hisat2_idx,
            '-1', R1_path,
            '-2', R2_path,
            '-S', samfile_path]

        logging.info('Launching subprocess "%s"', ' '.join(cmd))
        process = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        with utils.timing('hisat2'):
            stdout, stderr = process.communicate()
            retcode = process.wait()

        if retcode:
            # Logs have a 15k character limit, so log hisat2 error line by line.
            logging.error('hisat2 failed. hisat2 stderr begin:')
            for error_line in stderr.decode('utf-8').splitlines():
                if error_line.strip():
                    logging.error(error_line.strip())
            logging.error('hisat2 stderr end')
            raise RuntimeError('hisat2 failed\nstdout:\n%s\n\nstderr:\n%s\n' % (
                stdout.decode('utf-8'), stderr[:500_000].decode('utf-8')))

        raw_output = dict(
            output=stdout,
            stderr=stderr,
            samfile_path=samfile_path
        )
        return [raw_output]
