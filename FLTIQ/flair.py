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

"""Library to run FLAIR from Python."""

import glob
import os
import subprocess
from typing import Any, List, Mapping, Optional, Sequence

from absl import logging

from FLTIQ import utils


class FLAIR:
    """Python wrapper of the FLAIR binary."""

    def __init__(self,
                 *,
                 # binary_path: str,
                 reference: str,
                 annotation: str,
                 fl_reads: Sequence[str],
                 n_cpu: int = 4,
                 is_native_rna=True,
                 output_dir: str = "result_flair"):
        """Initializes the Python FLAIR wrapper.

    Args:
        binary_path: The path to the FLAIR executable.
        reference: Reference genome in fasta format.
            Flair will minimap index this file unless there already is a .mmi file in the same location.
        annotation: GTF annotation file.
        fl_reads: Raw reads in fasta or fastq format. This argument accepts multiple (comma/space separated) files.
        n_cpu: The number of CPUs to give FLAIR.

    Raises:
        RuntimeError: If FLAIR binary not found within the path.
    """
        # self.binary_path = binary_path
        self.reference = reference
        self.annotation = annotation
        self.fl_reads = fl_reads
        for reads_path in self.fl_reads:
            if not glob.glob(reads_path + '_*'):
                logging.error('Could not find FLAIR reads %s', reads_path)
                raise ValueError(f'Could not find FLAIR reads {reads_path}')
        self.n_cpu = n_cpu
        self.output_dir = output_dir
        self.is_native_rna = is_native_rna

    def align(self, ) -> List[Mapping[str, Any]]:
        """FLAIR align"""

        reads_cmd = []
        for db_path in self.fl_reads:
            reads_cmd.append('-r')
            reads_cmd.append(db_path)
        cmd = [
            # self.binary_path,
            "flair",
            'align',
            '-g', self.reference,
            '-t', str(self.n_cpu),
            '-o', self.output_dir]
        if self.is_native_rna:
            cmd += ['--nvrna']
        cmd += reads_cmd

        logging.info('Launching subprocess "%s"', ' '.join(cmd))
        process = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        with utils.timing('FLAIR align'):
            stdout, stderr = process.communicate()
            retcode = process.wait()

        if retcode:
            # Logs have a 15k character limit, so log FLAIR align error line by line.
            logging.error('FLAIR align failed. FLAIR align stderr begin:')
            for error_line in stderr.decode('utf-8').splitlines():
                if error_line.strip():
                    logging.error(error_line.strip())
            logging.error('FLAIR align stderr end')
            raise RuntimeError('FLAIR align failed\nstdout:\n%s\n\nstderr:\n%s\n' % (
                stdout.decode('utf-8'), stderr[:500_000].decode('utf-8')))

        raw_output = dict(
            output=stdout,
            stderr=stderr)
        return [raw_output]

    def correct(self, ) -> List[Mapping[str, Any]]:
        """FLAIR correct"""

        cmd = [
            self.binary_path,
            'correct',
            '-q', self.reference,
            '-g', self.reference,
            '--gtf', self.annotation,
            '-t', str(self.n_cpu),
            '-o', self.output_dir]

        logging.info('Launching subprocess "%s"', ' '.join(cmd))
        process = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        with utils.timing('FLAIR correct'):
            stdout, stderr = process.communicate()
            retcode = process.wait()

        if retcode:
            # Logs have a 15k character limit, so log FLAIR correct error line by line.
            logging.error('FLAIR correct failed. FLAIR correct stderr begin:')
            for error_line in stderr.decode('utf-8').splitlines():
                if error_line.strip():
                    logging.error(error_line.strip())
            logging.error('FLAIR correct stderr end')
            raise RuntimeError('FLAIR correct failed\nstdout:\n%s\n\nstderr:\n%s\n' % (
                stdout.decode('utf-8'), stderr[:500_000].decode('utf-8')))

        raw_output = dict(
            output=stdout,
            stderr=stderr)
        return [raw_output]

    def collapse(self, ) -> List[Mapping[str, Any]]:
        """FLAIR collapse"""

        reads_cmd = []
        for db_path in self.reads:
            reads_cmd.append('-r')
            reads_cmd.append(db_path)
        cmd = [
            self.binary_path,
            'collapse',
            '-g', self.reference,
            '-q', self.reference,
            '--gtf', self.annotation,
            '-t', str(self.n_cpu),
            '-o', self.output_dir]
        cmd += reads_cmd

        logging.info('Launching subprocess "%s"', ' '.join(cmd))
        process = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        with utils.timing('FLAIR collapse'):
            stdout, stderr = process.communicate()
            retcode = process.wait()

        if retcode:
            # Logs have a 15k character limit, so log FLAIR collapse error line by line.
            logging.error('FLAIR collapse failed. FLAIR collapse stderr begin:')
            for error_line in stderr.decode('utf-8').splitlines():
                if error_line.strip():
                    logging.error(error_line.strip())
            logging.error('FLAIR collapse stderr end')
            raise RuntimeError('FLAIR collapse failed\nstdout:\n%s\n\nstderr:\n%s\n' % (
                stdout.decode('utf-8'), stderr[:500_000].decode('utf-8')))

        raw_output = dict(
            output=stdout,
            stderr=stderr)
        return [raw_output]
