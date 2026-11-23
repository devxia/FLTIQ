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

"""Common utilities for data pipeline tools."""

import contextlib
import logging
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile
import time
from typing import Optional


@contextlib.contextmanager
def tmpdir_manager(base_dir: Optional[str] = None):
    """Context manager that deletes a temporary directory on exit."""
    tmpdir = tempfile.mkdtemp(dir=base_dir)
    try:
        yield tmpdir
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


@contextlib.contextmanager
def timing(msg: str):
    logging.info('Started %s', msg)
    tic = time.time()
    yield
    toc = time.time()
    logging.info('Finished %s in %.3f seconds', msg, toc - tic)


def run_tool(command, *, output=None, env=None, cwd=None):
    """Run one pipeline tool and check its promised output, including silent failures."""
    args = list(map(str, command))
    logging.info('Running %s', shlex.join(args))
    with timing(args[0]):
        result = subprocess.run(args, check=True, stdout=subprocess.PIPE, env=env, cwd=cwd)
    if output is not None and not (Path(output).is_file() and Path(output).stat().st_size):
        raise RuntimeError(f'{args[0]} produced no output at {output}')
    return [{'output': result.stdout, 'stderr': result.stderr}]
