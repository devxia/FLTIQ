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

from absl import app
from absl import flags
from absl import logging

from FLTIQ.flair import FLAIR
from FLTIQ.sqanti import SQANTI
from FLTIQ.hisat2 import Hisat2

flags.DEFINE_string('proj_name', "myProject", 'Project will be subdirectory name store the results.')
flags.DEFINE_string('sqanti_dir', "./", 'Path to the SQANTI directory.')
flags.DEFINE_string('output_dir', "./", 'Path to a directory that will store the results.')
flags.DEFINE_string('reference', "./", 'Reference genome (Fasta format)')
flags.DEFINE_string('annotation', "./", 'GTF annotation file')
flags.DEFINE_string('fl_reads', "./", 'Raw reads in fasta or fastq format. This argument accepts '
                                      'multiple (comma/space separated) files.')
flags.DEFINE_boolean('is_native_rna', True, 'Whether is native RNA.')
flags.DEFINE_integer('n_cpu', 4, 'Set number of CPUs')
flags.DEFINE_string('hisat2_idx', "./", 'Path to hisat2 index.')
flags.DEFINE_string('r1', "./", 'Path to NGS R1.')
flags.DEFINE_string('r2', "./", 'Path to NGS R2.')

FLAGS = flags.FLAGS

FLAIR_wrapper = FLAIR(FLAGS.reference, FLAGS.annotation, FLAGS.fl_reads, FLAGS.n_cpu, FLAGS.is_native_rna,
                      FLAGS.output_dir)
SQANTI_wrapper = SQANTI(FLAGS.sqanti_dir, f"{FLAGS.proj_name}.isoforms.gtf", FLAGS.annotation, FLAGS.reference,
                        FLAGS.n_cpu, FLAGS.output_dir)
Hisat2_wrapper = Hisat2(FLAGS.hisat2_idx, FLAGS.n_cpu, FLAGS.output_dir)


def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')
    logging.set_verbosity(logging.INFO)
    FLAIR_wrapper.align()
    FLAIR_wrapper.correct()
    FLAIR_wrapper.collapse()
    SQANTI_wrapper.qc()
    Hisat2_wrapper.run(FLAGS.r1, FLAGS.r2)


if __name__ == '__main__':
    app.run(main)
