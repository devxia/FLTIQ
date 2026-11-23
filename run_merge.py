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

from FLTIQ.merge import Cleaner

FLAGS = flags.FLAGS
flags.DEFINE_string('sqanti_classification', "./",'Path to the SQANTI classification')
flags.DEFINE_string('sqanti_gtf', "./", 'Path to the SQANTI output GTF annotation')
flags.DEFINE_string('known_annotation', "./", 'GTF annotation file from database')
flags.DEFINE_string('output_dir', "./", 'Path to the output directory')

def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')
    Cleaner_wrapper = Cleaner(FLAGS.sqanti_classification, FLAGS.sqanti_gtf, FLAGS.known_annotation, FLAGS.output_dir)
    logging.set_verbosity(logging.INFO)
    Cleaner_wrapper.del_FSM()
    Cleaner_wrapper.merge()
    Cleaner_wrapper.convert()
    Cleaner_wrapper.rm()
    Cleaner_wrapper.clean()


if __name__ == '__main__':
    app.run(main)
