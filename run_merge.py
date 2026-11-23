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

"""Merge SQANTI3 corrected isoforms with a reference GTF."""

import argparse
import logging

from FLTIQ.merge import merge_annotations


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--sqanti_classification', required=True)
    parser.add_argument('-g', '--sqanti_gtf', required=True)
    parser.add_argument('-a', '--known_annotation', required=True)
    parser.add_argument('-d', '--output_dir', required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    print(merge_annotations(args.sqanti_classification, args.sqanti_gtf, args.known_annotation, args.output_dir))


if __name__ == '__main__':
    main()
