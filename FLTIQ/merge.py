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

"""Library to run merge gtf from Python."""

import glob
import os
import subprocess
from typing import Any, List, Mapping, Optional, Sequence

import pandas as pd
from absl import logging
from gtfparse import read_gtf
from tqdm import tqdm

from FLTIQ import utils


class Cleaner:
    """Python wrapper of the merge gtf."""

    def __init__(self,
                 sqanti_classification: str,
                 sqanti_gtf: str,
                 known_annotation: str,
                 output_dir: str = "result_SQANTI"):
        """Initializes the Python wrapper.

    Args:
        sqanti_classification: The path to the SQANTI classification.
        sqanti_gtf: The path to the SQANTI output GTF annotation.
        known_annotation: GTF annotation file of database.
    Raises:
        RuntimeError: If file not found within the path.
    """
        self.sample_name = None
        self.noFSM_gtf = None
        self.sqanti_classification = sqanti_classification
        self.sqanti_gtf = sqanti_gtf
        self.known_annotation = known_annotation
        self.output_dir = output_dir[:-1] if output_dir.endswith("/") else output_dir

        self.sqanti_table = pd.read_csv(self.sqanti_classification, sep="\t")
        self.FSM_new_isoform_name_all = self.sqanti_table[
            self.sqanti_table.structural_category == "full-splice_match"].isoform.tolist()
        new_col_isoform = []
        for each_isoform in self.sqanti_table["isoform"]:
            if each_isoform in self.FSM_new_isoform_name_all:
                new_col_isoform.append(each_isoform)
            elif "ENS" in each_isoform:
                new_col_isoform.append(f"novel_{each_isoform}")
            else:
                new_col_isoform.append(each_isoform)
        self.sqanti_table["isoform"] = new_col_isoform
        sqanti_table_novel = self.sqanti_table[self.sqanti_table.structural_category != "full-splice_match"]
        self.isoform2geneid_dict = {}
        for each in zip(sqanti_table_novel.isoform, sqanti_table_novel.associated_gene):
            self.isoform2geneid_dict[each[0]] = each[1]

        gtf_df = read_gtf(self.known_annotation).to_pandas()
        self.trans_dict = {}
        self.gene_dict = {}
        self.gene_coding_num_dict = {}
        for i in tqdm(range(len(gtf_df))):
            # transcript_id, transcript_name, transcript_type, gene_id, gene_name, gene_type
            feature = gtf_df.feature.values[i]
            transcript_id = gtf_df.transcript_id.values[i]
            transcript_name = gtf_df.transcript_name.values[i]
            transcript_type = gtf_df.transcript_type.values[i]
            gene_id = gtf_df.gene_id.values[i]
            gene_name = gtf_df.gene_name.values[i]
            gene_type = gtf_df.gene_type.values[i]
            if feature == "transcript" and transcript_type == "protein_coding":
                self.gene_coding_num_dict[gene_id] = self.gene_coding_num_dict.get(gene_id, 0) + 1
            if transcript_id not in self.trans_dict.keys() and feature == "transcript":
                self.trans_dict[transcript_id] = [transcript_name, transcript_type, gene_id, gene_name, gene_type]
            if gene_id not in self.gene_dict.keys() and feature == "gene":
                self.gene_dict[gene_id] = [gene_name]
                self.gene_dict[gene_id].append(gene_type)

    def del_FSM(self,):
        self.sample_name = self.sqanti_gtf.split("/")[-1].rstrip(".gtf")
        self.noFSM_gtf = f"{self.output_dir}/{self.sample_name}_noFSM.gtf"
        with open(self.noFSM_gtf, "w") as f:
            for eachline in tqdm(open(self.sqanti_gtf)):
                each_list = eachline.strip().split("\t")
                p1 = each_list[:8]
                p2 = each_list[8]
                transcript_id = p2.split(";")[0].strip().split(" ")[1].strip('"')
                # 如果是FSM，不保留
                if transcript_id in self.FSM_new_isoform_name_all:
                    continue
                elif "ENS" in transcript_id:
                    new_transcript_id = f"novel_{transcript_id}"
                gene_id = p2.split(";")[1].strip().split(" ")[1].strip('"')

                # f.write(eachline)
                new_p1 = "\t".join(p1)
                f.write(f'{new_p1}\ttranscript_id "{transcript_id}"; gene_id "{gene_id}";\n')


    def merge(self, ) -> List[Mapping[str, Any]]:
        """SQANTI ATAG Control"""
        cmd = [
            "agat_sp_merge_annotations.pl",
            '--gff', self.known_annotation,
            '--gff', self.noFSM_gtf,
            '--out', f"{self.output_dir}/merged_init.gff"]

        logging.info('Launching subprocess "%s"', ' '.join(cmd))
        process = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        with utils.timing('AGAT'):
            stdout, stderr = process.communicate()
            retcode = process.wait()

        if retcode:
            # Logs have a 15k character limit, so log SQANTI error line by line.
            logging.error('AGAT failed. AGAT stderr begin:')
            for error_line in stderr.decode('utf-8').splitlines():
                if error_line.strip():
                    logging.error(error_line.strip())
            logging.error('AGAT stderr end')
            raise RuntimeError('AGAT failed\nstdout:\n%s\n\nstderr:\n%s\n' % (
                stdout.decode('utf-8'), stderr[:500_000].decode('utf-8')))

        raw_output = dict(
            output=stdout,
            stderr=stderr)
        return [raw_output]

    def convert(self, ) -> List[Mapping[str, Any]]:
        """SQANTI ATAG Control"""
        cmd = [
            "gffread",
            f"{self.output_dir}/merged_init.gff",
            "-T",
            '-o', f"{self.output_dir}/merged_init.gtf"]

        logging.info('Launching subprocess "%s"', ' '.join(cmd))
        process = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        with utils.timing('AGAT'):
            stdout, stderr = process.communicate()
            retcode = process.wait()

        if retcode:
            # Logs have a 15k character limit, so log SQANTI error line by line.
            logging.error('AGAT failed. AGAT stderr begin:')
            for error_line in stderr.decode('utf-8').splitlines():
                if error_line.strip():
                    logging.error(error_line.strip())
            logging.error('AGAT stderr end')
            raise RuntimeError('AGAT failed\nstdout:\n%s\n\nstderr:\n%s\n' % (
                stdout.decode('utf-8'), stderr[:500_000].decode('utf-8')))

        raw_output = dict(
            output=stdout,
            stderr=stderr)
        return [raw_output]


    def rm(self, ) -> List[Mapping[str, Any]]:
        """SQANTI ATAG Control"""
        cmd = ["rm", f"{self.output_dir}/merged_init.gff"]

        logging.info('Launching subprocess "%s"', ' '.join(cmd))
        process = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        with utils.timing('AGAT'):
            stdout, stderr = process.communicate()
            retcode = process.wait()

        if retcode:
            # Logs have a 15k character limit, so log SQANTI error line by line.
            logging.error('AGAT failed. AGAT stderr begin:')
            for error_line in stderr.decode('utf-8').splitlines():
                if error_line.strip():
                    logging.error(error_line.strip())
            logging.error('AGAT stderr end')
            raise RuntimeError('AGAT failed\nstdout:\n%s\n\nstderr:\n%s\n' % (
                stdout.decode('utf-8'), stderr[:500_000].decode('utf-8')))

        raw_output = dict(
            output=stdout,
            stderr=stderr)
        return [raw_output]


    def clean(self,):
        with open(f"{self.output_dir}/{self.sample_name}.final.gtf", "w") as f:
            for eachline in tqdm(open(f"{self.output_dir}/merged_init.gtf")):
                # 去除包含FSM误差的列
                if "IDmodified-transcript" in eachline:
                    continue
                each_list = eachline.strip().split("\t")
                p1 = each_list[:8]
                p2 = each_list[8]
                transcript_id = p2.split(";")[0].strip().split(" ")[1].strip('"')
                # 如果是FSM的新转录本，不保留
                if transcript_id in self.FSM_new_isoform_name_all:
                    continue
                gene_id = p2.split(";")[1].strip().split(" ")[1].strip('"')
                try:
                    gene_name = p2.split(";")[2].strip().split(" ")[1].strip('"')
                except IndexError:
                    # print(f"{transcript_id}不存在gene_name，依据gene_id补充")
                    # 新转录本依据gene_ID搭配gene_name
                    if ":" in transcript_id:  # 新转录本
                        if ":" in gene_id:  # 新基因
                            gene_id = self.isoform2geneid_dict[transcript_id]
                            gene_name = gene_id
                        else:  # 只是新转录本
                            gene_name = self.gene_dict[gene_id][0]

                new_p1 = "\t".join(p1)
                if gene_id.startswith("ENS") and ("_" not in gene_id):
                    gene_name = self.gene_dict[gene_id][0]
                f.write(f'{new_p1}\ttranscript_id "{transcript_id}"; gene_id "{gene_id}"; gene_name "{gene_name}";\n')
