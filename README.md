# FLTIQ

FLTIQ（Full-Length Transcript Identification and Quantification）从长读长 RNA 测序数据鉴定转录本，生成合并参考注释，并可用配对短读长 RNA-seq 定量。主入口为 `run_FLTIQ.sh`，已有 SQANTI3 结果时可单独运行 `run_merge.sh`。输入和输出均可放在项目目录以外。

## 软件环境

驱动脚本使用 Python 3.9 或更新版本的标准库。外部工具按分析步骤配置。

| 步骤 | 验证时使用的软件 |
| --- | --- |
| 长读长比对、剪接位点校正、转录本折叠 | FLAIR 2.0.0 及其依赖，包括 minimap2、samtools 和 bedtools |
| 转录本结构分类 | SQANTI3 5.2.2 及其官方运行环境，包括报告所需的 R 包 |
| 合并 GTF | AGAT 的 `agat_sp_merge_annotations.pl`、gffread 0.12.7 |
| 短读长比对和计数 | HISAT2 2.2.1、samtools、Subread 2.0.6 的 featureCounts |
| 转录本丰度估计 | gffread、Kallisto 0.50.1 |

将所选步骤的可执行文件加入 `PATH`。`-s` 指向包含 `sqanti3_qc.py` 的安装目录，`--sqanti_python` 指定其环境的 Python。下面是在大亨上验证过的设置，可按部署位置调整。

```bash
envs=/software_users/xiabinbin/miniconda3/envs
export PATH="/software_users/xiabinbin/kallisto:$envs/flair_env/bin:$envs/sequence/bin:$PATH"
export PYTHON="$envs/sequence/bin/python"
sqanti_dir=/software_users/xiabinbin/SQANTI3-5.2.2
sqanti_python="$envs/SQANTI3/bin/python"
```

## 从长读长生成合并注释

参考基因组 FASTA 和参考 GTF 应来自同一组装，染色体命名须一致。`-f` 可接收一个或多个长读长 FASTA/FASTQ，多个文件将共同生成一套转录本。直接 RNA 测序使用默认值 `--is_native_rna=True`，cDNA 数据使用 `False`。以下示例中的输入和输出路径需换成实际位置。

```bash
./run_FLTIQ.sh \
  -r /refs/genome.fa -a /refs/annotation.gtf \
  -f /reads/long_1.fastq.gz /reads/long_2.fastq.gz \
  -s "$sqanti_dir" --sqanti_python "$sqanti_python" \
  --is_native_rna=False -p 12 -o /results -n heart
```

该命令依次运行 FLAIR `align`、`correct`、`collapse` 和 SQANTI3 `qc`。合并步骤保留参考注释中的已知转录本，从 SQANTI3 校正结果中取出分类不等于 `full-splice_match` 的转录本，再用 AGAT 和 gffread 合并。若全部为完整剪接匹配，则以参考注释生成最终 GTF。每次运行使用新的 `<输出目录>/<项目名>`。

## 加入短读长定量

在同一条主流程命令中追加下列参数。单个样本使用 `-1 /reads/R1.fastq.gz -2 /reads/R2.fastq.gz`。多样本使用 `--samples /path/to/samples.tsv`，每个样本一行并以制表符分隔。相对路径从样本表所在目录解析。

```tsv
sample	r1	r2
sham_1	reads/sham_1_R1.fastq.gz	reads/sham_1_R2.fastq.gz
mi_1	reads/mi_1_R1.fastq.gz	reads/mi_1_R2.fastq.gz
```

| 追加参数（以样本表为例） | 生成的结果 |
| --- | --- |
| `--samples samples.tsv -i /refs/hisat2/genome` | 默认定量。HISAT2 比对到基因组，samtools 排序，featureCounts 按基因和转录本各生成一张计数矩阵 |
| `--samples samples.tsv --quant_method kallisto` | gffread 从同一份最终 GTF 提取转录本序列，Kallisto 建索引并逐样本估计丰度和 TPM |
| `--samples samples.tsv -i /refs/hisat2/genome --quant_method both` | 同时运行两种定量方法 |

`-i` 填 HISAT2 索引前缀。需要新建索引时，运行 `hisat2-build /refs/genome.fa /refs/hisat2/genome`。Kallisto 默认运行 100 次 bootstrap，可用 `--bootstraps` 调整。

featureCounts 以外显子为计数特征。历史分析用 Subread 2.0.6 按 `gene_id` 计单条 read，用 2.0.1 按 `transcript_id` 计片段。本流程在 2.0.6 中对转录本显式指定 `--countReadPairs` 以保留这两个计数单位，并对两张表统一使用本次生成的最终 GTF。多个转录本共享外显子时，无法唯一归属的片段会出现在转录本计数的 `.summary` 中。

## 主要输出

下表中的路径均相对于 `<输出目录>/<项目名>/`。

| 路径 | 内容 |
| --- | --- |
| `flair/<项目名>.isoforms.gtf` | 长读长支持的转录本结构 |
| `sqanti/<项目名>_classification.txt`、`sqanti/<项目名>_corrected.gtf` | 转录本分类及校正注释 |
| `merge/<项目名>_corrected.final.gtf` | 已知和新转录本的合并注释，也是两个短读长定量分支的注释来源 |
| `hisat2/<样本名>.bam` | 按坐标排序的短读长比对 |
| `counts/geneCounts.tsv`、`counts/transcriptCounts.tsv` | 多样本原始计数矩阵，另附分配情况 `.summary` |
| `kallisto/samples/<样本名>/abundance.tsv` | 每个转录本的估计片段数及 TPM |

## 单独合并已有 SQANTI3 结果

```bash
./run_merge.sh \
  -i /results/heart/sqanti/heart_classification.txt \
  -g /results/heart/sqanti/heart_corrected.gtf \
  -a /refs/annotation.gtf -d /results/heart_merge
```

合并入口与主流程调用同一份源码。此示例的最终文件为 `/results/heart_merge/heart_corrected.final.gtf`。

### Cite FLTIQ <a name="cite"></a>

If you use or discuss FLTIQ, please cite the following [paper](https://www.nature.com/):
> Binbin Xia, Jianghua Shen, et al. 2024. The alternative splicing landscape of infarcted mouse heart identifies isoform-level targets for therapeutic intervention. to be published
