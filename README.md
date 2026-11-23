# FLTIQ

FLTIQ（Full-Length Transcript Identification and Quantification）从长读长 RNA 测序数据鉴定转录本，生成合并参考注释，并可用配对短读长 RNA-seq 定量。安装后用 `fltiq` 运行主流程，已有 SQANTI3 结果时可用 `fltiq-merge` 单独合并注释。输入和输出均可放在项目目录以外。

## 本地安装

在项目根目录使用 Python 3.9 或更高版本运行以下命令。激活环境后，可从任意工作目录调用两个命令。

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install .
fltiq --version
fltiq --help
fltiq-merge --help
```

开发时可使用 `python -m pip install -e .`。将流程软件安装到另一台机器时，可用 `python -m pip wheel --no-deps -w dist .` 构建 wheel，复制 `.whl` 文件后在目标环境中用 `python -m pip install /path/to/fltiq-0.1.0-py3-none-any.whl` 安装。源码目录中的 `./run_FLTIQ.sh` 和 `./run_merge.sh` 也可继续使用。

## 软件环境

pip 安装 FLTIQ 命令和 Python 模块。驱动程序只依赖 Python 标准库，运行分析时需另行配置所选步骤对应的外部工具。

| 步骤 | 验证时使用的软件 |
| --- | --- |
| 长读长比对、剪接位点校正、转录本折叠 | FLAIR 2.0.0 及其依赖，包括 minimap2、samtools 和 bedtools |
| 转录本结构分类 | SQANTI3 5.2.2 及其官方运行环境，包括报告所需的 R 包 |
| 合并 GTF | AGAT 的 `agat_sp_merge_annotations.pl`、gffread 0.12.7 |
| 短读长比对和计数 | HISAT2 2.2.1、samtools、Subread 2.0.6 的 featureCounts |
| 转录本丰度估计 | gffread、Kallisto 0.50.1 |

原始 FASTQ 复现模式另需 Pychopper 2.7.10（含 HMMER）、Chopper 0.8.0 和 Trimmomatic 0.39（含 Java），注释合并使用 AGAT 1.1.0。运行时将所选步骤的工具加入 `PATH`，源码目录的 shell 入口支持用 `PYTHON` 指定解释器。`-s` 指向包含 `sqanti3_qc.py` 的目录，`--sqanti_python` 指向 SQANTI3 环境的 Python。输入、参考资源、样本表和运行结果均可保存在仓库外。

## 从原始 FASTQ 复现历史分析

`--historical-raw` 接受已完成 basecalling 和拆分文库的 Nanopore cDNA FASTQ。以下命令同时生成合并注释、featureCounts 矩阵和 Kallisto 丰度。环境变量由运行者设为本机的实际路径。

```bash
fltiq --historical-raw \
  -f "$LONG_FASTQ" -r "$GENOME_FASTA" -a "$REFERENCE_GTF" \
  -s "$SQANTI_DIR" --sqanti_python "$SQANTI_PYTHON" \
  --cage-peaks "$CAGE_BED" --polya-motifs "$POLYA_MOTIFS" \
  --samples "$LANES_TSV" --adapters "$ADAPTER_FASTA" \
  --quant_method both -i "$HISAT2_INDEX_PREFIX" \
  --gene-featurecounts "$FEATURECOUNTS_206" \
  --transcript-featurecounts "$FEATURECOUNTS_201" \
  --transcript-count-gtf "$TRANSCRIPT_COUNT_GTF" \
  --reference-transcripts "$REFERENCE_TRANSCRIPTS_FASTA" \
  -o "$OUTPUT_DIR" -n analysis -p 20 --pychopper-threads 40
```

| 步骤 | 复现模式的输入与行为 |
| --- | --- |
| 全长 cDNA 预处理 | `-f` 后按预期合并顺序列出各文库的 FASTQ。一个文库分成多个 FASTQ 时可改用 `--raw-library chunk1.fq.gz chunk2.fq.gz`，每个文库分别指定一次。每个文库独立运行 Pychopper（PCS109），输出报告、统计、未分类和 rescued reads。将主输出经 Chopper（平均质量至少 7、长度至少 500 nt）过滤，按文库顺序合成一个 FASTQ。 |
| FLAIR 和 SQANTI3 | FLAIR `align` 使用 `--nvrna`，`correct` 使用参考 GTF，`collapse` 使用校正后的 BED。SQANTI3 使用外部 CAGE BED 和 poly(A) motif 表，保留默认的 ID 检查和 HTML 报告。极小的验证集缺少某些分类时，SQANTI3 5.2.2 的 R 报告可能失败，可显式选择 `--sqanti-report skip` 来验证分类和定量。 |
| 注释合并 | 先从校正 GTF 去掉所有 `full-splice_match`，保留其余记录原有的 ID，与参考注释经 AGAT 和 gffread 合并。最后过滤 `IDmodified-transcript` 和残留的带竖线 FSM，按参考基因名及 SQANTI3 的 `associated_gene` 补全属性，输出 `transcript_id`、`gene_id`、`gene_name`。独立合并入口支持 `--historical`。 |
| 双端短读长 | 样本表含制表符分隔的 `sample`、`r1`、`r2` 三列，每行一批。相同 `sample` 的行按表内顺序处理，初测在前、补测在后。每批分别用 Trimmomatic PE、TruSeq3-PE 接头、`ILLUMINACLIP:<接头文件>:2:30:10:2:true SLIDINGWINDOW:10:20 MINLEN:150` 清洗，按样本分别拼接所有批次的成对 R1、R2。 |
| 基因与转录本计数 | HISAT2 比对合并批次并按坐标排序。基因矩阵用本次的最终 GTF 和 Subread 2.0.6 的 `-p -t exon -g gene_id`。转录本矩阵用外部指定的独立 GTF 和 Subread 2.0.1 的 `-p -t exon -g transcript_id`。这两个版本的 `-p` 分别对应 read 和片段计数。要重建补测前的矩阵，可用仅含首批记录的样本表单独运行。 |
| Kallisto | 按 SQANTI3 分类保留校正 FASTA 中的非 FSM 转录本，再接上原版参考转录本 FASTA。以此建索引，每个样本使用首批清洗后的配对 reads，默认 100 次 bootstrap。 |

Pychopper 默认按文库自动抽样调参。复跑已有分析时，`--pychopper-cutoffs` 可按文库顺序接收各自 `statistics.tsv` 中记录的 `Cutoff(q)`，避免重新抽样改变入选 reads。要对照历史数值，还需使用原版基因组、参考 GTF、CAGE BED、motif 表、接头、HISAT2 索引，以及所选定量分支使用的转录本 GTF 或参考转录本 FASTA。独立转录本 GTF 属于历史定量的外部输入。

并行预处理可能改变同一文库内的 read 顺序。若需要复现历史 FASTQ 的排列，可以在仓库外从当时合并的四行式 FASTQ 导出 ID，再传入 `--long-read-order "$ORDER_TSV"`。流程会检查每个过滤后的 read ID 恰好出现一次，并在运行 FLAIR 前按给定顺序重排。顺序表只保存 ID，序列和质量值仍从本次的原始 FASTQ 重新生成。

```bash
awk 'NR % 4 == 1 {print substr($1, 2)}' "$ARCHIVED_POOLED_FASTQ" > "$ORDER_TSV"
```

只运行长读长鉴定和合并时，保留长读长、参考基因组、参考 GTF、SQANTI3、CAGE 和 motif 的参数，省略短读长及定量参数。短读长单独选择 `--quant_method featurecounts` 或 `kallisto` 时，仅需对应分支的参考文件和工具。

## 从已质控的长读长生成合并注释

参考基因组 FASTA 和参考 GTF 应来自同一组装，染色体命名须一致。`-f` 可接收一个或多个长读长 FASTA/FASTQ，多个文件将共同生成一套转录本。直接 RNA 测序使用默认值 `--is_native_rna=True`，cDNA 数据使用 `False`。以下示例中的输入和输出路径需换成实际位置。

```bash
fltiq \
  -r /refs/genome.fa -a /refs/annotation.gtf \
  -f /reads/long_1.fastq.gz /reads/long_2.fastq.gz \
  -s "$SQANTI_DIR" --sqanti_python "$SQANTI_PYTHON" \
  --is_native_rna=False -p 12 -o /results -n analysis
```

此模式接受已质控的长读长。该命令依次运行 FLAIR `align`、`correct`、`collapse` 和 SQANTI3 `qc`。合并步骤保留参考注释中的已知转录本，从 SQANTI3 校正结果中取出分类不等于 `full-splice_match` 的转录本，再用 AGAT 和 gffread 合并。若全部为完整剪接匹配，则以参考注释生成最终 GTF。每次运行使用新的 `<输出目录>/<项目名>`。

## 为已质控的 reads 加入短读长定量

在同一条主流程命令中追加下列参数。单个样本使用 `-1 /reads/R1.fastq.gz -2 /reads/R2.fastq.gz`。多样本使用 `--samples /path/to/samples.tsv`，每个样本一行并以制表符分隔。相对路径从样本表所在目录解析。

```tsv
sample	r1	r2
sample_A	reads/sample_A_R1.fastq.gz	reads/sample_A_R2.fastq.gz
sample_B	reads/sample_B_R1.fastq.gz	reads/sample_B_R2.fastq.gz
```

| 追加参数（以样本表为例） | 生成的结果 |
| --- | --- |
| `--samples samples.tsv -i /refs/hisat2/genome` | 默认定量。HISAT2 比对到基因组，samtools 排序，featureCounts 按基因和转录本各生成一张计数矩阵 |
| `--samples samples.tsv --quant_method kallisto` | gffread 从同一份最终 GTF 提取转录本序列，Kallisto 建索引并逐样本估计丰度和 TPM |
| `--samples samples.tsv -i /refs/hisat2/genome --quant_method both` | 同时运行两种定量方法 |

`-i` 填 HISAT2 索引前缀。需要新建索引时，运行 `hisat2-build /refs/genome.fa /refs/hisat2/genome`。Kallisto 默认运行 100 次 bootstrap，可用 `--bootstraps` 调整。

featureCounts 以外显子为计数特征。历史分析用 Subread 2.0.6 按 `gene_id` 计单条 read，用 2.0.1 按 `transcript_id` 计片段。已质控模式在 2.0.6 中对转录本显式指定 `--countReadPairs` 以保留这两个计数单位，并对两张表统一使用本次生成的最终 GTF。多个转录本共享外显子时，无法唯一归属的片段会出现在转录本计数的 `.summary` 中。

## 主要输出

下表中的路径均相对于 `<输出目录>/<项目名>/`。

| 路径 | 内容 |
| --- | --- |
| `long_qc/pooled.fq`、`short_qc/<样本名>/` | 原始 FASTQ 模式下的全长过滤结果和分批短读长清洗结果 |
| `flair/<项目名>.isoforms.gtf` | 长读长支持的转录本结构 |
| `sqanti/<项目名>_classification.txt`、`sqanti/<项目名>_corrected.gtf` | 转录本分类及校正注释 |
| `merge/<项目名>_corrected.final.gtf` | 已知和新转录本的合并注释，用于基因计数；原始 FASTQ 复现模式按各分支使用历史指定的参考资源 |
| `hisat2/<样本名>.bam` | 按坐标排序的短读长比对 |
| `counts/geneCounts.tsv`、`counts/transcriptCounts.tsv` | 多样本原始计数矩阵，另附分配情况 `.summary` |
| `kallisto/samples/<样本名>/abundance.tsv` | 每个转录本的估计片段数及 TPM |

## 单独合并已有 SQANTI3 结果

```bash
fltiq-merge \
  -i /results/analysis/sqanti/analysis_classification.txt \
  -g /results/analysis/sqanti/analysis_corrected.gtf \
  -a /refs/annotation.gtf -d /results/analysis_merge
```

合并入口与主流程调用同一份源码。此示例的最终文件为 `/results/analysis_merge/analysis_corrected.final.gtf`。

### Cite FLTIQ <a name="cite"></a>

If you use or discuss FLTIQ, please cite the following [paper](https://www.nature.com/):
> Binbin Xia, Jianghua Shen, et al. 2024. The alternative splicing landscape of infarcted mouse heart identifies isoform-level targets for therapeutic intervention. to be published
