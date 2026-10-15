# 需要手动处理的数据集

这些数据没有可直接下载的处理后表格，需要自行下载和定量。处理完后按下面的**统一格式**放到指定目录，工作流会直接使用（已存在的文件不需要对应的 Snakemake 规则）。

## 统一输出格式

每个数据集一个目录，放在 `benchmark/results/data/real/<名字>/`（参与方法评测）或 `benchmark/results/data/mock/<名字>/`（只用于校准）：

| 文件 | 内容 |
|---|---|
| `counts.tsv` | 行 = 样本，第一列 `sample_id`，其余列 = taxa，**整数 reads 计数**（不要相对丰度） |
| `meta.tsv` | 第一列 `sample_id`（顺序与 counts 一致），必需列 `batch`、`phenotype`（二分类用 0/1，0 = 对照）；技术重复数据可加 `replicate` 列（同一生物样本的重复取相同值，RUV 类方法会用） |
| `taxonomy.tsv` | 可选。行 = taxa（与 counts 列名一致），列 = 分类等级（kingdom…genus/species），缺失填 `unclassified` |

过滤规则与其它数据集一致：taxa 流行度 ≥ 10%（全体样本），去掉全零样本。16S 请聚合到**属**，宏基因组用**种**。

参与评测的数据集，把名字加到 `config/config.yaml` 的 `real:` 列表里即可。

---

## 1. HIVRC：17 个 HIV 16S 研究 ✅ 已完成（2026-10-15）
- 来源：Synapse [syn18406854](https://www.synapse.org/#!Synapse:syn18406854)（由用户下载）
- 原始文件已上传到服务器 `benchmark/results/raw/hivrc/`：
  - `insight.merged_genus.txt` — 属水平计数，已抽平到每样本 20,000 reads（1026 属 × 1032 样本）
  - `insight.merged_otus.txt` — Resphera Insight 原始计数，行是物种名（38,124 行），备用
  - `SupplementaryMaterial.xlsx` — 样本清单（Study、SeqID、hivstatus 及年龄/性别/ART/CD4 等）+ 各研究的 NCBI accession
- 构建脚本 `data/hivrc_16s.R`，输出 `results/data/real/hivrc_16s/`，已加入 `config/config.yaml` 的 `real:` 列表
- 结果：**1032 样本 × 160 属 × 17 个研究**，与 DEBIAS-M 基准规模一致

## 2. ~~Salim 2022：猪粪宏基因组技术重复~~ ✅ 已自动化（2026-10-15，**无需自己处理原始 reads**）
原来以为必须从 ENA 下 184 个 run 自己定量，实际不用：**作者已发布 MGmapper 属级计数表**。

- 数据论文：Poulsen, Kaas, Aarestrup, Pamp. *Microbiol Spectr* 2021;9:e01387-21, [doi:10.1128/Spectrum.01387-21](https://doi.org/10.1128/Spectrum.01387-21)（PMC8510183）。我们原先引的 [doi:10.1038/s41598-022-26141-x](https://doi.org/10.1038/s41598-022-26141-x) 是 Fachrul/Méric/Inouye/Pamp/Salim *Sci Rep* 2022 的**方法学论文**，复用的就是这批数据（同一 ENA accession）。
- 三个来源都是小文件，`data/salim_pig.R` 会自动下载（规则 `salim_pig`，原始文件进 `results/raw/salim/`）：
  | 来源 | 文件 | 大小 | 内容 |
  |---|---|---|---|
  | figshare [10.6084/m9.figshare.12010989](https://doi.org/10.6084/m9.figshare.12010989) | `Item_A_TaxonomyRaw.txt` | 2.6 MB | 属级**原始整数计数**，2640 属 × 438 样本 |
  | 同上 | `Item_C_FeatureShort.txt` | 229 KB | 属的 Domain…Genus 谱系 |
  | Europe PMC `PMC8510183/supplementaryFiles` | `spectrum.01387-21-s0001.xlsx`（Table S1） | 6.2 MB（zip） | 逐样本条件：基质、时间、温度、冻融次数、重复组、spike、建库、lane、reads 数 |
- Table S1 只能从 **Europe PMC** 取：PMC 本体给 curl 返回 proof-of-work 挑战页，ASM 返回 403。
- 设计比原先记录的更丰富：**4 种基质**（Pig_feces_1/2、Sewage_1/2，而非 2 头猪）、19 种时间×温度组合（0h 到 8640h ≈ 1 年；−80/−20/5/22 °C）、冻融 0–4 次、**3 种建库**（Nextflex 306 / Kappa 94 / Nextera 38）、12 个 flow-cell lane、**125 组技术重复**、151 spiked + 233 unspiked + 24 纯 mock + 30 提取空白。
- 输出：`results/data/mock/salim_pig/`，`batch` = 储存条件（时间_温度），`phenotype` = 基质，`replicate` = `Replicate_Groups`；另存 `spiked`、`library_prep`、`lane`、`time_h`、`temp_c`、`freeze_thaw`，任一列都可当批次用。`truth.tsv` = 24 个纯 mock 样本的已知细胞比例（Table S3）。

## 3. Tourlousse 2021：日本微生物组联盟（JMBC）提取/建库方法比较
- 论文：Microbiome 9:95, [doi:10.1186/s40168-021-01048-3](https://doi.org/10.1186/s40168-021-01048-3)（全文 [PMC8082873](https://pmc.ncbi.nlm.nih.gov/articles/PMC8082873/)）
- 原始数据：NCBI SRA BioProject [PRJNA650228](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA650228)；样本与 run 的对应关系在论文补充表 **Table S8**，mock 群落组成在 **Table S1**。
- 设计：同一粪便样本和 mock 群落，用多种 DNA 提取方法和建库方法处理（宏基因组）。
- 定量：同上（MetaPhlAn4 / sylph），种水平计数。
- 输出：`results/data/mock/jmbc/`，`batch` = 提取方法（或 提取 × 建库组合），`phenotype` = 样本（粪便 = 0，mock = 1，或按具体样本编号），mock 样本的真实组成另存 `truth.tsv`（行 = mock 样本，列 = 种，值 = 真实比例）。

## 4. MBQC：多实验室 16S 质量控制
- 整合 OTU 表所在服务器（downloads.ihmpdcc.org）已下线。
- 原始数据：SRA **SRP047083**（约 16,500 个样本，15 个处理实验室）。
- 样本信息（仍可下载）：
  - 样本清单 https://www.hmpdacc.org/data/MBQC/mbqc_specimens.xls
  - 分装 https://www.hmpdacc.org/data/MBQC/mbqc_aliquots.xls
  - 处理实验室协议 https://www.hmpdacc.org/data/MBQC/mbqc_handling_protocols.xls
  - mock 群落组成 https://www.hmpdacc.org/data/MBQC/mbqc_artificial.xls
- 建议：不必全部重做。用一条流程（如 DADA2 + SILVA）处理**所有处理实验室**的一个子集，这样批次只来自湿实验差异。
- 输出：`results/data/mock/mbqc/`，`batch` = 处理实验室，`phenotype` = 标本编号（或 粪便 vs 人工群落），`replicate` = 标本编号；mock 标本的真实组成存 `truth.tsv`。
