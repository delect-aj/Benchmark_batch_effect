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

## 3. Tourlousse 2021：日本微生物组联盟（JMBC）提取/建库方法比较 ⚠️ 确认**没有**处理好的表，但可行性比原先估计的好
- 论文：Tourlousse et al.，*Microbiome* 9:95，"Validation and standardization of DNA extraction and library construction methods for metagenomics-based human fecal microbiome measurements"，[doi:10.1186/s40168-021-01048-3](https://doi.org/10.1186/s40168-021-01048-3)（全文 [PMC8082873](https://pmc.ncbi.nlm.nih.gov/articles/PMC8082873/)）
- **2026-10-15 复核**：全文只有 **1 个附件**（`40168_2021_1048_MOESM1_ESM.pdf`，54 页：方法 p2–7、图 p8–37、表 p38–53），**没有任何逐样本丰度表**。Data availability 只给 SRA PRJNA650228；文中提到数据也上传到 MOSAIC Standards Challenge 平台（platform.mosaicbiome.com），但该域名现在**解析不了**，平台似已下线。
  （PDF 可从 Springer 的 CDN 直接取，不走 Springer 本站的 cookie 墙：`https://static-content.springer.com/esm/art%3A10.1186%2Fs40168-021-01048-3/MediaObjects/40168_2021_1048_MOESM1_ESM.pdf`）
- 附件 PDF 里**有的是设计和真值**：
  - **Table S1**（p38）= **18 个菌株的 mock 群落真实组成** —— DNA mock 每株 5.0%，cell mock 不等量（0–19.1%），并给出每株的 GenBank 组装号、基因组大小、GC、革兰染色
  - **Table S2**（p39）= 11 种**建库**试剂盒（代号 A–K）
  - **Table S4**（p41）= 9 种**DNA 提取**方案（代号 L–T，含 IHMS 协议 Q）
  - **Table S8** = run 与样本/方案的对应关系（ENA 的 `library_name` 是随机码如 `FLYVT2013F`，**对应关系只在这张表里**）
- 实际规模（ENA 查询，2026-10-15）：**543 个 run / 721 GB**，其中 205 粪便、**163 DNA mock**、**111 cell mock**、64 个 MOSAIC 挑战赛样本。
- **如果要做，可行的做法**（只针对 274 个 mock run，粪便样本没有真值）：
  1. ENA 支持 **HTTP Range**（实测返回 206），mock run 平均 8.8M reads / 1.46 GB；18 个菌株的群落取 ~20 万 read pair 就足够稳定，即每个 run 只需前 ~2%（约 33 MB）→ **合计约 9 GB，按实测 242 KB/s 约 10 小时**，而不是 401 GB
  2. 不需要 MetaPhlAn/GTDB 大数据库：Table S1 给了 18 个菌株的组装号，自建 ~65 Mb 的小参考用 bwa/bowtie2 比对计数即可
  3. 成本在于：**Table S8 得从 PDF 里抠出来**（该 PDF 用子集字体，简单文本抽取失败，需要 pdftotext/OCR）
- 价值：这是唯一能提供 **DNA 提取方案**作为批次变量的数据（9 提取 × 11 建库，人粪 + 已知真值）。`salim_pig` 覆盖了储存条件和建库，但**没有提取方案**，而提取恰恰是文献公认最主要的偏倚来源（Costea 2017、DEBIAS-M）。
- 输出（若实施）：`results/data/mock/jmbc/`，`batch` = 提取方案（或 提取 × 建库），`phenotype` = 样本类型，`truth.tsv` = Table S1 的已知比例。

## 4. MBQC：多实验室 16S 质量控制 ⚠️ **原判断"不可行"是错的**（2026-10-15 复核）
整合 OTU 表确实没了（downloads.ihmpdcc.org 下线，Wayback 也没存档），但**原始数据完全可用，而且很小**。

- 数据类型：**16S 扩增子**（ENA：`library_strategy = AMPLICON`，`library_selection = PCR`），不是宏基因组。
- 原始数据：SRA **SRP047083** —— ENA 实际查询（2026-10-15）：**2,117 个 run / 28.2 GB / 1.54 亿 reads**，平均 13.3 MB、3.3 万 reads/run。不是之前以为的"约 16,500 个样本"那种规模（那是标本×分装的网格数）。按实测 242 KB/s 全量下载约 **33 小时**。
- **已按样本拆分好**（median 3.3 万 reads/run，不是多重混样），不需要自己拆 barcode。
- **批次变量直接编码在 `library_name` 里**：格式为 `<标本>.<处理实验室>.<id>`，例如 `sample10.agoodman.6874060577`。
- 设计（从 ENA 元数据还原）：**13 个处理实验室 × 22 个标本**，接近全交叉：

  | 实验室 | run 数 | | 实验室 | run 数 |
  |---|---|---|---|---|
  | jpetrosino | 470 | | dmills | 176 |
  | pschloss | 282 | | rburk / pturnbaugh / dlittman | 各 94 |
  | rknight | 276 | | ggloor | 92 |
  | rflores / kjones | 各 188 | | agoodman 78 / jravel 66 / dgevers 19 | |

- 还可下载的样本信息（补充协议细节与 mock 真值）：
  - 样本清单 https://www.hmpdacc.org/data/MBQC/mbqc_specimens.xls
  - 分装 https://www.hmpdacc.org/data/MBQC/mbqc_aliquots.xls
  - 处理实验室协议 https://www.hmpdacc.org/data/MBQC/mbqc_handling_protocols.xls
  - mock 群落组成 https://www.hmpdacc.org/data/MBQC/mbqc_artificial.xls

### 实施（2026-10-15，脚本已入库）
全部样本都是 **Illumina、16S V4 区**（2023 MiSeq + 94 HiSeq2500，无 454），所以一条流程即可。22 个标本中：11 Fresh、7 冻干、2 Robogut、**2 个人工菌群**（sample103 口腔、sample104 粪便），并带健康状态/性别/年龄。

| 步骤 | 脚本 | 说明 |
|---|---|---|
| 1. 元数据 | `data/mbqc_manifest.sh` | 从 ENA 生成 `mbqc_manifest.tsv`（run→标本/实验室/URL）和 `mbqc_samples.tsv`（22 个标本属性）。实验室来自 `library_name` 的 `<标本>.<实验室>.<id>` |
| 2. 下载 | `data/mbqc_download.sh` | 只取 forward，2117 个文件 **13.1 GB**，4 条流、可断点续传、按字节数校验 |
| 3. 建库 | `data/mbqc.R`（规则 `mbqc`） | DADA2 + SILVA 138.2 → 属水平 |

**关键设计：按实验室分别去噪。** 各实验室读长差别很大（126 / 151 / 174 / 243 / 251 / 294 bp），统一 truncLen 会把短读长的实验室整个过滤掉（`filterAndTrim` 会丢弃短于 truncLen 的 read）。因此每个实验室各自：
- truncLen 由质量曲线自动决定（抽 3 个文件，取每 cycle 中位 Q ≥ 25 的最后一个 cycle，下限 100）
- 各自 `learnErrors` / `dada` / 去嵌合（这本就是 DADA2 对不同测序批次的规范做法）

代价是跨实验室 ASV 序列长度不同、不可按序列合并，所以**统一聚合到属**（与 MicrobiomeHD 几个数据集一致）；分类只做一次（在全部 ASV 的并集上，加载 45 万条参考是最贵的一步）。每个实验室的 seqtab 和分类结果都存 RDS 断点，中断可续跑。

- 输出：`results/data/mock/mbqc/`，`batch` = 处理实验室，`phenotype` = 标本编号，`replicate` = 标本编号（同一材料被重复测量，RUV 类方法需要），另存 `specimen_type`、`health_status`、`sex`、`age`、`artificial`。
- **不写 `truth.tsv`**：两个人工菌群的已发表组成只有"大致接种环数"（`mbqc_artificial.xls` 的 approx. loopfuls），不是定量比例；而且公开元数据里没有标本编号与 MBQC 内部 ID 的对应，无法可靠匹配。所以本数据集当作**带技术重复的批次效应数据**用，不当已知组成的 mock。
- 该规则**不在 `snakemake data` 的默认目标里**，避免误触发 13 GB 下载。
