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

## 1. HIVRC：17 个 HIV 16S 研究（DEBIAS-M 的主基准）
- 位置：Synapse [syn18406854](https://www.synapse.org/#!Synapse:syn18406854)（DEBIAS-M 论文引用）；项目 wiki [syn18406803](https://www.synapse.org/Synapse:syn18406803/wiki/589668)（Park 2025 引用，名为 HIVRC）。需要免费注册 Synapse 账号，可能要接受数据使用条款。
- 需要的文件：合并后的 taxa 计数表（OTU 或属水平）+ 样本元数据（研究 ID、HIV 状态；如有，16S 区域、DNA 提取试剂盒）。
- 参考规模：DEBIAS-M 使用 17 个研究、1,032 名受试者。
- 输出：`results/data/real/hivrc_16s/`，`batch` = 研究，`phenotype` = HIV 阳性为 1。每名受试者保留一个样本。建议在元数据中额外保留 `region` 和提取试剂盒列（DEBIAS-M 用它们解释偏倚）。
- 下载（python）：`pip install synapseclient`，然后 `synapse get -r syn18406854`（需要 Personal Access Token）。

## 2. Salim 2022：猪粪宏基因组技术重复
- 论文：Sci Rep 2022, [doi:10.1038/s41598-022-26141-x](https://www.nature.com/articles/s41598-022-26141-x)
- 原始数据：ENA [PRJEB31650](https://www.ebi.ac.uk/ena/browser/view/PRJEB31650)（184 个宏基因组）
- 样本信息（储存温度/时长、冻融次数、建库试剂盒、测序平台、猪编号、重复组）：论文的补充表。
- 设计：2 头猪（P1/P2），21 种技术条件组合，每种 2–3 个技术重复；一半样本加了 8 种 spike-in。
- 定量：MetaPhlAn4 或 sylph，输出种水平 reads 计数。
- 输出：`results/data/mock/salim_pig/`。建议 `batch` = 技术条件（如储存条件），`phenotype` = 猪（P1 = 0，P2 = 1，这是真实生物信号），`replicate` = 重复组。另存一列 spike-in 与否。

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
