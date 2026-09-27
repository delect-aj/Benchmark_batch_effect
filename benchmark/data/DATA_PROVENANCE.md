# 各数据集的处理流程

每个数据集分两段：**上游**（原作者从测序数据到特征表做了什么）和**本项目**（`data/` 下脚本做了什么）。所有数据可用 `snakemake data` 重建。

---

## 宏基因组：crc_mgx / ibd_mgx / ici_mgx

**上游** — curatedMetagenomicData 3.14.0：所有研究用**同一条流程**重新处理，因此跨研究可比。
1. 原始 WGS reads（各研究 SRA/ENA）
2. **MetaPhlAn 3**（mpa_v30 ChocoPhlAn 标记基因库）→ 种水平**相对丰度**
3. 功能谱用 HUMAnN 3（本项目未使用）
4. 元数据由 cMD 策展人统一标准化（`study_condition`、`subject_id`、`sequencing_platform`、`DNA_extraction_kit` 等）

**本项目**（`data/cmd.R`）
1. `counts = TRUE`：cMD 用 相对丰度 ÷ 100 × 该样本 reads 数 得到**估计计数**，再取整
2. 选择队列与样本：CRC/IBD 用 `study_condition`（病例 vs 对照），ICI 用 `ORR`（应答 vs 无应答）
3. 每个受试者只保留一个样本：`days_from_first_collection` 最早的一个（HMP2 等纵向队列需要）
4. 过滤：taxa 流行度 ≥ 10%，去掉全零样本
5. 元数据保留 `platform` 和 `extraction_kit`（DEBIAS-M 发现提取试剂盒解释了大部分处理偏倚），分类表来自 `rowData`，缺失等级填 `unclassified`

> 📌 crc_mgx 中**提取试剂盒完全嵌套在研究里**（每个队列只用一种：Gnome 3 个、MoBio 2 个、Qiagen 2 个、YachidaS_2019 未记录），因此无法把试剂盒效应与研究效应分开；但可以按试剂盒给研究分组，作为另一种批次定义来检验"偏倚是否由提取方案驱动"（DEBIAS-M 的核心论点）。

> ⚠️ **注意**：这不是真实测序计数，而是由相对丰度反推的估计计数。对假设真计数分布的方法（ComBat-seq、ConQuR、RUV-III-NB）这是已知的近似，论文中需说明。

---

## 16S：crc_16s / ibd_16s / hiv_16s

**上游** — MicrobiomeHD（Duvallet et al. 2017, Zenodo 569601）：28 个病例-对照研究用**同一条流程**重新处理。
1. 原始 16S reads（各研究 SRA）
2. 质量过滤、去重
3. **100% 相似度 de novo OTU 聚类**（文件名中的 `otu_table.100.denovo`）
4. **RDP classifier** 赋予分类（`RDP/*.rdp_assigned`，行名即完整谱系）

**本项目**（`data/microbiomehd.R`）
1. 每个研究取 `*.otu_table.100.denovo.rdp_assigned`（部分研究还有 dbOTU 表，不用）
2. 样本过滤沿用原作者配置：`hiv_lozupone` 只取 `time_point = 1`；`hiv_noguerajulian` 只取 `cohort ∈ {BCN0, STK}`；`ibd_gevers_2014` 的 Zenodo 表本身只含粪便样本
3. 病例/对照映射：CRC = CRC vs H（丢弃腺瘤 `nonCRC`）；IBD = CD+UC vs H/nonIBD；HIV = HIV vs H
4. 每样本 ≥ 100 reads（原作者的阈值）
5. **按 RDP 谱系聚合到属**：各研究的 OTU 是各自 de novo 聚类的，编号不可跨研究比较，属水平才可比
6. 过滤：属流行度 ≥ 10%，去掉全零样本

> ⚠️ 各研究的 16S 区域、测序平台、引物不同（如 hiv_dinh 是 454 V3-V5，crc_baxter 是 MiSeq V4）——这正是要校正的批次效应来源。

---

## seed_16s（混杂案例）

**上游** — Foxx & Rivers 2025（Mendeley Data 10.17632/5xrfg5dym6.1）
1. 8 个种子微生物组研究的原始 16S reads
2. QIIME2 流程得到 **ASV 特征表**（作者代码依赖 `qiime2R`/`biomformat` 读取 QIIME2 artifact；细节在其 Methods S1）
3. 作者发布 `F&R_original_count.csv`（ASV 计数）、`F&R_original_abundance.csv`（相对丰度）、`F&R_metadata.csv`

**本项目**（`data/seed_16s.R`）
1. 用发布的**计数表**，按 `SampleID` 与元数据合并
2. `batch` = 研究，`phenotype` = 宿主植物（**多分类**，非病例-对照），并保留 `region`（扩增子区域）
3. 过滤：**在至少一个研究中**流行度 ≥ 10%（不是全体样本）——各研究扩增区域不同（V4 / V3-V4 / V5-V7 / V6-V8），ASV 几乎不重叠，用统一规则只剩 2/9 个研究

> ⚠️ ASV 为序列哈希，无分类注释；已发表的计数表中没有 Barret 2015 和 Liu 2019。

---

## brooks mock（仅用于校准）

**上游** — Brooks et al. 2015 的细胞 mock 群落，由 McLaren et al. 2019 (eLife) 整理进 `metacal` 包：7 种阴道菌种按**已知比例**混合，分 6 块板测序。

**本项目**（`data/brooks_mock.R`）：`counts.tsv` = 观测 reads，`truth.tsv` = 真实比例，`meta.tsv` 中 `batch` = 板号。用于检验"每批次乘性偏倚"这一模拟假设是否符合真实协议效应。

---

## 未能获取的数据（2026-09-27 复核）

| 数据集 | 是否有处理好的表 | 原因 |
|---|---|---|
| HIVRC（17 个 HIV 研究） | **有**（Resphera Insight 统一注释） | 在 Synapse syn18406854，需登录账号；服务器直接访问返回 403 |
| Salim 2022 猪粪 | 无 | 论文补充材料只有样本元数据和 spike-in 汇总表；需从 ENA 下载 184 个 run（461 GB）自行定量 |
| Tourlousse 2021 JMBC | 无 | 补充材料只有一个 PDF；原始数据 SRA PRJNA650228（543 run，721 GB） |
| MBQC | **已失效** | 整合 OTU 表所在服务器下线，Wayback 也没有存档（只存了 5 个 xlsx 元数据表）；只剩原始数据 SRP047083 |

**带宽实测**（从计算节点，2026-09-27）：ENA FTP ≈ 153 KB/s，ENA https ≈ 242 KB/s，NCBI S3 ≈ 27 KB/s，Zenodo ≈ 80 KB/s。按最快的 242 KB/s 计算，461 GB 需要约 22 天——因此**在本集群上重新处理原始数据不可行**，除非换用国内镜像或在其它网络环境下载后拷贝进来。
