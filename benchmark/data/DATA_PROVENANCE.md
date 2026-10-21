# 各数据集的处理流程

每个数据集分两段：**上游**（原作者从测序数据到特征表做了什么）和**本项目**（`data/` 下脚本做了什么）。所有数据可用 `snakemake data` 重建。

> 📋 本文件讲**怎么来的**；一行一个数据集的结构化清单（路径、规模、批次/表型、真值、协变量、注意事项）在 **[`data/DATASETS.tsv`](DATASETS.tsv)**，由 `data/inventory.sh` 生成。

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

## hivrc_16s（HIV，17 个研究）

**上游** — HIV Re-analysis Consortium（HIVRC，Tuddenham et al. 2020 CID；Synapse syn18406854）：22 个 HIV 肠道微生物组研究中的 **17 个**用**同一条流程**重新处理。
1. 各研究原始 16S reads（accession 见 `SupplementaryMaterial.xlsx` 的 "Raw Data Accession Numbers" 页）
2. **Resphera Insight** 统一做物种级分类注释
3. 发布两张表：物种级**原始计数**（`insight.merged_otus.txt`，行是物种名或无法区分的物种组）和属级表（`insight.merged_genus.txt`）
4. 属级表**每个样本抽平到 20,000 reads**（各列合计恰为 20000），行名是完整谱系

**本项目**（`data/hivrc_16s.R`）
1. 用**属级抽平计数表**（与 crc_16s/ibd_16s/hiv_16s 的属水平一致），按 `SeqID` 与样本清单合并（1032 个样本全部匹配）
2. `batch` = 研究（17 个），`phenotype` = `hivstatus`（1 = HIV 阳性），并保留 Age / gender / msm / ARTuse / vl / cd4 等临床协变量
3. 过滤：属流行度 ≥ 10%，去掉全零样本 → **1032 样本 × 160 属**
4. 分类表由谱系拆分为 kingdom…genus

> 📌 **天然部分混杂**：Perez-Santiago（0/13）和 Villar-Garcia（0/24）两个研究**没有 HIV 阴性对照**，Serrano-Villar 2017 只有 2 个对照。因此 Percentile Normalization 这类要求每批次都有对照的方法会在此数据上失败——这正是文献中"混杂设计"问题的真实案例。
>
> ⚠️ 属级表已抽平到 20,000 reads，样本间测序深度差异被移除；假设原始计数的方法（ComBat-seq、ConQuR）拿到的是抽平后的计数。物种级原始计数表（未抽平）也已上传，如需保留深度信息可改用它。
>
> 📊 未校正时批次效应很强：batch R²（Aitchison）= **0.305**，与文献中 HIVRC 的 0.12–0.30 一致，是目前所有真实数据集里最强的。

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

## salim_pig（猪粪 + 污水储存实验，仅用于校准）

**上游** — Poulsen et al. 2021 *Microbiol Spectr* 9:e01387-21（丹麦 DTU，Pamp 组；ENA PRJEB31650）。Fachrul et al. 2022 *Sci Rep*（"Salim 2022"）复用的就是这批数据。
1. 4 种基质（猪粪 ×2、污水 ×2）分装后按 19 种时间×温度组合储存（0h–8640h，−80/−20/5/22 °C），冻融 0–4 次，每种条件 2–3 个技术重复；一半分装加入 8 个菌株的合成 mock 群落，另有 24 个纯 mock 和 30 个提取空白
2. 测序：Illumina HiSeq 4000，**3 种建库试剂盒**（Nextflex / Kappa / Nextera），12 个 flow-cell lane
3. BBduk2 质量剪切 → **MGmapper**（BWA-MEM，best mode，比对 11 个参考库：bacteria / bacteria_draft / archaea / fungi / virus / human / HumanMicrobiome / 动植物 / 寄生虫等）
4. 计数**除以 2**（按 proper pair 计数），**聚合到属**，发布在 figshare 10.6084/m9.figshare.12010989

**本项目**（`data/salim_pig.R`，规则 `salim_pig`）
1. 用 figshare 的**原始整数计数表** `Item_A`（2640 属 × 438 样本），分类谱系用 `Item_C`（两张表属的排序不同，按属名 join，1:1 对上 2640 个）
2. 逐样本条件取论文 **Table S1**（从 Europe PMC 的 supplementaryFiles zip 里取；PMC 本体对 curl 返回 proof-of-work 页，ASM 返回 403），按 `Sample` 列 join
3. 去掉 30 个提取空白（`SpikedUnspiked == "EC"`）；计数表里还有 1 个来自后续批次（`DTU2018_…`）、Table S1 中没有的样本，一并丢弃
4. 只保留 **Bacteria + Archaea**，外加 mock 的两个真核成员（*Saccharomyces*、*Cryptosporidium*）—— MGmapper 还会比对人、动植物参考库，那些不属于群落
5. `batch` = 储存条件（`时间h_温度`），`phenotype` = 基质（4 类），`replicate` = `Replicate_Groups`；另存 `spiked`、`library_prep`、`lane`、`time_h`、`temp_c`、`freeze_thaw`
6. 过滤：属流行度 ≥ 10%，去掉全零样本 → **408 样本 × 1537 属，19 种储存条件，125 组技术重复**
7. `truth.tsv` = 24 个纯 mock 样本的已知**细胞比例**（Table S3 的 "added cells to samples"；Pig_feces_1 的 *Propionibacterium* 和 *Staphylococcus* 只加了其它三种基质的一半）。与 brooks mock 一致，基因组大小的影响留在被测的偏倚里 —— 这正是 McLaren 乘性模型要刻画的东西

> 📊 **直接验证了模拟用的批次模型**。纯 mock 样本里 per-taxon 的 log2（观测/真值）偏倚：
>
> | 属 | log2 偏倚 | 基因组 (Mb) | | 属 | log2 偏倚 | 基因组 (Mb) |
> |---|---|---|---|---|---|---|
> | *Salmonella* | +2.34 | 4.96 | | *Cryptosporidium* | −3.41 | 9.10 |
> | *Bacteroides* | +1.95 | 5.24 | | *Staphylococcus* | −4.00 | 2.82 |
> | *Escherichia* | +1.61 | 5.20 | | *Saccharomyces* | −5.88 | 12.16 |
> | *Fusobacterium* | +1.58 | 2.17 | | *Propionibacterium* | −6.19 | 2.65 |
>
> 跨度约 **360 倍**，方向符合裂解难易（革兰阳性的 *Propionibacterium*/*Staphylococcus* 和带厚壁的真核严重低估，易裂解的革兰阴性高估）。关键是：**这个偏倚向量在 3 种建库方式之间高度一致（r = 0.992–0.993）**，说明"每 taxon 一个乘性因子"的假设成立，且偏倚主要由提取/裂解决定、与建库方式基本无关（DEBIAS-M 的核心论点）。基因组大小只解释其中约 1 个 log2 单位（细菌内部大基因组确实偏高），真核成员基因组最大却最低估，说明裂解效应远大于基因组大小效应。

> 📊 **实测（2026-10-15，cu01）**：`phenotype`（基质）R² = 0.640（Aitchison）/ 0.692（Bray），压倒性 —— 猪粪与污水本就是完全不同的群落，过校正会毁掉它，正好用来检测过校正。技术批次相对温和：储存条件 R² = 0.068 / 0.059，建库 0.032 / 0.014，lane 0.039 / 0.009。**技术重复确实更接近**：同一 `Replicate_Groups` 内平均 Aitchison 距离 35.4，跨组 67.2（940 对同组），所以"重复应一致"这个真值可用。
>
> ⚠️ 这是 MGmapper 的参考库比对计数，不是 MetaPhlAn 标记基因丰度，与 crc_mgx 等 cMD 数据集不可直接比较；只用于校准，不进方法评测列表。作者另提供按基因组大小归一化并总和标度的表（figshare `Item_B`），本项目不用。

---

## brooks mock（仅用于校准）

**上游** — Brooks et al. 2015 的细胞 mock 群落，由 McLaren et al. 2019 (eLife) 整理进 `metacal` 包：7 种阴道菌种按**已知比例**混合。

**本项目**（`data/brooks_mock.R`）：`counts.tsv` = 观测 reads，`truth.tsv` = 真实比例，`meta.tsv` 中 `batch` = 板号，`phenotype` = 混入的菌种数（1/2/3/4/7），另存 `species_list`。用于检验"每批次乘性偏倚"这一模拟假设是否符合真实协议效应。

> ⚠️ **只有 2 块板、80 个样本**（此前文档写的"6 块板"有误）：metacal 打包的是细胞混合子集，不是 Brooks 2015 的全部数据。`Mixture_type` 在这个子集里恒为 "Cells"，因此不保留。批次数少意味着"每批次偏倚"只能在两块板之间比较，结论强度有限 —— 宏基因组一侧的同类检验用 salim_pig（3 种建库 × 24 个纯 mock）更有力。

---

## mbqc（多实验室 16S 质控，仅用于校准/批次研究）

**上游** — MBQC（Microbiome Quality Control，Sinha et al. 2017）：同一批标本分装后寄给多个"处理实验室"，各自做 DNA 提取与 16S V4 测序。**联盟发布的整合 OTU 表所在服务器已下线且无存档**，因此上游只剩原始 reads（SRA SRP047083）。

**本项目**（`data/mbqc_manifest.sh` → `data/mbqc_download.sh` → `data/mbqc.R`，规则 `mbqc`）
1. 从 ENA 还原设计：`library_name` 的格式是 `<标本>.<处理实验室>.<id>`，由此得到 13 个实验室 × 22 个标本；标本属性（Fresh/冻干/Robogut/人工菌群、健康状态、性别、年龄）取自 ENA 的 BioSample 属性
2. 只下载 **forward reads**（2117 个文件 / 13 GB）；各实验室测序配置不同（单端到 2×294），反向读无法一致重叠
3. **按实验室分别去噪**（DADA2 对不同测序批次的规范做法）：truncLen 由质量曲线自动决定（中位 Q ≥ 25 的最后一个 cycle），并以**读长 10 分位数封顶**；`maxEE = 2` 全局统一，通过率 < 50% 的实验室自动缩短 truncLen 重试
4. 各实验室独立 `learnErrors` / `dada` / 去嵌合；ASV 跨实验室长度不同、不可按序列合并，故**统一聚合到属**（与 MicrobiomeHD 一致）
5. 分类只做一次：在全部 ASV 的并集（49,765 条）上用自建的 SILVA 138.2 训练集，约 1 小时
6. 过滤：属流行度 ≥ 10%，去掉全零样本 → **2034 样本 × 181 属**

> 📊 **完全交叉设计**：22 个标本每一个都被全部 13 个实验室测过，批次与表型**零混杂**（Cramér's V = 0）——这是混杂梯度的零点参照。批次 R²（Aitchison）= **0.099**，标本 R² = **0.481**；同标本跨实验室距离 26.0 vs 不同标本 38.3（10 万组同标本样本对）。

> ⚠️ 不提供 `truth.tsv`：两个人工菌群标本（sample103 口腔、sample104 粪便）的已发表组成只有"大致接种环数"，不是定量比例，且公开元数据没有标本编号与 MBQC 内部 ID 的对应。259 个人工菌群样本以 `artificial = 1` 标记，可按需使用。

> ⚠️ 6 个 fastq 在源头即格式损坏，已跳过并在 data/MANUAL_DATASETS.md 记名。

---

## jmbc_extraction / jmbc_library（提取与建库方案对比）

**上游** — Tourlousse et al. 2021 *Microbiome* 9:95（日本微生物组联盟，SRA PRJNA650228）。论文**只发布了一个 54 页 PDF**，没有任何特征表；ENA 的元数据也没有方案信息（`library_name` 是随机码，BioSample 属性全是占位值）。因此全部自行重建。

**本项目**（`data/jmbc_tables.py` → `jmbc_ref.sh` → `jmbc_download.sh` → `jmbc_quantify.sh` → `jmbc.R`）
1. 用 `pdftotext -layout` 从补充 PDF 抠出两张表并入库（集群上无需 poppler 即可复现）：
   - **Table S1** → 20 个菌株的 GenBank 组装号、基因组大小、革兰型、**真值比例**（DNA mock 均等 5%；cell mock 不等量，其中 *Megamonas funiformis* 为 0）
   - **Table S8** → 每个 run 属于哪个实验、哪种方案。注：pdftotext 偶尔把两字符方案码拆开（`C protocol L` = CL），已在解析时归一
2. 只取 phase I 的两个方案对比实验（164 run），其余是 phase II 的操作员/批次/实验室变化
3. **参考基因组**：按 GCA 号从 ENA 直接取完整组装（含质粒），20 个共约 75 Mb，建 bowtie2 索引 —— 群落成员已知，不需要 MetaPhlAn/GTDB 大库
4. **只下载每个 run 的前 20 万条 read**（ENA 支持 HTTP Range）：对 20 菌株 mock 绰绰有余，下载量从 100 GB 降到 **2.3 GB**
5. 比对后只计 **MAPQ ≥ 10** 的 read：两个 *B. longum* 亚种会产生无法区分的比对，丢弃比让 bowtie2 随机断连更干净
6. 按实验拆成两个数据集（批次轴不同，见下）

> 📊 **观测 vs 真值（cell mock）**：per-taxon log2 偏倚 **−2.09 到 +1.50（约 12 倍）**，方向符合裂解难易 —— *Bifidobacterium*（两亚种）、*Lactobacillus*、*Staphylococcus* 严重低估，*Blautia*、*Bacillus*、*Parabacteroides* 高估。批次（提取方案）R² = **0.143**，材料（mock vs 粪便）R² = **0.976**。
>
> 📌 **一个与 salim_pig 对照的重要结果**：偏倚向量在 9 种提取方案中有 7 种高度一致（r ≈ 0.92–0.98），但**方案 O 和 R 明显偏离**（r ≈ 0.56–0.84）。而 salim_pig 只改变建库方式时 r = 0.992–0.993。也就是说，"每 taxon 一个乘性因子"的 McLaren 模型在建库层面很稳，在提取层面**多数成立但并非普适** —— 这正是模拟里需要同时保留"系统性"与"非系统性"批次效应两种情形的实证依据。

> ⚠️ `jmbc_library` 的 84 个样本是**同一份 DNA mock**，生物信号完全不变，因此是纯 null 场景，只用于检测过校正与假信号，不进方法评测列表。
>
> ⚠️ 163/164 个 run 成功；SRR12996048 多次重试后仍返回空，已跳过。

---

## 未能获取的数据（2026-10-20 复核）

> ✅ **Salim 2022 猪粪已不在此列**：之前判断"只有样本元数据和 spike-in 汇总表"是看错了论文 —— 我们引的 *Sci Rep* 2022 是复用这批数据的方法学论文，真正的数据论文 Poulsen 2021 把属级计数表发在了 figshare（2.6 MB）。**不需要下载 461 GB 原始 reads。**

| 数据集 | 是否有处理好的表 | 原因 |
|---|---|---|
| Tourlousse 2021 JMBC | 无（已复核） | 全文只有 1 个 PDF 附件，含设计与 mock 真值但无逐样本丰度；数据另一个去处 MOSAIC 平台已下线。原始数据 SRA PRJNA650228（543 run / 721 GB，其中 274 个是 mock）。**但只取 mock run 的前 ~2% reads（ENA 支持 Range）约 9 GB / 10 小时即可自行定量**，见 data/MANUAL_DATASETS.md |
| ~~MBQC~~ | ✅ **已自行重建（2026-10-19）** | 整合 OTU 表确已下线，但原始数据可用：取 2117 个 run 的 forward reads（13 GB），用 SILVA 138.2 跑统一 DADA2 流程 → **2034 样本 × 181 属，13 实验室 × 22 标本**。详见下方 mbqc 一节与 data/MANUAL_DATASETS.md |

**带宽实测**（从计算节点，2026-09-27）：ENA FTP ≈ 153 KB/s，ENA https ≈ 242 KB/s，NCBI S3 ≈ 27 KB/s，Zenodo ≈ 80 KB/s。按最快的 242 KB/s 计算，461 GB 需要约 22 天——因此**在本集群上重新处理原始数据不可行**，除非换用国内镜像或在其它网络环境下载后拷贝进来。
