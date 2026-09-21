# 微生物组数据批次效应去除方法研究现状：系统性综述（2018–2026）

## TL;DR
- 微生物组批次校正已从借用基因表达工具（ComBat/limma/SVA/RUV）演进到为微生物数据"零膨胀、过度离散、组成性"特征量身定制的专用方法（ConQuR、MMUPHin、PLSDA-batch、DEBIAS-M、MetaDICT等）；截至2026年**尚无公认的金标准方法或统一评价体系**，方法选择必须由研究目标（可视化/关联检验/跨研究预测）和实验设计（是否混杂、是否平衡）决定。
- **性能结论在不同基准间相互矛盾且高度依赖评价指标**：同一方法在Bray-Curtis与Aitchison距离下、PERMANOVA R²与预测AUC之间排名会翻转；例如复合分位数回归在Bray-Curtis下R²最低（最好），但在Aitchison下反被ComBat超过。因此报告多指标结果、避免单一指标结论至关重要。
- 当批次与处理/表型**混杂（confounded）**时（跨研究meta分析中最常见），所有校正方法都会残留批次效应或删除真实生物信号——这是当前最核心的争议与研究空白；推荐优先在实验端做协变量平衡设计、纳入mock community/spike-in，并在混杂场景下把批次作为协变量建模或用meta分析而非强行"校正"。

## Key Findings

**1. 方法生态已分层。** 微生物组批次校正方法可分为五类：(a) 借自基因组学的通用方法（ComBat/ComBat-seq、limma removeBatchEffect、SVA、RUV系列）；(b) 微生物专用统计模型（Percentile Normalization、BDMMA、MMUPHin、ConQuR、PLSDA-batch、复合分位数回归）；(c) 机器学习/领域自适应与深度学习（DEBIAS-M、MetaDICT，以及来自单细胞的autoencoder/VAE/GAN思路）；(d) 基于阴性对照/spike-in/mock community的偏倚校准方法（RUV-III-NB、McLaren一致性偏倚模型/metacal）；(e) 从单细胞迁移的整合方法（Harmony、scVI/scANVI、Scanorama、fastMNN）。

**2. 专用方法的核心突破在于处理"批次不需均匀作用于所有分类单元"（非系统性批次效应）与零膨胀。** ConQuR用两部（logistic + 分位数回归）非参数模型直接在read count上校正整个条件分布；DEBIAS-M用领域自适应学习每个批次每个taxon的乘性偏倚因子；MetaDICT用因果加权+共享字典学习，专门应对未观测混杂与完全混杂设计。

**3. 定量基准存在显著且可解释的矛盾。** 详见Details。关键是评价维度不同：批次分离度（PERMANOVA R²/ANOSIM/silhouette）vs 生物信号保留 vs 跨研究预测AUC；距离度量（Bray-Curtis vs Aitchison vs UniFrac）；数据集混杂程度。

**4. 混杂设计是共识性的局限。** Foxx & Rivers 2025（iMetaOmics）用种子微生物组混杂数据集测试5种算法，全部残留批次效应。

**5. 实验端源头控制正成为并行主线。** 2026年ISME J提出标准化、批次效应无关技术以支持全球协作；mock community与DNA提取/文库标准化（如Tourlousse等2021的日本微生物组联盟JMBC研究）是computational校正的必要补充。

## Details

### 一、通用（借自基因组学）方法

- **ComBat / ComBat-seq**（Johnson et al. 2007, Biostatistics；ComBat-seq Zhang et al. 2020）：经验贝叶斯估计每个特征在每批次的位置-尺度参数。假设批次效应系统性、对数正态（ComBat）或负二项计数（ComBat-seq）。微生物组应用需加伪计数。批次与生物效应不混杂时最有效。
- **limma removeBatchEffect**（Ritchie et al. 2015）：线性模型回归掉批次项，输出用于可视化/聚类，不适合下游正式推断。单变量、不考虑taxa间相关性。
- **SVA（surrogate variable analysis）**（Leek & Storey 2007）：估计代理变量捕捉未观测批次/混杂，纳入线性模型校正下游显著性检验。可处理未观测混杂，但微生物组缺乏"管家基因"参照。
- **RUV系列（RUVg/RUVs/RUV-4/RUV-III）**（Gagnon-Bartsch & Speed 2013等）：用阴性对照特征或技术重复估计并移除未知不需要的变异。需指定阴性对照集与不需要因子数K。

### 二、微生物专用统计模型

- **Percentile Normalization**（Gibbons, Duvallet & Alm 2018, PLoS Comput Biol 14:e1006102）：病例样本每个taxon转换为同研究对照分布的百分位，使对照成为均匀分布后跨研究合并。非参数、模型无关。**仅限病例-对照设计**，需相似的病例/对照定义。Python脚本 + QIIME2插件。原文显示ComBat后研究间仍显著（PERMANOVA p<0.01），而百分位归一化后批次效应消失（PERMANOVA p>0.5）。
- **BDMMA（Bayesian Dirichlet-multinomial regression meta-analysis）**（Dai, Wong, Yu & Wei 2019, Bioinformatics 35(5):807）：贝叶斯meta分析，用Dirichlet-multinomial回归同时建模批次效应并检测表型相关taxa。输入计数/组成型数据。显式建模taxa间依赖与过度离散——作者明确指出"既有校正方法不考虑变量（微生物taxa）间的交互和微生物组数据的过度离散，因此不适用于微生物组数据"。用于关联检测，非通用输出。原文称"BDMMA可成功校正批次效应并大幅减少微生物meta分析中的假发现"，在四个NCBI结直肠癌数据集上验证，"优于所有竞争方法"（混合效应模型、加权Z检验合并P值、ComBat）。
- **MMUPHin**（Ma et al. 2022, Genome Biology 23:208）：把ComBat扩展到零膨胀微生物profile，协变量控制下的批次/研究效应校正 + meta分析差异丰度 + 群体结构发现。假设零膨胀高斯，只适用于相对丰度的某些变换。R/Bioconductor，后端用MaAsLin2。可指定生物协变量以区分技术效应与生物效应。
- **ConQuR**（Ling et al. 2022, Nature Communications 13:5418）：首个"comprehensive"微生物批次去除工具，两部条件分位数回归（存在-缺失用logistic，非零计数用quantile regression for counts），直接在taxa read count上工作，输出校正后的整数计数可用于任意下游分析。非参数、可校正均值/方差/高阶批次效应。需指定参考批次、可加协变量。R包（wdl2459/ConQuR）。
- **PLSDA-batch**（Wang & Lê Cao 2023, Briefings in Bioinformatics 24(2):bbac622）：基于偏最小二乘判别分析的多变量、非参数方法，估计与处理和批次相关的潜在成分后从数据中减去批次成分。变体：weighted PLSDA-batch处理非平衡batch×treatment设计，sparse PLSDA-batch做变量选择避免过拟合。需要处理/表型标签。要求数据预先归一化（如CLR）。输出可用于任意下游分析。R/Bioconductor（PLSDAbatch）。
- **复合分位数回归（Composite Quantile Regression）**（Park & Park 2025, Front Microbiol 16:1484183）：两部模型——负二项回归处理系统性批次效应 + 复合分位数回归（CQR，k=19个分位数）处理OTU层面非系统性批次效应；用Kruskal-Wallis检验+稳健CV（MAD/中位数）选参考批次。作者指出局限：高零膨胀使极端分位数估计不稳，在比率型度量（Aitchison、Canberra）上效果下降。

### 三、机器学习/领域自适应与深度学习

- **DEBIAS-M**（Austin et al. 2025, Nature Microbiology 10(4):897–911；bioRxiv 2024）：领域自适应框架（Domain adaptation with phenotype Estimation and Batch Integration Across Studies of the Microbiome），学习每个taxon每个批次的乘性偏倚校正因子，同时最小化批次差异并最大化跨研究表型关联。考虑组成性。可解释（推断的偏倚因子与具体实验协议强相关）。适合机器学习/跨研究预测；用表型标签但对测试集样本仅用于最小化跨批次差异（避免过拟合）。Python（korem-lab/DEBIAS-M）。
  - 跨研究预测基准（留一研究法，逻辑回归，中位auROC[IQR]）：**HIV**（肠道16S，17研究，1,032受试者）DEBIAS-M 0.70 [0.61–0.75] vs 原始0.50 [0.45–0.53]、ComBat 0.53 [0.49–0.58]、ConQuR 0.54 [0.49–0.57]、voom-SNM 0.53 [0.51–0.63]（所有两两比较DeLong检验p<0.01）；**CRC**（宏基因组）DEBIAS-M 0.76 [0.69–0.78] vs 原始0.58、ComBat 0.62、ConQuR 0.66、voom-SNM 0.56（p<0.01）；**宫颈上皮内瘤变**（5研究，322受试者）DEBIAS-M 0.65 [0.55–0.65] vs 原始0.50、ComBat 0.52、ConQuR 0.56、voom-SNM 0.49。可解释性：DNA提取试剂盒解释了学到偏倚的43%方差（Adonis PERMANOVA p=0.002），16S区域和样本类型额外解释27%和14%（p=0.021、0.039）。
  - **范围说明**：DEBIAS-M论文正文主分类基准明确对比"ComBat、voom-SNM、ConQuR（加原始数据）"三种"流行或新近开发"的方法；MMUPHin、PLSDA-batch在方法附录/第三方基准（如mBatchNet）中与之并列，但非DEBIAS-M主分类基准的对比方法。
- **MetaDICT**（Yuan & Wang 2025, Nature Communications 16:8147）：两阶段——先用因果推断加权法初估批次效应，再用共享字典学习（graph Laplacian利用系统发育/分类学平滑度）精修。相比现有方法，能在存在未观测混杂变量、数据集高度异质、或批次与协变量完全混杂时更好地避免过校正、保留生物变异。可生成taxa和sample层面嵌入用于揭示隐藏结构。用CRC宏基因组和免疫治疗meta分析验证。
- **单细胞可迁移的深度学习**：autoencoder-based batch correction（ABC，Danino et al. 2024, Bioinformatics Advances）、scGen（VAE + 潜空间算术）、scVI/scANVI、CLEAR（对比学习）、DeepBID（负二项autoencoder + 双KL散度）、scMEDAL（混合效应深度autoencoder，分离批次不变与批次特异表示）。这些主要为scRNA-seq开发，在微生物组上有潜力但直接迁移受限于微生物组样本量远小于单细胞、零膨胀与组成性。微生物组内部的深度学习多用于插补/降维（mbSparse用CVAE插补、DeepMicro、VTrans），专门用于批次校正的成熟微生物组深度学习工具仍稀缺。

### 四、基于阴性对照/spike-in/mock community的偏倚校准

- **McLaren一致性偏倚模型 / metacal**（McLaren, Willis & Callahan 2019, eLife 8:e46923）：数学证明并用mock community验证——一致的分类学偏倚产生的fold error在样本间变化，会严重扭曲跨样本比较，甚至导致方向错误的推断；但偏倚可测量、可校正（乘性偏倚结构）。用已知组成的对照样本估计每taxon的相对检测效率。metacal R包实现；fido/pibble贝叶斯扩展。这是"源头可校准"思路的理论基石。
- **RUV-III-NB**（Salim et al. 2022, Nucleic Acids Research 50(16):e96）：负二项/零膨胀负二项GLM直接建模count，用伪重复和阴性对照特征捕捉不需要的变异，输出百分位校正计数（PAC）。原为scRNA-seq开发，已被证明可用于微生物组（用"经验阴性对照taxa"即样本间恒定丰度的taxa）。Salim等2022（Sci Rep, s41598-022-26141-x）在"184份猪粪便微生物组样本、多达21种技术变异组合"上基准测试，识别储存条件/冻融为主要不需要变异来源，并证明RUV-III-NB相较ComBat、ComBat-Seq、RUVg、RUVs"在保留生物信号和移除不需要效应两方面均具比较优势"。
- **实验端定量校准（QMP/spike-in）**：绝对丰度校准（quantitative microbiome profiling）在16S中常需用于使跨研究比较有意义（ISME J 2026提及）。

### 五、从单细胞迁移的整合方法

Harmony（Korsunsky et al. 2019）、scVI/scANVI（Lopez et al. 2018）、Scanorama、fastMNN、Seurat CCA/MNN在单细胞跨数据集整合中表现优异（Harmony/scVI/Scanorama可扩展到百万细胞、运行时间<1小时），在微生物组meta分析中被作为对照或探索性工具使用。MetaDICT论文的比较中，PLSDA-batch和scANVI用欧氏距离、其余方法用Bray-Curtis——这本身说明整合方法输出的是嵌入而非校正后的丰度表，与传统校正方法的可比性有限。

### 六、假设–适用场景对照表

| 方法 (原始出处) | 输入类型 | 分布假设 | 需协变量/表型标签? | 能否处理混杂/非平衡? | 未观测混杂? | 输出可用于任意下游? | 实现 |
|---|---|---|---|---|---|---|---|
| ComBat/ComBat-seq (2007/2020) | 计数(+伪计数)/log | 经验贝叶斯正态/负二项 | 可选协变量 | 弱（假设系统性、不混杂） | 否 | 是 | R (sva/ComBat-seq) |
| removeBatchEffect (2015) | log/CLR | 线性、正态 | 可选 | 弱 | 否 | 限可视化/聚类 | R (limma) |
| SVA (2007) | log | 线性+因子分析 | 表型用于保护 | 部分 | 是（代理变量） | 用于推断校正 | R (sva) |
| RUV系列 (2013+) | 计数/log | 因子模型 | 需阴性对照/重复 | 部分 | 是 | 视变体 | R (RUVSeq) |
| Percentile Norm (Gibbons 2018) | 相对丰度 | 非参数 | 需病例-对照标签 | 仅病例-对照 | 部分 | 是（合并后） | Python/QIIME2 |
| BDMMA (Dai 2019) | 计数/组成 | Dirichlet-multinomial贝叶斯 | 需表型+协变量 | 部分 | 有限 | 仅关联检测 | R |
| MMUPHin (Ma 2022) | 相对丰度 | 零膨胀高斯 | 协变量控制 | 中（协变量调整） | 有限 | 是 | R/Bioc |
| ConQuR (Ling 2022) | read count | 非参数(两部分位数) | 需关键变量+可选协变量 | 中（可加协变量） | 否 | 是（输出计数） | R |
| PLSDA-batch (Wang 2023) | CLR/归一化 | 非参数多变量 | 需处理标签 | 强（weighted变体） | 否 | 是 | R/Bioc |
| 复合分位数回归 (Park 2025) | count | 负二项+CQR | 需批次(参考选择) | 中 | 否 | 是 | R (GitHub) |
| DEBIAS-M (Austin 2025) | 计数/相对丰度 | 乘性偏倚、组成性 | 用表型(防过拟合设计) | 中 | 部分 | 嵌入+校正，适ML | Python |
| MetaDICT (Yuan 2025) | 计数/相对丰度 | 字典学习+因果加权 | 用观测协变量 | 强（完全混杂/未观测） | 是 | 嵌入+校正 | R |
| RUV-III-NB (Salim 2022) | count | (零膨胀)负二项GLM | 需阴性对照+伪重复 | 部分 | 是 | 是（PAC计数） | R |
| metacal/McLaren (2019) | 计数(需mock) | 乘性一致偏倚 | 需mock/已知组成 | 校准非校正 | 否 | 是（校准丰度） | R (metacal) |

### 七、基准评测结论与相互矛盾之处

- **mBatchNet案例研究**（Bioinformatics 2026, btag538；集成12种方法）：在苯酚扰动厌氧消化16S数据（231 OTU、75样本、5个处理日期批次）上，校正前批次ANOSIM≈0.401（p=0.001）；校正后ConQuR、DEBIAS-M、MMUPHin的Bray-Curtis ANOSIM接近0或为负，但MMUPHin和DEBIAS-M按PERMANOVA仍保留可检测的残余批次关联；MetaDICT也降低了ANOSIM统计量（0.098，p=0.013）。结论：需跨多面板解读，无单一分数可判优劣。
- **复合分位数回归论文**（Park & Park 2025）在HIVRC数据上（PERMANOVA R²，越低越好）：Bray-Curtis下本法R²=0.0128（最低/最好）> ConQuR 0.0149 > ComBat 0.0637 > MMUPHin 0.0822 > 百分位归一化0.1191 > 原始0.1199；但Aitchison下本法0.0854反被ComBat 0.0564（最好）超过，ConQuR 0.0901最差；Manhattan下本法0.0065最好；Canberra下MMUPHin 0.0779最好、本法0.1493略优于ConQuR 0.1628。在HPV/MOUTH数据集上本法在所有度量下最低（Bray-Curtis 0.002、Aitchison 0.0246、Canberra 0.0264、Manhattan 0.0029）。**同一数据集、不同距离度量导致方法排名翻转**——这是矛盾结论最清晰的例证，根源是比率型度量（Aitchison/Canberra）对零膨胀敏感。
- **DEBIAS-M基准**用跨研究预测auROC为指标，得出ConQuR/ComBat改善有限甚至不如原始，而以PERMANOVA R²为指标的基准（如ConQuR原文）则显示ConQuR近乎完全移除批次。这说明**"移除批次"与"提升跨研究泛化"是不同目标**：过度移除可能连同真实可迁移信号一起删除。
- **ConQuR原文**（Ling 2022）：CARDIA数据批次效应从5.66%降到0.10%（减少98%）；HIVRC研究变异（原始计数尺度）从30.39%降到1.94%（减少94%），相对丰度尺度从18.66%降到1.63%；HIV随机森林5折CV预测AUC从0.75提升到0.92；且是"唯一在所有场景控制FDR在0.05左右"的方法。

### 八、核心争议

1. **过校正风险**：MetaDICT、Foxx & Rivers等反复强调，当批次与生物变异纠缠时，激进校正会删除真实信号、扭曲效应量、增加假阳/假阴、降低ML准确性。
2. **混杂设计下的局限**：Foxx & Rivers 2025（iMetaOmics 2(3):e70025）用种子微生物组，测试ZMC、Ratio-A、ConQuR、PLSDA、weighted PLSDA五种算法在"目标物种（植物种）不出现在所有批次"的混杂数据上，全部残留批次效应；呼吁"微生物组合成研究的未来需发展稳健处理批次-协变量混杂数据的方法"。
3. **评价缺乏金标准**：真实数据无"无批次的ground truth"（Foxx & Rivers明确指出这是其分析的短板）；PERMANOVA R²、ANOSIM、silhouette、预测AUC测量不同侧面且常给出冲突排名。
4. **校正后统计推断有效性**：removeBatchEffect等把校正后数据当作固定输入再做检验，会低估不确定性、膨胀I类错误；DEBIAS-M批评基于结果变量的校正有过拟合风险；voom-SNM在肿瘤微生物reads上给稀疏特征引入非零值的问题被反复引用为警示（Gihawi et al. 2023指出重大数据分析错误使癌症微生物组发现失效）。
5. **校正 vs 作为协变量建模**：两大流派——"correct"（从数据移除批次变异，输出灵活可做任意下游）vs "account/adjust"（把批次作为协变量放进统计模型，仅限差异丰度检验但推断更诚实）。混杂/小样本时后者更稳健；需可视化/聚类/预测时前者更实用。

### 九、特殊场景

- **纵向/时序数据**：长周期实验中实验室环境本身漂移；RUVg、PLSDA-batch、ConQuR被列为已在微生物组验证的可用方法（RUVg/PLSDA-batch为需预先归一化的线性潜因子模型，ConQuR直接建模原始计数）；有大规模meta纵向研究用Harman校正+MicrobiomeAnalyst评估功能富集（PICRUSt）。
- **多组学联合校正**：Nature Communications 2026多组学Perspective指出即便标准化协议，缓冲液pH、孵育时间、试剂批次、质谱响应的微小差异仍产生批次效应；元数据协调与FAIR是前提。
- **功能通路层面**：MMUPHin可对基因/通路做meta分析；HUMAnN输出的KO/EC/通路丰度可作为校正对象；功能层面校正研究相对薄弱。
- **环境/农业/非人类 vs 人体**：种子/植物微生物组meta分析常缺乏重叠物种、混杂更严重（Foxx & Rivers）；建议用模型物种（拟南芥、水稻、二穗短柄草）作跨研究标准化锚点。
- **宏基因组 vs 16S**：DEBIAS-M、复合分位数回归在两种数据上均验证；16S受引物区、拷贝数影响，宏基因组受数据库/比对影响；ISME J 2026强调即便同一16S协议，储存时间、试剂批号、测序平台也会造成显著差异。

### 十、实验端源头控制与标准化

- **ISME J 2026**（Ge et al. 2026, ISME J 20(1):wrag122）："标准化和批次效应无关技术使微生物组研究的全球协作成为可能"——系统梳理从采样到下游的标准化，指出不同DNA提取和文库制备方法对推断组成的影响可能超过生物效应本身。
- **mock community与标准化**：Tourlousse等2021（Microbiome 9:95, 10.1186/s40168-021-01048-3），由产业主导的日本微生物组联盟（JMBC）推动，用新开发的mock community对九种DNA提取协议（7种商业试剂盒、一种自建酚/氯仿法、以及IHMS国际人类微生物组标准推荐的QIAamp DNA Stool Mini Kit）共21种条件做head-to-head比较，通过MOSAIC Standards Challenge验证；原文呼吁"需建立常规监测分析性能和测试新方法的指南，包括可达性能的目标值，以确保跨方法和实验室/研究的可重复性与可比性"。
- **协变量平衡分配设计**（bioRxiv, 2025年3月24日, 10.1101/2025.03.21.644523）：在样本分配到批次时做协变量平衡，遵循通常归于George Box的箴言"能blocking就blocking，不能就randomize"（Box, Hunter and Hunter 2005）。该工具用模拟退火优化"平衡分数"（Kruskal-Wallis/Fisher精确检验p值的调和平均），发布为R包，并论证单纯随机化"是较差的策略"，因为它消除统计偏倚但不消除统计变异，可能残留批次变量与目标协变量间的混杂——从源头减少批次-生物混杂，是校正无法替代的一步。

## Recommendations

**阶段0（实验设计，最高优先级）**：用协变量平衡分配工具（如上述模拟退火R包）做批次分配，而非单纯随机化；每批次纳入mock community/spike-in（用于McLaren/metacal校准）与阴性对照（用于RUV类）。若能做到，很多下游校正问题可避免。**基准转变点**：若批次与关键表型的关联PERMANOVA R²可通过设计压到接近0，则无需激进校正。

**阶段1（诊断）**：先用mBatchNet或MBECS生成批次-表型mosaic图与PERMANOVA/silhouette诊断。**若批次与表型高度混杂**（一个表型组只出现在部分批次）→ 不要强行校正，转meta分析（Percentile Normalization或MMUPHin_MetaDA/BDMMA把批次作协变量）。

**阶段2（按目标选方法）**：
- 目标=**跨研究/跨批次预测建模**：首选DEBIAS-M（跨研究auROC最优、可解释）；MetaDICT在存在未观测混杂/完全混杂时更稳。
- 目标=**可视化/beta多样性/任意下游**且设计较平衡：ConQuR（输出计数、控FDR好）；非平衡batch×treatment设计用weighted/sparse PLSDA-batch。
- 目标=**病例-对照meta分析**：Percentile Normalization（简单稳健）或MMUPHin。
- 有**mock/spike-in或阴性对照**：metacal（McLaren模型）校准 / RUV-III-NB。
- 目标=**差异丰度并要诚实推断**：把批次作协变量建模（MaAsLin2/BDMMA），而非先校正再检验。

**阶段3（验证）**：多指标交叉验证——同时报告Bray-Curtis与Aitchison下的PERMANOVA R²、silhouette、以及生物信号保留（阳性对照taxa是否留存）。**基准转变点**：若某方法在Aitchison下R²升高或阳性对照信号丢失，说明过校正，应回退到更保守方法或改为协变量建模。

## Caveats
- **无金标准**：真实数据无"无批次ground truth"，所有基准都依赖代理指标；本报告的方法排名随指标/数据集变化，不应外推为普适结论。
- **性能数字的可比性有限**：不同论文用不同数据集、距离、指标，DEBIAS-M的auROC与ConQuR的PERMANOVA R²不能直接比较；引用的具体数字来自各方法的原文，存在作者自评偏向。
- **深度学习迁移仍属早期**：单细胞方法（Harmony/scVI/scANVI/autoencoder）在微生物组上多为潜力/探索性，缺乏大规模独立基准；微生物组样本量小、零膨胀使直接迁移受限。
- **mBatchNet（2026）与ISME J（2026）为最新发表**，部分结论尚待更广泛复现。
- 报告聚焦分类丰度层面；功能通路、病毒组/真菌组、绝对定量层面的批次校正证据相对薄弱，是明确空白。