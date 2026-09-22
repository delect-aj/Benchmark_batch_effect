# 微生物组批次校正工具中立Benchmark方案

## Context
文献调研（`文献调研_batch_effect_microbiome.md`）指出三大空白：①无金标准评价体系，排名随距离/指标翻转（Park 2025 BC vs Aitchison）；②"移除批次"(PERMANOVA R²)与"跨研究泛化"(LOSO AUC)结论矛盾（ConQuR原文 vs DEBIAS-M）；③混杂设计下所有方法失效，但缺乏系统量化。现有基准多为方法作者自评、单数据集、单指标。
目标：做一个**中立、多维、带ground truth、覆盖16S+宏基因组、系统梯度化混杂程度**的benchmark，产出"按研究目标×实验设计"的选型指南。MicroCoRE（`script/`）不纳入（用户选择纯中立）。

## 1. 纳入方法（按输出类型分轨，避免不公平比较）

| 轨道 | 方法 |
|---|---|
| A. 输出校正丰度表 | ComBat、ComBat-seq、limma removeBatchEffect、MMUPHin、ConQuR（及ConQuR-libsize）、PLSDA-batch / weighted / sparse、Percentile Normalization、复合分位数回归(Park 2025)、DEBIAS-M、MetaDICT、RUV-III-NB（仅有重复/对照的数据集） |
| B. 输出嵌入 | Harmony、fastMNN、scANVI、MetaDICT嵌入 |
| C. 不校正、建模批次（DA参照） | MaAsLin2/ANCOM-BC2 + batch协变量、BDMMA、SVA代理变量、MMUPHin_MetaDA |
| D. 校准（仅mock数据） | metacal |
| 基线/边界 | Raw（下界）、Oracle（模拟中无批次的真值数据，上界）、仅CLR/TSS |

中立规则：每方法用作者推荐默认参数 + 作者文档推荐的调参方式各跑一次；实施前把wrapper发给原作者确认（邮件/GitHub issue）；预注册评价指标与汇总权重（OSF），数据跑完前冻结。

## 2. 数据

### 2.1 模拟（有完整ground truth，主力）
- 生成器：以真实数据为模板的 **MIDASim** / **SparseDOSSA2**（16S模板：HMP/AGP；宏基因组模板：curatedMetagenomicData健康肠道），保留真实稀疏度、相关结构。
- 生物信号：植入已知DA taxa（比例5%/10%，效应量log2FC 1/2/3），同时植入一个连续协变量。
- 批次效应模型（分别与组合）：
  1. 乘性per-taxon偏倚（McLaren模型，最贴近真实提取/引物偏倚）
  2. log尺度加性+尺度（ComBat假设，作为"对ComBat友好"对照）
  3. 批次特异dropout/零膨胀改变
  4. 非系统性（仅部分taxa受影响，比例10%–100%）
  5. 测序深度差异
- **混杂梯度**（核心创新点）：batch×表型关联 Cramér's V = 0 / 0.3 / 0.6 / 0.9 / 1.0（完全混杂）；批次数 2/5/10；每批样本量 20/50/200；批次大小不平衡。
- 零信号null场景：无生物信号 → 检测"方法制造假信号"（过校正/泄漏）。
- 规模：全因子网格用拉丁超立方抽样约200个参数组合 × 20重复 ≈ 4000数据集/方法，集群跑。

### 2.2 带实验性ground truth的真实数据
- Mock community / 技术重复：Tourlousse 2021 JMBC（9种提取协议）、MBQC、Brooks 2015 mock、Salim 2022猪粪便（184样本/21种技术组合）。真值=已知组成或同一样本的重复应一致。
- 用途：验证乘性偏倚模型是否成立，校准模拟参数（让模拟的批次强度落在真实范围）。

### 2.3 真实meta分析数据（无真值，用代理指标）
- 16S：HIV（17研究，DEBIAS-M用）、IBD/CRC 16S（MicrobiomeHD, Duvallet 2017）、种子微生物组（Foxx & Rivers 2025，天然混杂案例）
- 宏基因组：CRC 8队列（curatedMetagenomicData，Wirbel 2019）、IBD（MMUPHin用）、免疫治疗应答（MetaDICT用）
- 统一预处理：16S统一到genus层（DADA2/SILVA重注释或采用原作者表），宏基因组MetaPhlAn4 species层；统一过滤（prevalence ≥10%）。

## 3. 评价指标（四维 + 实用性）

1. **批次去除**：PERMANOVA R²(batch)在 Bray-Curtis / Aitchison / Jaccard / (16S有树时)UniFrac 下分别报告；iLISI、kBET、batch ASW、PCR-batch。
2. **生物信号保留**：PERMANOVA R²(表型)、cLISI、bio ASW；模拟中额外：与Oracle的距离（Aitchison距离 / per-taxon Spearman）、alpha多样性保真、taxon-taxon相关结构保真（与真值相关矩阵Mantel r）。
3. **下游推断**：DA的FDR、power、AUPRC（模拟真值）；真实数据用跨研究签名一致性（留一研究DA结果重叠）。null场景的I类错误率。
4. **跨研究预测**：留一研究(LOSO) auROC，RF + L1-logistic两个分类器。**严格防泄漏**：方法只在训练批次上fit再transform测试批次；无法做到的（转导式）单独标注并报告两种设定。
5. **实用性**：运行时间、峰值内存、失败率、所需输入（是否需表型标签/参考批次/对照）。

过校正专项：完全混杂场景中"批次R²下降"与"生物R²下降"的比值；表型标签打乱后校正是否仍"保留"表型信号（泄漏检测）。

## 4. 汇总与报告
- 每维度内对指标做min-max标准化后取均值（scIB风格），维度间**不给单一总分**，而是按目标给权重方案：可视化/探索、DA推断、预测建模三套。
- 主图：方法×场景热图；混杂梯度折线图（x=Cramér's V，y=各维度得分）；距离度量翻转的bump chart（直接回应文献矛盾）。
- 产出决策树：研究目标 × 混杂程度 × 数据类型 → 推荐方法。

## 5. 工程实现
```
benchmark/
  envs/              # 每方法一个conda env / Singularity镜像（R与Python依赖冲突大）
  methods/<name>.R|py  # 统一接口：in counts.tsv + meta.tsv(batch, phenotype, covars) → out corrected.tsv 或 embedding.tsv + runtime.json
  simulate/          # MIDASim/SparseDOSSA2 + 批次注入
  data/              # 真实数据下载与预处理脚本（不存原始数据进git）
  metrics/           # 统一评价，输入方法输出+真值
  Snakefile          # 数据集×方法×重复 调度，SLURM profile
  report/            # 汇总与作图
```
- 工作流：Snakemake + SLURM profile；每任务设超时（如6h）与内存上限，超时记为失败而非丢弃。
- 固定随机种子；记录每方法版本号（sessionInfo / pip freeze）。
- 代码与结果公开（GitHub + Zenodo），便于后续新方法接入（提供接口模板即可）。

## 6. 阶段与时间（参考）
1. 第1–2周：方法wrapper + 容器；在1个模拟、1个真实数据集上跑通全部方法（pilot）
2. 第3–4周：用mock/技术重复数据校准模拟器参数；冻结并预注册指标与权重
3. 第5–8周：全量模拟网格 + 真实数据
4. 第9–10周：汇总、作图、决策树、联系作者复核异常结果
5. 第11–12周：撰写论文（目标：Genome Biology / Microbiome / Briefings in Bioinformatics）

## 7. 验证（确保benchmark本身可信）
- **边界检查**：每个场景Oracle须在所有生物保留指标上最优、Raw在批次去除上最差；否则指标或模拟有bug。
- **已知结论复现**：在ComBat友好的加性场景中ComBat应表现良好；ConQuR在HIVRC上批次方差下降≈94%（复现原文）；DEBIAS-M在HIV LOSO auROC≈0.70。复现不了先查wrapper再联系作者。
- null场景中所有方法DA的I类错误率应可计算且Raw≈名义水平。
- 同一参数不同种子结果方差报告，确认20重复足够（重复数加倍看排名是否稳定，Kendall τ > 0.9）。

## 待定/风险
- 部分方法（Park 2025 CQR、MetaDICT）实现成熟度未知 → pilot阶段确认能否运行，不能则记为"不可用"并报告。
- BDMMA等贝叶斯方法在大网格上可能过慢 → 只在子网格上跑并注明。

## Pilot 记录（2026-09-22，cu01）
实施中与原方案不同之处（均已写入代码/README）：
- **模拟器**：仅用 MIDASim（parametric 模式）。SparseDOSSA2 不在 bioconda，且 CentOS 7 无法源码编译，放弃。
- **复合分位数回归 (Park 2025)**：官方代码无法运行（NB 步骤引用不存在的系数、分位数步骤依赖未定义全局变量），按论文部分重实现：robust-CV 选参考批次 + ConQuR composite 分位数回归，省略 NB 步骤。结果需标注为"部分重实现"。
- **RUV-III-NB**：在当前 Bioconductor 下两个入口（ruvIII.nb / fastruvIII.nb + get.res）均有包内部错误，记为"不可用"，保留 wrapper。
- **ANCOM-BC2（轨道 C）**：需 CVXR < 1.0，conda-forge 无此版本，推迟到轨道 C 启动时从 CRAN archive 安装。
- **curatedMetagenomicData**：与主 R 环境 TBB 版本冲突，单独放在 bench-data 环境，仅用于下载真实数据。

Pilot（1 次重复；16S 混杂梯度 5 个 + null 3 个 + CRC 真实数据）初步观察，**待多重复确认，不作结论**：
- 完全混杂（conf=1）时 ComBat、ComBat-seq、MMUPHin、MetaDICT、Percentile 直接拒绝运行，这本身是结果。
- null 场景中 limma 在 conf=1 时产生 140 个假 DA；ConQuR/CQR 在平衡 null（conf=0）也产生 15/26 个假 DA，需核查是方法特性还是 wrapper/检验问题。
- CRC 上 ConQuR/CQR 使批次 R² 高于未校正数据，需核查参考批次选择。
- 嵌入类方法（scANVI、fastMNN、Harmony）的 oracle Mantel 值偏低，Aitchison 距离对非线性嵌入可能不公平，指标需再议。

## 跑前检查结果（2026-09-22）
**A2 批次强度**：未校正 batch R²（Aitchison）中位数 0.029、90% 分位 0.096、最大 0.256；Bray-Curtis 中位 0.082、最大 0.435。真实数据 CRC 为 0.077（Aitchison）/ 0.095（BC），文献 HIVRC BC ≈ 0.12。设计整体偏弱，MGX 模板尤甚（bias_sd=2 时仅 0.029）。建议 bias_sd 由 {0.5,1,2} 改为 {1,2,3}（待确认）。

**A3 DA 效应**：不是天花板而是地板——无批次的 Oracle 在 57% 场景中 power < 0.2，总样本量 ≤150 时 power ≈ 0。建议预注册：DA 主指标用 AP；FDR/power 只在 Oracle power ≥ 0.2 的场景报告（待确认）。

**B1 DA 真值**：Oracle 在所有混杂水平下平均 FDR ≈ 0.15（名义 0.05），来自组成性重归一化使非 DA taxa 的 CLR 也移动。sim.R 已额外输出每个 taxon 的真实 log2FC（truth_fc.tsv，不改变模拟数据），指标可在跑完后修正。

**A1 ConQuR**（null 场景 ×10 重复 + CRC）：
| 设置 | null 假 DA（均值） | null batch R² | CRC batch R² BC / Aitchison / Jaccard |
|---|---|---|---|
| 未校正 | 0 | 0.074 | 0.095 / 0.077 / 0.103 |
| ConQuR 默认（参考批次=字母序第一，保护真实表型） | 9.7 | 0.026 | 0.008 / 0.182 / 0.503 |
| ConQuR，协变量=打乱的表型 | 0.0 | 0.026 | 0.008 / 0.183 / 0.505 |
| Tune_ConQuR（作者推荐调参） | 2.2 | 0.016 | 0.010 / 0.090 / 0.198 |
结论（方法特性，非 wrapper 错误）：
1. 假 DA 完全来自以真实表型为条件：ConQuR 把每个值替换为"该批次+该表型"的拟合分位数，噪声中的偶然表型差异被系统化，之后用同一表型检验是循环论证。
2. CRC 上 Aitchison/Jaccard 反升来自零值处理：零值对应一段分位数，ConQuR 取参考分位数均值并取整；非参考批次零比例远高于参考批次时几乎所有零都变成正数（HanniganGD 零比例 82% → 2.3%）。BC 被高丰度 taxa 主导看不出，零敏感距离暴露出大量人为批次结构——这正是文献中"不同距离下排名翻转"的一个具体机制。
3. 参考批次选择影响巨大；Tune_ConQuR 选 YuJ_2015 + lasso 后各批次零比例均衡（46–57%）。
Tune_ConQuR 依赖已被移除的 vegan::adonis，需兼容 shim（pilot/conqur_tuned.R）。

## 冻结前的决定（2026-09-22，用户确认）
1. bias_sd 取值改为 {1, 2, 3}（只改这一列，其余场景参数不变）。
2. DA 维度预注册：主指标为 AP；FDR 与 power 只在 Oracle power ≥ 0.2 的场景中报告。
3. ConQuR 主结果用 Tune_ConQuR（参考批次池 = 最大的 3 个批次）；默认 ConQuR 和"打乱表型协变量"版本作为敏感性分析，只跑 sweep/null 场景。
