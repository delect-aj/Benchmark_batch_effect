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
