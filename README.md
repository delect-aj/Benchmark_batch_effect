# Benchmark_batch_effect

微生物组批次校正方法的中立基准测试（16S + 宏基因组）。本文件是**操作手册**：远端服务器在哪、每个文件夹和脚本
是干什么的、整条流程怎么跑起来。

- 方法 wrapper 的接口契约和每个 wrapper 的状态 → [`benchmark/README.md`](benchmark/README.md)
- 每个数据集怎么来的 → [`benchmark/data/README.md`](benchmark/data/README.md)
- 完整设计与历史决策 → [`benchmark/PLAN.md`](benchmark/PLAN.md)
- 论文草稿 → [`benchmark/paper/`](benchmark/paper/)

---

## 1. 远端服务器

### 登录

先经跳板机进登录节点 `mu01`，再从 `mu01` 跳进计算节点（`cu01`~`cu10`）。具体地址、端口和口令不写在仓库里，
见本地的私人记录。

```bash
ssh -J <跳板> <登录节点>        # -> mu01
ssh cu03                        # mu01 上再跳计算节点，sinfo 看哪个空
```

家目录是共享的（NFS），所以在哪个节点看到的 `~/Benchmark_batch_effect` 都是同一份。

**算力任务必须在计算节点跑，不要在 mu01 上跑。** `sinfo` 看空闲节点；上一次完整跑用的是 cu03，20 核。

### 项目位置

```
~/Benchmark_batch_effect/        # git 仓库（= 本地这个仓库）
├── README.md                    # 本文件
├── benchmark/                   # 全部代码 + 数据 + 结果
└── script/                      # MicroCoRE（你自己的方法，不参与本 benchmark）
```

### 运行环境

| 环境 | 路径 | 用途 |
|---|---|---|
| `bench-r` | `~/software/miniconda3/envs/bench-r/bin/Rscript` | **主力**：绝大多数方法、全部指标、全部出图 |
| `bench-py` | `~/software/miniconda3/envs/bench-py/bin/python` | Python 方法（DEBIAS-M、scANVI）+ **snakemake 本体** |
| `bench-data` | `~/software/miniconda3/envs/bench-data/bin/Rscript` | 只用于 curatedMetagenomicData 取数 |
| R 4.1.3 | `~/software/R-4.1.3/bin/R` | 只用于 MBQC 的 DADA2（新版 Rscript 编译进了旧前缀，必须用这个） |

路径写死在 [`benchmark/config/config.yaml`](benchmark/config/config.yaml) 顶部，**不需要 `conda activate`**，
Snakefile 直接用绝对路径调用。重建环境：`bash benchmark/envs/build.sh <放 GitHub clone 的目录>`
（GitHub 包要在登录节点 clone，计算节点出不去外网）。

外部大文件：SILVA 训练集在 `/beegfs/dongbiao/silva_138.2/silva_138.2_dada2_genus.fa.gz`。

---

## 2. 目录地图

### 代码（都在 `benchmark/` 下，跟着 git 走）

| 目录 | 内容 |
|---|---|
| `Snakefile` | **总调度**。所有 rule 都在这一个文件里 |
| `config/` | 配置三件套 |
| `data/` | 把各个公开数据集做成统一格式的脚本 |
| `simulate/` | 模拟数据生成 |
| `methods/` | 每个批次校正方法一个 wrapper |
| `metrics/` | 打分脚本 |
| `report/` | 汇总脚本 |
| `figures/` | 五张主图：`.R` 脚本出成图，`.ipynb` 按子图逐个跑 |
| `pilot/` | 冒烟测试（跑通性检查） |
| `envs/` | 环境构建 |
| `paper/` | 论文草稿和图规划 |

### 结果（`benchmark/results/`，**全部 gitignore，只有 `results/data/` 入库**）

| 路径 | 大小 | 是什么 |
|---|---|---|
| `results/data/` | 621 M | **入库**。整理好的数据集（counts/meta/taxonomy/truth），重建太贵所以进 git |
| `results/raw/` | 113 M | 原始下载件（hivrc 94M、seed 11M、mbqc、jmbc、salim） |
| `results/out/` | **12 G** | 每个「方法 × 数据集」的校正结果 `.tsv` + `.status` + `.score` |
| `results/bench/` | 88 M | snakemake 记的耗时和峰值内存 |
| `results/logs/` | 203 M | 每个任务的 stdout/stderr |
| `results/loso/` | 480 K | 跨研究预测，每个数据集每个方法一个文件 |
| `results/truth/` | 444 K | 对已知配方的打分 |
| `results/figures/` | 5.9 M | 五张图的 PDF + PNG |
| `scores.tsv` | 7 M | **汇总表**，20,840 行，图主要读它 |
| `loso.tsv` / `truth_summary.tsv` / `truth_taxon.tsv` | | 另外三张汇总表 |
| `.snakemake/` | 235 M | snakemake 的元数据，**别删**（见第 6 节） |

> 服务器上 `results/` 里还有几个历史遗留：`pilot/`、`test_sim/`、`conqur_diag.txt`、
> `design_check*.txt`。生成它们的两个一次性脚本已经删掉了，这些文件可以一起删。

---

## 3. 脚本清单

### `config/` — 改这里就能改跑什么

| 文件 | 作用 |
|---|---|
| `config.yaml` | 解释器路径、跑哪些数据集、跑哪些方法、重复数、超时、每方法线程数 |
| `methods.tsv` | 方法清单：显示名、track（A=校正表 / B=嵌入）、输出类型、需不需要表型、角色。**图的坐标轴顺序和标签都来自这里** |
| `scenarios.tsv` | 213 个模拟场景的参数表。**已冻结**，有结果之后不要重新生成 |

### `data/` — 数据准备（每个脚本对应一批数据集）

| 脚本 | 产出 |
|---|---|
| `cmd.R` | `crc_mgx` `ibd_mgx` `ici_mgx`（curatedMetagenomicData，要 bench-data 环境） |
| `microbiomehd.R` | `crc_16s` `ibd_16s` `hiv_16s`（MicrobiomeHD，归并到属） |
| `hivrc_16s.R` | `hivrc_16s`（HIV 再分析联盟，源文件要手动从 Synapse 下） |
| `seed_16s.R` | `seed_16s`（种子微生物组） |
| `brooks_mock.R` | `brooks` 模拟群落（metacal 自带） |
| `salim_pig.R` | `salim_pig` 猪粪储存实验 |
| `mbqc.R` + `mbqc_download.sh` + `mbqc_manifest.sh` | `mbqc`。**最重的一步**：13 个实验室分开跑 DADA2，13 G 下载 + 十几小时 |
| `jmbc.R` + `jmbc_download.sh` / `jmbc_ref.sh` / `jmbc_quantify.sh` / `jmbc_tables.py` | `jmbc_extraction` + `jmbc_library`。从原始 reads 重建：建 20 株参考 → bowtie2 比对 → 计数 |
| `silva_dada2_trainset.sh` | 把 SILVA 138.2 做成 DADA2 训练集 |
| `inventory.sh` / `DATASETS.tsv` | 数据集索引 |

### `simulate/` — 模拟

| 脚本 | 作用 |
|---|---|
| `design.R` | 生成 `scenarios.tsv`：10 个 sweep（混杂梯度）+ 3 个 null（无信号）+ 200 个 LHS（拉丁超立方） |
| `make_templates.R` | 从真实数据抽 16S 和宏基因组两个模板 |
| `sim.R` | **核心**。MIDASim 拟合模板 → 种 DA taxa → 叠五种批次效应 → 抽 counts，同时输出无批次的 `oracle.tsv` |
| `toy_sim.R` | 冒烟测试用的小玩具数据 |

### `methods/` — 被测方法，接口完全一致

```
Rscript methods/<name>.R  counts.tsv meta.tsv out.tsv
python  methods/<name>.py counts.tsv meta.tsv out.tsv
```

Track A（输出校正后的表）：`limma` `combat` `combatseq` `mmuphin` `conqur` `conqur_default` `cqr`
`plsdabatch` `wplsdabatch` `percentile` `metadict` `debias_m`
Track B（输出嵌入）：`harmony` `fastmnn` `scanvi`
基线：`raw`（不校正）`oracle`（无批次真值，仅模拟数据）
对照：`conqur_permcov`（打乱表型的 ConQuR，用来隔离「条件于表型」带来的假阳性）
留着但跑不了：`ruviiinb`（当前 Bioconductor 下包内部报错，原因记在 benchmark/README.md）

公共函数在 `_common.R` / `_common.py`。详细契约见 [`benchmark/README.md`](benchmark/README.md)。

### `metrics/` — 打分

| 脚本 | 算什么 |
|---|---|
| `compute.R` | 主指标：三种距离下的批次/表型 PERMANOVA R²、ASW、iLISI/cLISI、对 oracle 的 Mantel r 和逐 taxon rho、DA 的 FDR/power/AP（对两套真值） |
| `truth.R` | 对**实际配制的配方**打分：Aitchison/Bray 距离、逐 taxon log2 偏差、同材料间的假结构 |
| `loso.R` | 留一研究交叉预测 auROC（L1 logistic + 随机森林）。**是上界不是估计**，校正时见过被留出研究的标签 |

### 其他

- `report/aggregate.R` — 把几万个 `.score` 拼成 `scores.tsv`
- `figures/fig1.R` … `fig5.R` — 五张主图，`_panels.R` 是共享的配色/排序/主题
- `figures/fig1.ipynb` … `fig5.ipynb` — 同样五张图的 notebook 版：**每个子图一个单元格，可以单独跑单独看**，
  最后一个单元格把它们拼成投稿用的组图并写 PDF。notebook 只 `source()` 对应的 `.R`，面板代码仍然只有一份，
  所以在 notebook 里调好的效果和 `snakemake` 出的图是同一个。
  `.R` 脚本末尾的 `main()` 加了守卫（`length(commandArgs(TRUE)) >= 2`），source 时不会触发出图
- `pilot/run_pilot.sh` + `pilot/check.R` — 冒烟测试：造玩具数据，把每个 wrapper 跑一遍，报 ok / FAIL

---

## 4. 运行逻辑

一句话：**准备数据 → 每个方法在每个数据集上跑一遍 → 每个结果打一次分 → 汇总成四张表 → 画五张图。**

```
         data/*.R                    simulate/sim.R
    （真实 + mock 数据集）          （213 场景 × 重复）
              │                              │
              └──────────┬───────────────────┘
                         ▼
              results/data/<ds>/counts.tsv + meta.tsv
                         │
                         ▼  rule run   ← methods/<m>.{R,py}
              results/out/<ds>/<m>.tsv + .status
                         │
        ┌────────────────┼────────────────┐
        ▼                ▼                ▼
   rule score        rule loso        rule truth
  (metrics/        (metrics/        (metrics/
   compute.R)        loso.R)          truth.R)
        │                │                │
        ▼                ▼                ▼
   .score 文件      loso/*.tsv      truth/*.tsv
        │                │                │
        ▼                └────────┬───────┘
  rule aggregate                  ▼
  (report/aggregate.R)       rule collect
        │                         │
        └────────► scores.tsv  loso.tsv  truth_summary.tsv  truth_taxon.tsv
                                  │
                                  ▼  rule figures
                           results/figures/fig1-5.pdf
```

几个关键设计：

- **一个任务一个文件。** 每个「方法 × 数据集」是独立的 job，失败只影响自己，整条流程可以断点续跑。
- **失败是结果，不是异常。** wrapper 非零退出 → `.status` 记 `fail`，`.score` 写一行空记录。方法的脆弱性是要被测量的，不是被绕过的。超时 6 小时同样记为失败。
- **打分不重跑校正。** `score` / `loso` / `truth` 都只读 `run` 已经产出的表，所以换指标不用重算方法。
- **哪个方法跑哪些数据集是有规则的**（Snakefile 里 `methods_for()`）：`oracle` 只在模拟数据上有意义；`conqur`（调参版）太贵，不上 200 个 LHS 场景，只跑真实数据和 13 个梯度场景；`conqur_permcov` 只在 sweep/null 场景上跑。
- **重复数分层**：sweep/null 场景 20 次重复（Fig 4 的曲线要精度），LHS 场景 5 次（要覆盖面）。

规模参考：完整一轮 **33,250 个 job，20 核在 cu03 上跑了 1 天 9 小时 43 分**，产出 20,840 条打分记录
（19,465 成功，1,375 失败）。

---

## 5. 常用操作

全部在 `cd ~/Benchmark_batch_effect/benchmark` 下执行。

```bash
# 完整跑 / 断点续跑（后台，断开 ssh 也不会停）
nohup bash ~/full_run.sh > ~/full.log 2>&1 &
tail -f ~/full.log

# 看进度
~/software/miniconda3/envs/bench-py/bin/snakemake --cores 20 --rerun-triggers mtime -n | tail -20

# 停掉
bash ~/stop_run.sh

# 开 jupyter 调图（kernel 选 "R (bench-r)"）
~/software/miniconda3/envs/metagenome_humann3/bin/jupyter lab --no-browser --port 8888
# 本地再开一条隧道转发 8888，浏览器打开给出的 token 链接

# 只重画图（最常用，几分钟）
for i in 1 2 3 4 5; do
  ~/software/miniconda3/envs/bench-r/bin/Rscript figures/fig$i.R results results/figures/fig$i.pdf
done

# 只准备数据，不跑方法
~/software/miniconda3/envs/bench-py/bin/snakemake --cores 10 data

# 小范围试跑（一个方法 + 一类场景）
~/software/miniconda3/envs/bench-py/bin/snakemake --cores 4 \
  --config methods_only='limma|raw|oracle' only=sweep_16s

# 冒烟测试：每个 wrapper 还能不能跑
RSCRIPT=~/software/miniconda3/envs/bench-r/bin/Rscript \
PYTHON=~/software/miniconda3/envs/bench-py/bin/python bash pilot/run_pilot.sh

# snakemake 被强杀后解锁
~/software/miniconda3/envs/bench-py/bin/snakemake --unlock
```

服务器家目录下有两个现成的脚本：

- `~/full_run.sh` — 三行：`cd` 到 benchmark，然后
  `snakemake --cores 20 --keep-going --rerun-triggers mtime`。`--keep-going` 是必须的：
  单个方法失败不应该停掉整条流程。
- `~/stop_run.sh` — 先收集 PID 再按 PID 杀。不能直接 `pkill -f`：那条命令自己的命令行里带着匹配模式，
  会把自己所在的 shell 一起杀掉。

---

## 6. 必须知道的坑

**① `--rerun-triggers mtime` 不能省。** snakemake 默认还会把「rule 的 shell 命令变了」当成重跑理由。
改过 Snakefile 之后不加这个参数，它会把已经跑完的几万个 job 全部重做。加了就只看文件时间。
反过来，当你**确实想让它重新展开 DAG**（比如新增了数据集或方法）时要去掉它，否则看不到待办任务。

**② 不要删 `.snakemake/`（235 M）。** 里面是 provenance 元数据。删了之后 snakemake 认不出哪些 job 是最新的，
默认触发器会要求重跑全部 33,250 个 job。真要删，之后每次都必须带 `--rerun-triggers mtime`。

**③ BLAS 线程要锁。** conda 的 R 链的是 OpenBLAS，不设 `OPENBLAS_NUM_THREADS` 的话，
一个申请 1 核的 job 会吃掉整个节点。Snakefile 的 `run`/`score`/`truth`/`loso` 四条 rule 里都已经 export 了。

**④ `scenarios.tsv` 冻结。** 已经有结果之后重新跑 `simulate/design.R` 会改掉场景定义，之前的结果全部作废。

**⑤ `rule template_mgx` 用了 `ancient()`。** 模拟模板依赖 `crc_mgx`，但模拟设计已冻结，
所以哪怕 crc_mgx 被重建，也不能让它触发重新拟合模板、重抽所有宏基因组场景。

**⑥ 被 `kill -9` 打断后可能留下孤儿 `.score`。** 如果 `aggregate` 在 100% 处报错，
通常是有 `.score` 文件存在但对应的 `.status` 被删了。删掉这些孤儿 `.score`，让 snakemake 重建这一对。

**⑦ 大跑之前先 `sinfo`。** 挑 idle 的节点，别和别人抢。

**⑧ R 的 jupyter kernel 是手工注册的。** kernelspec 在 `~/.local/share/jupyter/kernels/ir/kernel.json`，
直接指向 `bench-r` 的 R。`bench-py` 里的 jupyter 是残的（没有 `kernelspec` / `nbconvert` 子命令），
所以开 notebook 用 `metagenome_humann3` 那个 jupyter，kernel 选 "R (bench-r)"。
IRkernel 是用 `conda install --freeze-installed` 装进 bench-r 的 —— 这个参数保证只新增包、不动已有的，
benchmark 环境没被改过。

**⑨ 图里不要用非 ASCII 字符。** 默认 pdf 设备的字体没有 `→` 和 `é` 这些字形，会渲染成 `...` 和 `Cram..r`。
要用的话得把输出设备换成 `cairo_pdf`。

---

## 7. 当前状态

- 完整一轮已经跑完，四张汇总表和五张图都在 `results/` 里
- 五张图已调到投稿质量：180 mm 宽，225–240 mm 高，单页，矢量 PDF
- 论文的 Introduction / Methods 草稿在 `benchmark/paper/`，**参考文献里有 4 条标了 ⚠️ 需要核对**
- 还没做：Track C（MaAsLin2 / ANCOM-BC2 / SVA / BDMMA，输出的是差异丰度结果而非表，需要另一套契约）、
  metacal wrapper、kBET 和 PCR-batch 指标、inductive 版本的 LOSO（目前只有 transductive 上界）
- 5 个 ComBat-seq 任务撞上 6 小时超时（exit 124），记为失败，没有深究
