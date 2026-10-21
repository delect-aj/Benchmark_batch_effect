#!/bin/bash
# Regenerate data/DATASETS.tsv: one row per prepared dataset, measured from the files themselves and joined
# with the provenance below. Run it from benchmark/ after any dataset is rebuilt, so the table cannot go stale.
#   data/inventory.sh
# Everything that can be derived (sample and feature counts, batch and phenotype levels, median depth, which
# files exist, which covariates are present, whether it is in the evaluation list) is measured here. Only what
# a file cannot tell you -- where the data came from, what batch and phenotype actually mean, and the caveats
# that matter when reading results -- is curated, and that lives in the block directly below.
set -eu
OUT=data/DATASETS.tsv

# dataset | assay | feature_level | batch_meaning | phenotype_meaning | source | build_script | rule | raw_input | notes
curated() {
cat <<'PROV'
crc_mgx	metagenome	species	研究队列	CRC vs 对照	curatedMetagenomicData 3.14（MetaPhlAn 3 统一处理）	data/cmd.R	cmd	-	相对丰度反推的估计计数，非真实测序计数；提取试剂盒完全嵌套在研究内
ibd_mgx	metagenome	species	研究队列	IBD vs 对照	curatedMetagenomicData 3.14	data/cmd.R	cmd	-	同上：估计计数
ici_mgx	metagenome	species	研究队列	免疫治疗应答 vs 无应答	curatedMetagenomicData 3.14	data/cmd.R	cmd	-	同上：估计计数；样本量最小
crc_16s	16S	genus	研究	CRC vs 健康	MicrobiomeHD（Zenodo 569601，100% de novo OTU + RDP）	data/microbiomehd.R	microbiomehd	results/raw/microbiomehd	各研究扩增区域与平台不同；OTU 不可跨研究比较，故聚合到属
ibd_16s	16S	genus	研究	CD/UC vs 健康/nonIBD	MicrobiomeHD（Zenodo 569601）	data/microbiomehd.R	microbiomehd	results/raw/microbiomehd	同上
hiv_16s	16S	genus	研究	HIV vs 健康	MicrobiomeHD（Zenodo 569601）	data/microbiomehd.R	microbiomehd	results/raw/microbiomehd	同上；含 454 与 MiSeq 混合
hivrc_16s	16S	genus	研究	HIV 阳性 vs 阴性	HIVRC / Synapse syn18406854（Resphera Insight 统一注释，用户手动下载）	data/hivrc_16s.R	hivrc_16s	results/raw/hivrc	已抽平到每样本 20000 reads，深度差异被移除；2 个研究无阴性对照（天然混杂），Percentile Normalization 会失败
seed_16s	16S	ASV	研究	宿主植物（多分类）	Foxx & Rivers 2025（Mendeley 10.17632/5xrfg5dym6.1，QIIME2 ASV）	data/seed_16s.R	seed_16s	results/raw/seed	天然混杂案例；ASV 为序列哈希无分类注释；流行度过滤改为"至少在一个研究中 ≥10%"
mbqc	16S	genus	处理实验室（湿实验）	标本编号	MBQC SRP047083 原始 reads 自行重建（整合 OTU 表已下线）	data/mbqc.R	mbqc	results/raw/mbqc	完全交叉、零混杂（Cramér V = 0），混杂梯度的零点参照；10 万组同标本重复；按实验室分别去噪后聚合到属；6 个源头损坏的 fastq 已跳过
jmbc_extraction	metagenome	strain	DNA 提取方案（9 种）	cell mock vs 粪便	Tourlousse 2021 Microbiome 9:95 原始 reads 自行重建	data/jmbc.R	jmbc_extraction	results/raw/jmbc	全基准唯一以提取方案为批次；每 run 只取前 20 万条 read；MAPQ≥10 以排除两个 B. longum 亚种的歧义比对
brooks	16S	species	测序板	混入菌种数	Brooks 2015，metacal 包自带（细胞混合子集）	data/brooks_mock.R	brooks_mock	-	仅 2 块板，批次数少；仅用于校准乘性偏倚模型
salim_pig	metagenome	genus	储存条件（19 种）	基质（猪粪×2/污水×2）	Poulsen 2021 Microbiol Spectr / figshare 12010989（MGmapper 计数表）	data/salim_pig.R	salim_pig	results/raw/salim	125 组技术重复；另有 library_prep / lane / freeze_thaw 可作批次；MGmapper 参考库计数，与 MetaPhlAn 不可直接比较
jmbc_library	metagenome	strain	建库方案（28 种）	无（同一份 DNA mock）	Tourlousse 2021，同 jmbc_extraction	data/jmbc.R	jmbc_library	results/raw/jmbc	纯 null 场景：生物信号完全不变，任何残留结构都是批次、任何差异都是假信号
PROV
}

eval_list=$(grep '^real:' config/config.yaml | sed 's/.*\[//;s/\].*//;s/,/ /g')

{
printf 'dataset\ttrack\tevaluated\tpath\tassay\tfeature_level\tsamples\tfeatures\tbatch_var\tn_batch\t'
printf 'phenotype_var\tn_phenotype\tmedian_depth\thas_truth\thas_taxonomy\tcovariates\tsource\tbuild_script\t'
printf 'snakemake_rule\traw_input\tnotes\n'

for d in results/data/real/*/ results/data/mock/*/; do
  [ -s "$d/counts.tsv" ] || continue
  name=$(basename "$d"); track=$(echo "$d" | cut -d/ -f3)
  prov=$(curated | awk -F'\t' -v n="$name" '$1 == n')
  [ -n "$prov" ] || { echo "no provenance entry for $name" >&2; continue; }

  rows=$(( $(wc -l < "$d/counts.tsv") - 1 ))
  cols=$(( $(head -1 "$d/counts.tsv" | awk -F'\t' '{print NF}') - 1 ))
  hdr=$(head -1 "$d/meta.tsv")
  bi=$(echo "$hdr" | tr '\t' '\n' | grep -nx batch     | cut -d: -f1 || true)
  pi=$(echo "$hdr" | tr '\t' '\n' | grep -nx phenotype | cut -d: -f1 || true)
  nb=$(awk -F'\t' -v i="${bi:-0}" 'NR>1 && i>0 {print $i}' "$d/meta.tsv" | sort -u | wc -l)
  np=$(awk -F'\t' -v i="${pi:-0}" 'NR>1 && i>0 {print $i}' "$d/meta.tsv" | sort -u | wc -l)
  depth=$(awk -F'\t' 'NR>1 {s=0; for (i=2; i<=NF; i++) s+=$i; print s}' "$d/counts.tsv" \
          | sort -n | awk '{a[NR]=$1} END {print a[int((NR+1)/2)]}')
  cov=$(echo "$hdr" | cut -f2- | tr '\t' ',')
  ev=no; for e in $eval_list; do [ "$e" = "$name" ] && ev=yes; done
  [ -s "$d/truth.tsv" ]    && truth=yes    || truth=no
  [ -s "$d/taxonomy.tsv" ] && taxonomy=yes || taxonomy=no

  echo "$prov" | awk -F'\t' -v OFS='\t' -v track="$track" -v ev="$ev" -v p="${d%/}" -v r="$rows" -v c="$cols" \
      -v nb="$nb" -v np="$np" -v dep="$depth" -v t="$truth" -v tax="$taxonomy" -v cov="$cov" \
      '{print $1, track, ev, p, $2, $3, r, c, $4, nb, $5, np, dep, t, tax, cov, $6, $7, $8, $9, $10}'
done

# Simulation inputs: templates MIDASim is fitted to, and the scenario grid driving the simulated datasets
for t in 16s mgx; do
  f=data/templates/$t.tsv
  [ -s "$f" ] || continue
  printf '%s\tsim\tn/a\t%s\t%s\ttemplate\t%s\t%s\t-\t0\t-\t0\t-\tno\tno\t-\t%s\tsimulate/make_templates.R\ttemplate_%s\t-\t%s\n' \
    "template_$t" "$f" "$([ "$t" = 16s ] && echo 16S || echo metagenome)" \
    "$(( $(wc -l < "$f") - 1 ))" "$(( $(head -1 "$f" | awk -F'\t' '{print NF}') - 1 ))" \
    "$([ "$t" = 16s ] && echo 'HMP/AGP 健康肠道' || echo 'crc_mgx 健康样本')" "$t" \
    "MIDASim 参数拟合结果存为 $f.midasim.rds，sim 规则据此生成模拟数据"
done
s=config/scenarios.tsv
[ -s "$s" ] && printf 'scenarios\tsim\tn/a\t%s\t-\t-\t%s\t%s\t-\t0\t-\t0\t-\tno\tno\t-\t本项目设计\tsimulate/design.R\t-\t-\t%s\n' \
  "$s" "$(( $(wc -l < "$s") - 1 ))" "$(( $(head -1 "$s" | awk -F'\t' '{print NF}') - 1 ))" \
  "模拟场景网格；每个场景 × config 的 reps 次重复由 sim 规则生成到 results/data/sim/"
} > "$OUT"

echo "wrote $OUT: $(( $(wc -l < "$OUT") - 1 )) rows"
awk -F'\t' 'NR>1 {n[$2]++; s+=$7} END {for (k in n) printf "  %-5s %d datasets\n", k, n[k]; print "  samples total: " s}' "$OUT"
