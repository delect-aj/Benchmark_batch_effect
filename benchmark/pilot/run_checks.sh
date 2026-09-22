#!/bin/bash
# Pre-run checks A2/A3 (design) and A1 (ConQuR). Run from benchmark/ on the compute node.
E=~/software/miniconda3/envs; SMK=$E/bench-py/bin/snakemake; R=$E/bench-r/bin/Rscript
set -u
echo "== A2/A3: simulate all scenarios, score raw + oracle only"
$SMK --cores 24 --keep-going --config reps=1 'methods=["raw","oracle"]' 'real=[]' 2>&1 | grep -E "Error|steps \(|Nothing" | tail -3
$R pilot/design_check.R results/scores.tsv > results/design_check.txt 2>&1; echo DESIGN_DONE

echo "== A1: null scenario x 10 reps, raw/oracle/conqur via workflow"
$SMK --cores 24 --keep-going --config reps=10 'only=^null_16s_conf0$' 'methods=["raw","oracle","conqur"]' 'real=[]' 2>&1 | grep -E "Error|steps \(|Nothing" | tail -3
run() {  # <datadir> <variant>
  o=results/out/$1; mkdir -p $o results/logs/$1
  timeout 6h $R pilot/$2.R results/data/$1/counts.tsv results/data/$1/meta.tsv $o/$2.tsv > results/logs/$1/$2.log 2>&1 &&
    $R metrics/compute.R results/data/$1 $o/$2.tsv $o/$2.score || echo "FAIL $1 $2"
}
export -f run; export R
for r in $(seq 1 10); do for v in conqur_nocov conqur_tuned; do echo "sim/null_16s_conf0/rep$r $v"; done; done |
  xargs -P 10 -n 2 bash -c 'run "$0" "$1"'
run real/crc_mgx conqur_nocov & run real/crc_mgx conqur_tuned & wait
$R pilot/conqur_diag.R results/out > results/conqur_diag.txt 2>&1
echo CHECKS_END
