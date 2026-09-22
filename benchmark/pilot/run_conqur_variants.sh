#!/bin/bash
# Rerun only the two ConQuR diagnostic variants (A1 b/c) on existing data, then summarise.
E=~/software/miniconda3/envs; R=$E/bench-r/bin/Rscript; export R
run() {
  o=results/out/$1; mkdir -p $o results/logs/$1
  timeout 6h $R pilot/$2.R results/data/$1/counts.tsv results/data/$1/meta.tsv $o/$2.tsv > results/logs/$1/$2.log 2>&1 &&
    $R metrics/compute.R results/data/$1 $o/$2.tsv $o/$2.score || echo "FAIL $1 $2"
}
export -f run
{ for r in $(seq 1 10); do for v in conqur_permcov conqur_tuned; do echo "sim/null_16s_conf0/rep$r $v"; done; done
  echo "real/crc_mgx conqur_permcov"; echo "real/crc_mgx conqur_tuned"; } | xargs -P 12 -n 2 bash -c 'run "$0" "$1"'
$R pilot/conqur_diag.R results/out > results/conqur_diag.txt 2>&1
echo VARIANTS_END
