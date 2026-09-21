#!/usr/bin/env bash
# Pilot: toy data -> every wrapper -> status table + contract check.
# Usage: bash pilot/run_pilot.sh [confounding 0..1]   (run from benchmark/)
set -u
cd "$(dirname "$0")/.."
data=results/pilot/data; out=results/pilot/out; logs=results/pilot/logs
mkdir -p "$out" "$logs"
Rscript simulate/toy_sim.R "$data" "${1:-0}" || exit 1

for w in methods/[!_]*.R methods/[!_]*.py; do
  m=$(basename "${w%.*}")
  case $w in *.R) run=Rscript ;; *) run=python3 ;; esac
  if $run "$w" "$data/counts.tsv" "$data/meta.tsv" "$out/$m.tsv" >"$logs/$m.log" 2>&1; then
    echo "ok    $m"
  else
    rm -f "$out/$m.tsv" "$out/$m.tsv.json"
    echo "FAIL  $m   ($(grep -m1 -iE 'error|no module|there is no package' "$logs/$m.log" | cut -c1-90))"
  fi
done

Rscript pilot/check.R "$data" "$out"
