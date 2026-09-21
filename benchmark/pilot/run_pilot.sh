#!/usr/bin/env bash
# Pilot: toy data -> every wrapper -> metrics/compute.R -> boundary check.
# Usage: bash pilot/run_pilot.sh [confounding 0..1]   (run from benchmark/)
set -u
cd "$(dirname "$0")/.."
data=results/pilot/data; out=results/pilot/out; logs=results/pilot/logs
mkdir -p "$out" "$logs"
${RSCRIPT:-Rscript} simulate/toy_sim.R "$data" "${1:-0}" || exit 1

for w in methods/[!_]*.R methods/[!_]*.py; do
  m=$(basename "${w%.*}")
  case $w in *.R) run=${RSCRIPT:-Rscript} ;; *) run=${PYTHON:-python3} ;; esac
  if $run "$w" "$data/counts.tsv" "$data/meta.tsv" "$out/$m.tsv" >"$logs/$m.log" 2>&1; then
    echo "ok    $m"
  else
    rm -f "$out/$m.tsv" "$out/$m.tsv.json"
    echo "FAIL  $m   ($(grep -m1 -iE 'error|no module|there is no package' "$logs/$m.log" | cut -c1-90))"
  fi
done

# Oracle (batch-free truth) is scored like a method: it bounds what correction can achieve
cp "$data/oracle.tsv" "$out/oracle.tsv" && echo '{"kind": "counts"}' > "$out/oracle.tsv.json"
scores=results/pilot/scores.tsv; rm -f "$scores"
for f in "$out"/*.tsv; do
  m=$(basename "${f%.tsv}")
  if ${RSCRIPT:-Rscript} metrics/compute.R "$data" "$f" "$out/$m.score" 2>"$logs/$m.score.log"; then
    [ -f "$scores" ] && tail -n +2 "$out/$m.score" >> "$scores" || cp "$out/$m.score" "$scores"
  else
    echo "SCORE-FAIL  $m   ($(grep -m1 -i error "$logs/$m.score.log" | cut -c1-90))"
  fi
done
${RSCRIPT:-Rscript} pilot/check.R "$scores"
