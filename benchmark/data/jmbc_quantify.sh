#!/bin/bash
# Map each run against the 20 mock genomes and count uniquely-placed reads per strain.
# Two of the strains are Bifidobacterium longum subspecies, so reads that cannot be told apart would be
# assigned arbitrarily; MAPQ >= 10 drops them rather than letting bowtie2 break the tie at random.
# Usage (run from benchmark/): data/jmbc_quantify.sh results/raw/jmbc [threads]
set -u
D="${1:-results/raw/jmbc}"
P="${2:-4}"
BIN="${BIN:-$HOME/software/miniconda3/envs/metagenome_amr/bin}"
out="$D/counts_long.tsv"
: > "$out.part"
n=0
for f in "$D"/reads/*.fastq.gz; do
  run=$(basename "$f" .fastq.gz); n=$((n + 1))
  "$BIN/bowtie2" -x "$D/ref/mock" -U "$f" --no-unal -p "$P" 2>/dev/null \
    | "$BIN/samtools" view -q 10 -F 0x900 - \
    | awk -v r="$run" '{ g = $3; sub(/_[0-9]+$/, "", g); c[g]++ }
                       END { for (k in c) print r "\t" k "\t" c[k] }' >> "$out.part"
  [ $(( n % 20 )) -eq 0 ] && echo "$(date '+%T') mapped $n runs"
done
mv "$out.part" "$out"
echo "=== $n runs, $(wc -l < "$out") run x strain rows ==="
