#!/bin/bash
# Build a bowtie2 index of the 20 JMBC mock strains from the assembly accessions in Table S1.
# We know exactly what is in this material, so a ~75 Mb reference does the job and no MetaPhlAn or GTDB
# database has to be downloaded. ENA serves whole assemblies (chromosome plus plasmids) by GCA accession.
# Usage (run from benchmark/): data/jmbc_ref.sh results/raw/jmbc
set -eu
D="${1:-results/raw/jmbc}"
BIN="${BIN:-$HOME/software/miniconda3/envs/metagenome_amr/bin}"
MOCK="data/jmbc/jmbc_mock.tsv"

mkdir -p "$D/ref"
: > "$D/ref/all.fa"
while IFS=$'\t' read -r species strain assembly rest; do
  f="$D/ref/$assembly.fa"
  [ -s "$f" ] || curl -sS --retry 5 --connect-timeout 30 -o "$f" \
      "https://www.ebi.ac.uk/ena/browser/api/fasta/$assembly?download=true"
  # rename every contig to its assembly so reads can be tallied per strain
  awk -v a="$assembly" '/^>/ {print ">" a "_" (++n); next} {print}' "$f" >> "$D/ref/all.fa"
done < <(tail -n +2 "$MOCK")

echo "reference: $(grep -c '^>' "$D/ref/all.fa") contigs from $(( $(wc -l < "$MOCK") - 1 )) assemblies, \
$(du -h "$D/ref/all.fa" | cut -f1)"
"$BIN/bowtie2-build" --threads 4 -q "$D/ref/all.fa" "$D/ref/mock"
echo "index built: $D/ref/mock"
