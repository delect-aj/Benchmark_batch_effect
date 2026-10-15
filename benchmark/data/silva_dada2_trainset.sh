#!/bin/bash
# Turn the SILVA 138.2 release files into a DADA2 assignTaxonomy training set.
# taxmap gives every Bacteria/Archaea sequence a path that is exactly domain;phylum;class;order;family;genus;
# (checked: all 451,706 of them), which is the format assignTaxonomy wants, so no rank-designation file is
# needed and nothing has to be downloaded. The _DNA fasta is already T-based, so no U->T conversion either.
# Usage: data/silva_dada2_trainset.sh /beegfs/dongbiao/silva_138.2
set -eu
D="${1:-/beegfs/dongbiao/silva_138.2}"
OUT="$D/silva_138.2_dada2_genus.fa.gz"
awk -F'\t' '
  NR==FNR { if (FNR>1) { split($4,a,";"); if (a[1]=="Bacteria" || a[1]=="Archaea") tax[$1"."$2"."$3]=$4 } ; next }
  /^>/    { id=substr($0,2); sub(/ .*/,"",id); keep=(id in tax); if (keep) print ">" tax[id]; next }
  keep    { print }
' "$D/taxmap_slv_ssu_ref_nr_138.2.txt" "$D/SILVA_138.2_SSURef_NR99_DNA.fasta" | gzip -c > "$OUT.tmp"
mv "$OUT.tmp" "$OUT"
echo "built $OUT  ($(zcat "$OUT" | grep -c '^>') sequences, $(du -h "$OUT" | cut -f1))"
