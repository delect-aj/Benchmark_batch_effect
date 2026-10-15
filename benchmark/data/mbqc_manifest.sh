#!/bin/bash
# Build the two MBQC metadata tables from ENA, so the dataset is reproducible without the consortium's
# (now offline) integrated table:
#   mbqc_manifest.tsv  run -> specimen, handling lab, forward-read URL, size   [lab comes from library_name,
#                      which ENA stores as "<specimen>.<handling lab>.<id>"]
#   mbqc_samples.tsv   the 22 specimens with type / health / sex / age, and which two are artificial colonies
# Usage: data/mbqc_manifest.sh results/raw/mbqc
set -eu
D="${1:-results/raw/mbqc}"; mkdir -p "$D"; cd "$D"
API=https://www.ebi.ac.uk/ena/portal/api/filereport

curl -sS "$API?accession=SRP047083&result=read_run&fields=run_accession,library_name,fastq_ftp,fastq_bytes" \
 | awk -F'\t' 'BEGIN{OFS="\t"; print "run","specimen","lab","url","bytes"}
   NR>1 {
     n=split($3,f,";"); split($4,b,";"); split($2,p,".")
     url=""; sz=0
     for(i=1;i<=n;i++) if (f[i] ~ /_1\.fastq\.gz$/)                 { url=f[i]; sz=b[i]; break }
     if (url=="") for(i=1;i<=n;i++) if (f[i] ~ /\/[SE]RR[0-9]+\.fastq\.gz$/) { url=f[i]; sz=b[i]; break }
     if (url!="") print $1, p[1], p[2], "https://"url, sz
   }' > mbqc_manifest.tsv

curl -sS "$API?accession=SRP047083&result=read_run&fields=sample_accession,sample_alias" \
 | awk -F'\t' 'NR>1 && !seen[$2]++ {print $1"\t"$2}' | sort -k2,2 > .samples
mkdir -p .sxml
while IFS=$'\t' read -r acc alias; do
  [ -s ".sxml/$acc.xml" ] || curl -sS -m 90 --retry 3 -o ".sxml/$acc.xml" "https://www.ebi.ac.uk/ena/browser/api/xml/$acc"
done < .samples

python3 -I - <<'EOF' > mbqc_samples.tsv
import re, io
want = ["specimen type", "health_status_collapsed", "sex", "host age"]
print("\t".join(["specimen", "biosample", "specimen_type", "health_status", "sex", "age", "artificial"]))
for line in io.open(".samples", encoding="utf-8"):
    acc, alias = line.rstrip("\n").split("\t")
    a = dict(re.findall(r"<TAG>(.*?)</TAG>\s*<VALUE>(.*?)</VALUE>",
             io.open(f".sxml/{acc}.xml", encoding="utf-8", errors="replace").read(), re.S))
    v = [a.get(w, "NA").strip() for w in want]
    v = ["NA" if x in ("missing", "not applicable", "") else x for x in v]
    print("\t".join([alias, acc] + v + ["1" if "artificial" in v[0].lower() else "0"]))
EOF

echo "runs: $(( $(wc -l < mbqc_manifest.tsv) - 1 ))   specimens: $(( $(wc -l < mbqc_samples.tsv) - 1 ))"
awk -F'\t' 'NR>1{s+=$5} END{printf "forward reads to download: %.1f GB\n", s/1e9}' mbqc_manifest.tsv
