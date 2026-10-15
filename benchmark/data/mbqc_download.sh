#!/bin/bash
# Fetch the forward reads of every MBQC run (SRP047083) from ENA into <raw>/fastq.
# Resumable: a file whose size already matches the manifest is skipped, a partial one is continued.
# Four streams, because this waits on the network rather than the CPU (~13 GB, a few hours on this cluster).
# Usage: data/mbqc_download.sh results/raw/mbqc
set -u
D="${1:-results/raw/mbqc}"
[ -f "$D/mbqc_manifest.tsv" ] || { echo "no manifest at $D/mbqc_manifest.tsv" >&2; exit 1; }
mkdir -p "$D/fastq"; cd "$D/fastq" || exit 1

get_one() {
  run=$1; url=$2; want=$3; f="$run.fastq.gz"
  have=$(stat -c %s "$f" 2>/dev/null || echo 0)
  [ "$have" = "$want" ] && return 0
  curl -sS -L -C - --retry 5 --retry-delay 5 --connect-timeout 30 -o "$f" "$url" || return 1
  have=$(stat -c %s "$f" 2>/dev/null || echo 0)
  [ "$have" = "$want" ] || { echo "SIZE MISMATCH $run: got $have want $want" >&2; return 1; }
}
export -f get_one

tail -n +2 "../mbqc_manifest.tsv" | awk -F'\t' '{print $1"\t"$4"\t"$5}' \
 | xargs -P 4 -n 1 -d '\n' bash -c 'IFS=$'"'"'\t'"'"' read -r r u s <<< "$0"; get_one "$r" "$u" "$s"'

n=$(ls -1 *.fastq.gz 2>/dev/null | wc -l)
want=$(( $(wc -l < ../mbqc_manifest.tsv) - 1 ))
echo "=== $n / $want files, $(du -sh . | cut -f1) ==="
[ "$n" = "$want" ] && touch .complete && echo "marked complete"
