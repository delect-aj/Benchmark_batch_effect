#!/bin/bash
# Fetch only the head of each JMBC run. 200,000 reads is far more than enough to quantify a 20-strain mock,
# and ENA honours HTTP Range, so this costs a few GB instead of the 100 GB these files hold at full depth.
# The truncated gzip stream is decompressed as far as it goes and cut at the read count we want.
# Usage (run from benchmark/): data/jmbc_download.sh results/raw/jmbc [n_reads]
set -u
D="${1:-results/raw/jmbc}"
NREAD="${2:-200000}"
BYTES=$(( NREAD * 130 ))          # gzipped 150 bp reads run ~70 bytes each; ask for roughly double
MIN=$(( NREAD / 5 ))              # accept a run only if it yielded at least a fifth of the target
mkdir -p "$D/reads"

get() {
  run=$1; url=$2; out="$D/reads/$run.fastq.gz"
  [ -s "$out" ] && return 0
  curl -sS -r "0-$BYTES" --retry 3 --connect-timeout 30 "$url" 2>/dev/null \
    | zcat 2>/dev/null | head -n $(( NREAD * 4 )) | gzip -c > "$out.part"
  n=$(( $(zcat "$out.part" 2>/dev/null | wc -l) / 4 ))
  if [ "$n" -ge "$MIN" ]; then mv "$out.part" "$out"
  else rm -f "$out.part"; echo "SHORT $run: only $n reads" >&2; fi
}
export -f get; export D BYTES NREAD MIN

want=$(( $(wc -l < data/jmbc/jmbc_runs.tsv) - 1 ))
# A few runs come back empty under concurrency; each pass skips what is already there, so just go round again.
for attempt in 1 2 3; do
  have=$(ls -1 "$D/reads"/*.fastq.gz 2>/dev/null | wc -l)
  [ "$have" -ge "$want" ] && break
  echo "pass $attempt: $have / $want present"
  tail -n +2 data/jmbc/jmbc_runs.tsv | awk -F'\t' '{print $1"\t"$7}' \
    | xargs -P 4 -n 1 -d '\n' bash -c 'IFS=$'"'"'\t'"'"' read -r r u <<< "$0"; get "$r" "$u"'
done

have=$(ls -1 "$D/reads"/*.fastq.gz 2>/dev/null | wc -l)
echo "=== $have / $want runs, $(du -sh "$D/reads" | cut -f1) ==="
[ "$have" -ge "$want" ] && touch "$D/reads/.complete"
# partial success still feeds the rest of the pipeline; only a near-total failure is worth stopping for
[ "$have" -ge $(( want * 9 / 10 )) ]
