#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <sequence.fasta> (single line fasta only)" >&2
  exit 1
fi

FASTA="$1"

awk '
/^>/ { next }
{
  buf = buf $0
  while (length(buf) >= 100) {
    window = toupper(substr(buf, 1, 100))
    if (window ~ /^[ACGT]+$/) print window
    buf = substr(buf, 2)
  }
}
' "$FASTA" | ../build/estimate-entropy > entropy_windows.tsv


datamash --header-out min 1 max 1 mean 1 median 1 sstdev 1 < entropy_windows.tsv > entropy_summary.txt

