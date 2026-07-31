#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENTROPY_BIN="${ENTROPY_BIN:-$SCRIPT_DIR/../build/estimate-entropy}"
WINDOW=1000

usage() {
  cat >&2 <<EOF
Usage: $(basename "$0") <reference.fasta> <regions.bed>

Runs estimate-entropy on each region in the BED file and prints:
  CHROM  START  END  ENTROPY  (tab-separated, one line per region)

Regions wider than ${WINDOW} bp are tiled into non-overlapping windows;
the reported entropy is the mean across windows.

Override binary path: ENTROPY_BIN=/path/to/estimate-entropy $0 ...
EOF
  exit 1
}

[[ "$#" -eq 2 ]] || usage

REF="$1"
BED="$2"

[[ -f "$REF" ]] || { echo "Error: FASTA not found: $REF" >&2; exit 1; }
[[ -f "$BED" ]] || { echo "Error: BED not found: $BED" >&2; exit 1; }
[[ -x "$ENTROPY_BIN" ]] || { echo "Error: estimate-entropy not found: $ENTROPY_BIN" >&2; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "Error: samtools not in PATH" >&2; exit 1; }

while IFS=$'\t' read -r chrom start end rest; do
  # skip BED header/comment lines
  case "$chrom" in '#'*|track|browser) continue ;; esac
  [[ -z "$chrom" ]] && continue

  seq=$(samtools faidx "$REF" "${chrom}:$((start+1))-${end}" 2>/dev/null \
        | grep -v '^>' | tr -d '\n\r' | tr '[:lower:]' '[:upper:]')

  if [[ -z "$seq" ]]; then
    echo "Warning: no sequence for ${chrom}:${start}-${end}, skipping" >&2
    continue
  fi

  entropy=$(printf '%s\n' "$seq" \
    | awk -v w="$WINDOW" '{ for (i=1; i<=length($0); i+=w) print substr($0, i, w) }' \
    | "$ENTROPY_BIN" \
    | awk '{ sum += $1; n++ } END { printf "%.6f", (n > 0 ? sum/n : 0) }')

  printf '%s\t%s\t%s\t%s\n' "$chrom" "$start" "$end" "$entropy"
done < "$BED"
