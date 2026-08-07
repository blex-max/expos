#! /bin/bash

# Call variants on one of the synthetic pileups in generate-out/ and annotate
# them with expos. Writes the annotated VCF to stdout.
#
#   ./call.sh generate-out/clustered-variant.pileup.sam generate-out/ref.fa
#
# The reference must be faidx-indexed (generate.py writes ref.fa.fai).

set -euo pipefail

if [ "$#" -lt 2 ]; then
  echo "usage: $0 <pileup.sam> <ref.fa> [expos-binary]" >&2
  exit 2
fi

unsorted_sam=$1
ref=$2
expos=${3:-"$(cd "$(dirname "$0")/../.." && pwd)/build/expos"}

if [ ! -x "$expos" ]; then
  echo "no expos binary at $expos -- build it, or pass one as \$3" >&2
  exit 2
fi

# BSD mktemp needs the -t form; GNU's --directory is not portable to macOS.
td=$(mktemp -d "${TMPDIR:-/tmp}/expos-call.XXXXXX")
trap 'rm -rf "$td"' EXIT

sorted_bam=${td}/pileup.sorted.bam
samtools sort -Obam "$unsorted_sam" > "$sorted_bam"
samtools index "$sorted_bam"

vcf=${td}/calls.vcf
bcftools mpileup -f "$ref" "$sorted_bam" \
  | bcftools call -mv \
  | bcftools norm -a -m- > "$vcf"

# Positional order is VCF REF ALN. The VCF is passed as a file rather than on
# stdin because expos also needs the indexed alignment alongside it.
"$expos" -u "$vcf" "$ref" "$sorted_bam"
