#! /bin/bash

set -euxo pipefail

unsorted_sam=$1
ref=$2

td=$(mktemp --directory)
sorted_sam=${td}/pileup.sorted.sam
samtools sort -Obam $unsorted_sam > $sorted_sam
samtools index $sorted_sam

bcftools mpileup -f $ref $sorted_sam | bcftools call -mv | bcftools norm -a -m- | /Users/ab63/sanger/expos/expos-repo/build/expos --debug -I 0 -u --expected-read-len 50 - $sorted_sam
