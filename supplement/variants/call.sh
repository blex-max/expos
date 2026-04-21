#! /bin/bash

set -euxo pipefail

workdir=$1
echo $workdir

cd $workdir

rm pileup-ref.fa.fai pileup-sorted.sam calls.vcf || true
samtools sort -u -Osam pileup.sam > pileup-sorted.sam
bcftools mpileup -f pileup-ref.fa pileup-sorted.sam | bcftools call -mv | bcftools norm -a -m- > calls.vcf
