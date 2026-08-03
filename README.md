# expos

[![CI](https://github.com/blex-max/expos/actions/workflows/ci.yml/badge.svg)](https://github.com/blex-max/expos/actions/workflows/ci.yml)
[![Docs](https://img.shields.io/badge/docs-online-blue)](https://blex-max.github.io/expos/)
[![C++23](https://img.shields.io/badge/C%2B%2B-23-blue)](https://en.cppreference.com/w/cpp/23)
[![status: alpha](https://img.shields.io/badge/status-alpha-orange)](https://github.com/blex-max/expos)

**FULL DOCUMENTATION AVAILABLE [HERE](https://blex-max.github.io/expos/)**

Statistically-backed VCF flagging calculating effect sizes
and p-values for spatial properties of somatic mutations via
Monte Carlo simulation. Compares mutant reads (supporting observations)
with the set of all reads (background), including reads
from one or more "normal" samples.

Useful for inspecting and flagging false positive variants caused
by a variety of processes most commonly associated with a high
number of PCR cycles. Builds on Ellis et al. 2021, GATK ReadPosRankSum, bcftools RPBZ, amongst others.

Applicable to SNVs, small MNVs, and small indels.

## Installation

```bash
 # from the cloned repo
 mkdir build
 cd build
 cmake .. -DCMAKE_BUILD_TYPE=Release
 cmake --build .
 ./expos --help
```

## Usage

```
Usage: expos [--help] [--version] [--seed SEED] [--uncompressed] [--quiet]
             [--skip-filtered] [--additional-background-samples PATH...]
             [--debug] VCF REF ALN

Positional arguments:
  VCF                 input VCF/BCF of variants to annotate
  REF                 indexed reference genome FASTA
  ALN                 indexed alignment (BAM/CRAM) of the sample

Optional arguments:
  -h, --help          shows help message and exits
  -v, --version       prints version information and exits
  --seed SEED         random seed for the Monte-Carlo simulation [default: 24601]
  -u, --uncompressed  write uncompressed VCF (default: bgzip-compressed)
  -q, --quiet         suppress per-record warnings to stderr
  --skip-filtered     only analyse records where FILTER is PASS or . (unset)
  --bg, --additional-background-samples PATH...
                      additional indexed alignment file(s) whose reads are
                      merged into the Monte-Carlo background. Supporting reads
                      are always taken from the primary ALN only.
  --debug             enable debug logging to stderr
```

The three positional arguments are required and consumed in the order
`VCF REF ALN`. The annotated VCF is written to stdout; the reference FASTA
and the alignment must both be indexed (`.fai`, and `.bai`/`.crai`). Basic
usage then looks like:
```bash
expos my.vcf ref.fa my.bam > annotated.vcf
```
