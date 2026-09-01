# expos

[![CI](https://github.com/blex-max/expos/actions/workflows/ci.yml/badge.svg)](https://github.com/blex-max/expos/actions/workflows/ci.yml)
[![Docs](https://img.shields.io/badge/docs-online-blue)](https://blex-max.github.io/expos/)
[![C++23](https://img.shields.io/badge/C%2B%2B-23-blue)](https://en.cppreference.com/w/cpp/23)

**FULL DOCUMENTATION AVAILABLE [HERE](https://blex-max.github.io/expos/)**

Statistically-backed VCF flagging calculating effect sizes
and p-values for spatial properties of somatic mutations via
Monte Carlo simulation and local reference complexity via Lempel-Ziv 76 method.
Compares mutant reads with the set of all reads,
including reads from one or more additional background samples.

Useful for inspecting and flagging false positive variants caused
by a variety of processes most commonly associated with a high
number of PCR cycles. Builds on Ellis et al. 2021, GATK ReadPosRankSum, bcftools RPBZ, amongst others.

High-performance implementation, adding minimal overhead to existing workflows.

Applicable to SNVs, small MNVs, and small indels.

## Installation

### Requirements

- C++23-capable compiler (tested with gcc 12 and clang 21)
- CMake 3.22+
- [htslib](https://github.com/samtools/htslib) ≥ 1.14

[argparse](https://github.com/p-ranav/argparse) and [fmt](https://github.com/fmtlib/fmt) are vendored directly under `vendor/`. Tests additionally fetch [Catch2](https://github.com/catchorg/Catch2) via CMake FetchContent when `MAKE_TEST=ON` - making the test binary therefore requires network access for at least the configure step.

## Build from Source

```bash
# from the cloned repo
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
./expos --help
```

If htslib is not discoverable via `pkg-config`, point CMake to it directly:

```bash
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DHTSLIB_INCLUDE_DIR=/path/to/htslib/include \
  -DHTSLIB_LIBRARY=/path/to/libhts.so
```

Prebuilt docker images for linux/arm64 and linux/amd64 are also provided via the GitHub container registry. They can be found in the `packages` section of the GitHub page for the repo.

To also compile `estimate-entropy`, a daughter tool for complexity estimation of (genomic) strings:

```bash
cmake .. -DMAKE_DAUGHTER=ON -DCMAKE_BUILD_TYPE=Release
cmake --build .
```

## Usage

```
Usage: expos [--help] [--version] [--seed SEED] [--flank SIZE]
             [--max-frag-len LEN] [--quiet] [--skip-filtered]
             [--background-sample PATH]...
             VCF REF ALN

...

Positional arguments:
  VCF                           input VCF/BCF of variants to annotate
  REF                           indexed reference genome FASTA
  ALN                           indexed alignment (BAM/CRAM) of the sample

Optional arguments:
  -h, --help                    shows help message and exits
  -v, --version                 prints version information and exits
  --seed SEED                   random seed for the Monte-Carlo simulation [default: 24601]
  --flank SIZE                  Size of reference sequence flanks to retrieve from
                                either side of variant, for use in calcuating
                                reference complexity. It is only necessary to modify
                                from the default value if low complexity windows more
                                distant from variant loci are likley to correlate with
                                artefacts given the sequencing protocol. [default: 250]
  --max-frag-len LEN            Upper bound on fragment size for a fragment to be
                                included in analysis. Useful to avoid confounding
                                template TJAC statistic with ambiguously mapped
                                fragments with improbable length. [default: 2000]
  -q, --quiet                   suppress per-record warnings to stderr
  --skip-filtered               Only analyse records where FILTER is PASS or . (unset)
  -b, --background-sample PATH  additional indexed alignment file/s from which
                                to merge data into Monte-Carlo background.
                                Supporting reads are always taken from the
                                primary ALN only. [default: {}] [may be repeated]
```

The three positional arguments are required and consumed in the order
`VCF REF ALN`. The annotated VCF is written to stdout; the reference FASTA
and the alignment must both be indexed (`.fai`, and `.bai`/`.crai`).
Basic usage then looks like:
```bash
expos my.vcf ref.fa my.bam > annotated.vcf
```

## Development

To enable the repo's git hooks (clang-format on staged changes at commit
time), run this once after cloning:
```bash
git config core.hooksPath .githooks
```

## Acknowledgements

High-performance PRNG implementation follows [Vigna](https://prng.di.unimi.it/MWC192.c) and [Lemire](https://lemire.me/blog/2019/06/06/nearly-divisionless-random-integer-generation-on-various-systems/).
