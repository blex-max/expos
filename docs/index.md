# expos

Statistically-backed VCF flagging calculating effect sizes
and p-values for spatial properties of somatic mutations via
Monte Carlo simulation and local reference complexity via Lempel-Ziv 76 method.
Compares mutant reads with the set of all reads, including reads
from one or more background-only samples.

Useful for inspecting and flagging false positive variants caused
by a variety of processes most commonly associated with a high
number of PCR cycles. Builds on Ellis et al. 2021, GATK ReadPosRankSum, bcftools RPBZ, amongst others.

Applicable to SNVs, small MNVs, and small indels.

High-performance implementation, adding minimal overhead to existing workflows.

---

## **Quick Start**

```bash
# clone and build
git clone https://github.com/blex-max/expos.git
cd expos
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .

# run
./expos my.vcf ref.fa my.bam > annotated.vcf
```

See [Installation](installation.md) for full details and [Usage](usage.md) for all options.
