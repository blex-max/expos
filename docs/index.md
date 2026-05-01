# expos

Statistically-backed VCF flagging calculating effect sizes
and p-values for spatial properties of somatic mutations via
Monte Carlo simulation. Compares mutant reads (supporting observations)
with the set of all reads (background), including reads
from one or more "normal" samples.

Useful for inspecting and flagging false positive variants caused
by a variety of processes most commonly associated with a high
number of PCR cycles. Builds on Ellis et al. 2021, GATK ReadPosRankSum, bcftools RPBZ, amongst others.

Applicable to SNVs, small MNVs, and small indels.

For a full introduction see [Concepts & Theory](concepts.md)

---

## **Quick Start**

```bash
# clone and build
git clone https://github.com/blex-max/expos-repo.git
cd expos-repo
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .

# run
./expos my.vcf my.bam
```

See [Installation](installation.md) for full details and [Usage](usage.md) for all options.
