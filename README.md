# expos

Spatial information extraction for variant sites in alignment data,
plus comparison to background via monte-carlo simulation.
Useful for inspecting and flagging false postive variants.

Core functionality all present but niceties and guard rails are not.

## Installation

```bash
 # from the cloned repo
 mkdir build
 cd build
 cmake ..
 cmake --build .
 ./expos --help
```

## Example

For clarity commands are separated, but in practice these can all be piped together
```bash
 # calculate statistics and encode to VCF
 expos --ref ref.fa my.vcf my.bam > expos.vcf.gz

 # Add some soft flags in the FILTER column
 # (or alternately, subset entirely with bcftools view instead of filter)

 # monte-carlo simulation backed spatial distribution filters
 # filtering on effect size and pvalue
 bcftools filter \
 --mode + \
 -s SPATIAL \
 -e '(INFO/QM1NN[1] <= -1.0 & INFO/QM1NN[2] <= 0.05) | (INFO/TM1NN[1] <= -1.0 & INFO/TM1NN[2] <= 0.05)' \
 expos.vcf.gz > f1.vcf

 # heuristic for poorly aligned supporting reads in regions of low complexity
 bcftools filter \
 --mode + \
 -s LOW_CMPLX \
 -e 'INFO/MLAS <= 0.93 & INFO/KC < 150' \
 f1.vcf > f2.vcf
```
