# Usage

<!-- TODO: update this page with any new options as they are added -->

```
EXtract POSitional data and statistics from alignment at
VCF variant sites, and encode them as INFO fields to VCF.
Requires the presence of .(b/cr)ai indexes of the same name
as the relevant alignment. Annotated VCF to stdout. See
README or output VCF header for descriptions of fields
added.

Usage:
  expos [OPTION...] <VCF/BCF (- for stdin)> <ALN.(b/cr)am>

  -h, --help                   Print usage
  -l, --expected-read-len arg  Sequencing read length. Default: 150
  -r, --ref arg                Alignment Reference Fasta for optionally
                               adding reference complexity to statistics.
  -n, --normal arg             Alignment for use as additional background
                               data for simulation. May be passed multiple
                               times.
      --normal-only            Use only reads from the provided normal as
                               background data, excluding reads from the
                               sample
  -i, --include-records arg    Only operate on VCF records with this value
                               present in FILTER. e.g. -i PASS. May be
                               passed multiple times.
  -e, --exclude-records arg    Only operate on VCF records without this
                               value present in FILTER. May be passed
                               multiple times.
  -I, --include-reads arg      SAM flag: only include reads with all of
                               these bits set. Set 0 to disable. Default: 3
  -E, --exclude-reads arg      SAM flag: exclude reads with any of these
                               bits set. Default: 3852
  -t, --tsv arg                Write a tsv of extended statistics to file
                               specified.
  -u, --uncompressed           output uncompressed VCF
      --seed arg               Set random seed. Default: 24601
      --uniform                additionally simulate against uniform null
                               model for query position, and add result to
                               --tsv output. For assessment of correlation
                               with simulation against all-reads null.
      --assess-microhomology   additionally assess STR and homopolymer
                               content of reference regions, and add result
                               to --tsv output. For assessment of
                               correlation with drop in LZ.
      --debug                  Enable debug logging
```

Basic usage:

```bash
expos my.vcf my.bam
```
