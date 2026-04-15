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
                               data for simulation
      --normal-only            Use only reads from the provided normal as
                               background data, excluding non-supporting
                               reads from the sample
  -i, --include arg            Only operate on VCF records with this value
                               present in FILTER. e.g. -i PASS. May be
                               passed multiple times.
  -e, --exclude arg            Only operate on VCF records without this
                               value present in FILTER. May be passed
                               multiple times.
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
```

Basic usage:

```bash
expos my.vcf my.bam
```
