# Conceptual Walkthrough


expos looks for clustering in query position of variant supporting reads (and template endpoints).
Here's a simple example:

![A variant with clustered query position compared to reference](images/clustered-variant.svg)

what is the likelihood that this variant supporting subset of query positions could have been drawn
from the total set of reads?

<!-- eff-sz, pval -->
QRK=2.2006,0.00079968;
approximately 4x as clustered as the total set, with a >0.0001% probablility that this could be a subsample of the total set.

and in this case of a broad variant:

![A variant with well-distributed query position](images/broad-variant.svg)

QRK=0.0300375,0.72611;
The mutant reads show almost no clustering compared to the total set, and there is a 72% probability that they could be a subsample of the total set

and a complex example with mixed clustering:

![A more nuanced case](images/mixed-clustering.svg)

QRK=2.04051,0.00079968;

<!-- Because the total set of reads is used, rather than ref v alt, ref and alt in distinct clusters should be flagged -->
<!-- TODO: add test case for the above -->
