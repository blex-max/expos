# Conceptual Walkthrough


expos looks for clustering in query position of variant supporting reads (and template endpoints).
Here's a simple example:

![A clustered variant compared to reference](images/broad-variant.svg)

what is the likelihood that this variant supporting subset of query positions could have been drawn
from the total set of reads?


and in this case of a broad variant:


and a complex example with mixed clustering:



<!-- Because the total set of reads is used, rather than ref v alt, ref and alt in distinct clusters should be flagged -->
<!-- TODO: add test case for the above -->
