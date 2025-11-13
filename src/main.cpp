// clang-format off
// GET TEMPLATES AND DISTRIBUTION OF TEMPLATES SUPPORTING A VARIANT
// (optionally) merge templates within +-n bp of eachother
// write (merged) templates with representative qnames + distribution/s
// into VCF
// ---
// The above in itself is probably an advancement in method for assessing variant
// legitimacy for low yield seq. If theres very little template devation (compared to the average)
// then almost certainly it's an artefact
// ---
// WITH BCFTOOLS CONSENSUS
// fetch template region from reference.
// attempt to assemble sample template from ref base,
// i.e. with mutations if any, by merging good reads
// WITH RNAlfold
// evaluate foldiness of recovered templates
// clang-format on


int main (int argc,
          char *argv[]) {
    return 0;
}
