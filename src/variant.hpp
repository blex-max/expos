#pragma once

#include <cstring>
#include <format>
#include <stdexcept>
#include <string>

#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>


// in caller, to assess type
// bcf_has_variant_types (b, VCF_DEL | VCF_INS | VCF_SNP | VCF_MNP, bcf_variant_match::bcf_match_overlap)
// preconditions:
// -- normalised bcf1_t (single alt)
// -- bam_pileup1_t covers variant
inline bool evaluate_support (
    const bam_pileup1_t *pl,
    const bcf1_t        *b,
    int                  mtype
) {
    assert (b->n_allele == 2);

    // TODO (ask samteam/read source)
    // if (!b->unpacked)

    std::string ref, alt, mnv, qbase;
    size_t      qpos;
    int         vdiff;
    bool        indel_lmatch;

    switch (mtype) {
        case (VCF_DEL):
        case (VCF_INS):
        case (VCF_SNP):
        case (VCF_MNP):
            ref = b->d.allele[0];
            alt = b->d.allele[1];

            qpos = static_cast<size_t> (pl->qpos);
            qbase = seq_nt16_str[bam_seqi (bam_get_seq (pl->b), qpos)];
            // std::cout << std::format ("ref {}, alt {}, query {}", ref, alt, qbase)
            // << std::endl;
            if (qbase
                != alt.substr (
                    0,
                    1
                ))     // compare to first character of mutation
                return false;     // all mutation types must conform
            if (mtype
                == VCF_SNP)     // only relevant condition for SNP
                return true;

            vdiff = static_cast<int> (alt.length())
                    - static_cast<int> (
                        ref.length()
                    );     // negative for del, pos for ins
            indel_lmatch = (vdiff == pl->indel);

            // needs verification
            if ((mtype & (VCF_INS | VCF_DEL)) && !indel_lmatch)
                return false;
            if (mtype == VCF_DEL)
                return indel_lmatch;     // i.e. true due to the previous conditional

            // needs verification
            if (mtype == VCF_INS) {
                kstring_t ins;
                ks_initialize (&ins);
                auto rc = bam_plp_insertion (pl, &ins, NULL);
                if (rc != vdiff) {     // failure or wrong length
                    ks_free (&ins);
                    return false;
                }
                auto ins_match
                    = static_cast<int> ( // check bases match
                        strcmp (ks_c_str (&ins), alt.c_str()) == 0);
                ks_free (&ins);
                return (ins_match);
            }

            // only mnv remains
            for (size_t i = 0; i < alt.length(); ++i) {
                mnv.push_back (
                    seq_nt16_str[bam_seqi (bam_get_seq (pl->b), qpos + i)]
                );
            }

            return (mnv == alt);
        default:
            throw std::runtime_error (
                std::format (
                    "type {} of variant {} does not match one and "
                    "only one of VCF_DEL, "
                    "VCF_INS, VCF_SNP, VCF_MNP",
                    mtype,
                    b->d.id
                )
            );
    }
}
