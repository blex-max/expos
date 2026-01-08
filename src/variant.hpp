#pragma once

#include <cstring>
#include <format>
#include <stdexcept>
#include <string>
#include <vector>

#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>


inline std::vector<bool> has_filters (
    bcf_hdr_t                *hdr,
    bcf1_t                   *rec,
    std::vector<std::string> &flt
) {
    std::vector<bool> out;
    for (const auto &f : flt) {
        out.push_back (
            bcf_has_filter (hdr, rec, const_cast<char *> (f.c_str()))
            > 0
        );
    }
    return out;
}


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

    std::string ref, alt, mnv, qbase, abase;
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

            // at refskip/is_del bases
            // qpos will be -1.
            // These cannot support a variant
            // since all variants are anchored by a real base
            if (pl->qpos < 0) {
                // std::cerr << std::format ("no material base")
                // << std::endl;
                return false;
            }
            qpos = static_cast<size_t> (pl->qpos);
            qbase = seq_nt16_str[bam_seqi (bam_get_seq (pl->b), qpos)];
            abase = alt.substr (0, 1);
            // std::cerr << std::format (
            //     "qpos {} char at qpos {}, "
            //     "alt base {}",
            //     qpos,
            //     qbase,
            //     abase
            // ) << std::endl;
            if (qbase
                != abase) {     // compare to first character of mutation
                // std::cerr << "qbase does not match" << std::endl;
                return false;     // all mutation types must conform
            }
            if (mtype
                == VCF_SNP)     // only relevant condition for SNP
                return true;

            vdiff = static_cast<int> (alt.length())
                    - static_cast<int> (
                        ref.length()
                    );     // negative for del, pos for ins
            indel_lmatch = (vdiff == pl->indel);

            // std::cerr << std::format (
            //     "vdiff {} rlen {} -- {} qpos {} cigar {} indel {} lmatch {}",
            //     vdiff,
            //     b->rlen,
            //     bam_get_qname (pl->b),
            //     qpos,
            //     bam_cigar_opchr (bam_get_cigar (pl->b)[pl->cigar_ind]),
            //     pl->indel,
            //     indel_lmatch
            // ) << std::endl;

            if ((mtype & (VCF_INS | VCF_DEL)) && !indel_lmatch)
                return false;
            if (mtype == VCF_DEL)
                return indel_lmatch;     // i.e. true due to the previous conditional

            // verified, but know that the VCF spec is not strict enough
            // to allow for perfect detection of matching insertions
            if (mtype == VCF_INS) {
                kstring_t ins;
                ks_initialize (&ins);
                // TODO check del len
                auto rc = bam_plp_insertion (pl, &ins, NULL);
                if (rc != vdiff) {     // failure or wrong length
                    ks_free (&ins);
                    return false;
                }

                // on the basis that insertions seem to go
                // anchor-insertion-rest. And that this gets the most success
                // -- but check with samtools team
                const auto var_ins = alt.substr (1, rc);
                // std::cerr << std::format("ins seq {}, len {} - alt ins seq {}, len {}", ins.s, rc, var_ins, var_ins.size()) << std::endl;

                auto ins_match = strcmp (
                                     ks_c_str (&ins),
                                     var_ins.c_str()
                                 )
                                 == 0;     // check bases match
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
