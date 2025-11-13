#pragma once

#include <cstdint>
#include <cstring>
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>
#include <string>


struct Mut {
    int32_t rid;
    int64_t rpos0; // 0 indexed start
    std::string ref, alt;
    // bool imprecise; // for imprecise indels <>

    size_t rlen () const noexcept { return ref.size(); }

    size_t qlen () const noexcept { return alt.size(); }

    int diff () const { // TODO guard lengths in constructor between 0
                        // < x < some short length
        // if (imprecise)
        //     throw std::runtime_error (
        //         "Mut: attempt to call len_diff on imprecise
        //         indel");
        return (rlen() > qlen()) ? -(rlen() - qlen())
                                 : qlen() - rlen();
    }

    bool is_del () const { return (diff() < 0); }
    bool is_ins () const { return (diff() > 0); }
    bool is_indel () const { return (diff() != 0); }

    // TODO method of Mut
    // return 1 if bam_pileup1_t supports read
    // 0 otherwise. -1 on error
    bool evaluate_support (const bam_pileup1_t *pl) const {
        auto qpos = static_cast<size_t> (pl->qpos);

        auto qbase (seq_nt16_str[bam_get_seq (pl->b)[qpos]]);
        if (qbase != alt[0])
            return false; // all mutation types must conform

        auto indel_lmatch =
            diff() == pl->indel; // both + for ins, - for del
        if (is_indel() && !indel_lmatch)
            return false;
        if (is_del())
            return static_cast<int> (indel_lmatch);

        if (is_ins()) {
            kstring_t ins;
            ks_initialize (&ins);
            auto rc = bam_plp_insertion (pl, &ins, NULL);
            if (rc != diff()) { // failure or wrong length
                ks_free (&ins);
                return false;
            }
            auto ins_match = static_cast<int> ( // check bases match
                strcmp (ks_c_str (&ins), alt.c_str()) == 0);
            ks_free (&ins);
            return (ins_match);
        }

        if (qlen() == 1) // SNV
            return true;

        std::string mnv; // must be a MNV
        for (size_t i = 0; i < qlen(); ++i) {
            mnv.push_back (
                seq_nt16_str[bam_get_seq (pl->b)[qpos + i]]);
        }

        if (mnv == alt)
            return true;

        return false;
    }
};
