#pragma once

#include <algorithm>
#include <cstdint>
#include <format>
#include <htslib/vcf.h>
#include <iostream>
#include <unordered_set>
#include <vector>

#include <htslib/kstring.h>
#include <htslib/sam.h>

#include "hts_ptr_t.hpp"
#include "variant.hpp"

// nothing but C please
extern "C" {
struct pf_capture {
    htsFile   *fh;
    hts_itr_t *it;
    int        min_mapq = 0;
};
// NOTE can attach data to pileup members
// via bam_pileup_cd and bam_plp_constructor/destructor.
// NOTE could evalute support here
// but would need to make sure that function is C only.
inline int pileup_func (
    void   *data,
    bam1_t *b
) {
    pf_capture *d = static_cast<pf_capture *> (data);
    int         ret;
    uint16_t    flag;
    // find the next good read
    while (1) {
        ret = sam_itr_next (d->fh, d->it, b);
        if (ret < 0) {
            break;     // EOF/err
        }
        flag = b->core.flag;
        if (!(flag & 3852) && ((flag & 3) == 3)
            && b->core.qual >= d->min_mapq) {
            break;     // found good read
        };
    }
    return ret;
};
}
// end nothing but C

struct read_template_s {
    std::string qname;
    int64_t     start;
    int64_t     end;
    size_t      len;
};

// get estimated template coordinates for
// all templates supporting a variant.
// i.e. fetch pileup overlapping the variant;
// check pileup reads support that variant;
// then get fragment endpoints.
inline std::vector<read_template_s> get_templates (
    htsFile   *aln_fh,
    hts_idx_t *aln_idx,
    bcf1_t    *v
) {
    std::vector<read_template_s> tmpls;     // return

    // TODO enum
    auto mtype = bcf_has_variant_type (
        v,
        1,
        VCF_DEL | VCF_INS | VCF_SNP | VCF_MNP
    );     // TODO and if not one of these?
    // std::cout << std::format ("mtype {}", mtype) << std::endl;

    // prepare to pileup
    hts_itr_upt iter{
        sam_itr_queryi (aln_idx, v->rid, v->pos, v->pos + v->rlen),
        hts_itr_destroy
    };
    pf_capture pfc{aln_fh, iter.get()};     // not using mapq at present
    bam_plp_upt buf{bam_plp_init (pileup_func, &pfc), bam_plp_destroy};     // initialize pileup
    bam_plp_set_maxcnt (buf.get(), 10000);     // TODO max depth placeholder

    // loop vars
    // htslib
    int64_t                         plp_pos = -1;
    int                             plp_tid = -1;
    int                             n_plp   = -1;
    const bam_pileup1_t            *plarr;
    // this
    int64_t                         rpos;     // reference position
    std::unordered_set<std::string> qnames;
    std::string                     mc, qname;
    std::array<int64_t, 4>          endpoints;
    std::pair<int64_t *, int64_t *> tco;
    bam1_upt                        mateb{bam_init1(), bam_destroy1};

    while ((plarr = bam_plp64_auto (buf.get(), &plp_tid, &plp_pos, &n_plp))
           != 0) {
        if (n_plp < 0 || plp_tid < 0 || plp_pos < 0)
            throw std::runtime_error ("pileup failed");

        if (plp_pos < v->pos || plp_pos > v->pos + v->rlen) {
            // std::cout << "pos skip" << std::endl;
            continue;     // doesn't cover variant
        }
        // std::cout << std::format ("pos {}", plp_pos) << std::endl;

        for (size_t i = 0; i < static_cast<size_t> (n_plp); i++) {
            const bam_pileup1_t *pli = plarr + i;
            // These are error states rather than skips
            // because if either occurs then something is fundamentally wrong
            // i.e. neither can occur during correct application of this function
            qname       = std::string (bam_get_qname (pli->b));
            auto raw_mc = bam_aux_get (pli->b, "MC");
            if (raw_mc == NULL)
                throw std::runtime_error (
                    std::format ("no MC tag for qname {}", qname)
                );
            if (bam_aux_type (raw_mc) != 'Z')
                throw std::runtime_error (
                    std::format ("MC tag is not of type 'Z' for qname {}. Record data corrupt; type 'Z' is mandated for MC tag by SAM format spec.", qname)
                );
            mc = std::string (bam_aux2Z (raw_mc));
            if (bam_parse_cigar (mc.c_str(), NULL, mateb.get()) < 1) {
                throw std::runtime_error (
                    std::format ("unable to parse MC tag {} as cigar string for qname {}", mc, qname)
                );
            }

            // avoid reporting for the same qname more than once
            if (qnames.find (qname) != qnames.end()) {     // qname already seen
                // std::cout << "qname skip" << std::endl;
                continue;
            }
            qnames.insert (qname);

            // check mate mapped to same reference
            // NOTE: this would probably need to be adjusted
            // if single end data is to be accepted
            if (plp_tid != pli->b->core.mtid) {
                // std::cout << "mtid skip" << std::endl;
                continue;
            }

            // check variant support
            if (!evaluate_support (pli, v, mtype)) {
                // std::cout << "support skip" << std::endl;
                continue;
            }

            // TODO logging
            // std::cout << "retrieving template" << std::endl;

            //--- get template region ---//
            rpos         = pli->b->core.pos;     // leftmost coord
            endpoints[0] = rpos;
            endpoints[1] = rpos
                           + bam_cigar2rlen (
                               static_cast<int> (pli->b->core.n_cigar),
                               bam_get_cigar (pli->b)
                           );
            rpos = pli->b->core.mpos;     // leftmost mate coord
            endpoints[2] = rpos;
            endpoints[3] = rpos
                           + bam_cigar2rlen (
                               static_cast<int> (mateb->core.n_cigar),
                               bam_get_cigar (mateb)
                           );

            tco = std::minmax_element (
                endpoints.begin(),
                endpoints.end()
            );     // NOTE returns pair of *ptrs*

            tmpls.emplace_back (
                qname,
                *tco.first,
                *tco.second,
                *tco.second - *tco.first
            );
        }
    }

    return tmpls;
}


// inline read_template_s get_template (const bam_pileup1_t *p) {}
