#pragma once

#include <algorithm>
#include <cstdint>
#include <format>
#include <htslib/vcf.h>
#include <unordered_set>
#include <vector>

#include <htslib/kstring.h>
#include <htslib/sam.h>

#include "hts_ptr_t.hpp"
#include "stats.hpp"
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

auto inline get_aln_data (
    htsFile   *aln_fh,
    hts_idx_t *aln_idx,
    bcf1_t    *v
) {
    using templ_endpoints = std::vector<endpoints1D<uint64_t>>;
    using qposv           = std::vector<uint64_t>;
    using qdat            = struct qdat {
        qposv qpv;     // NOTE: equivalent to analysing read endpoints
        templ_endpoints tev;
    };
    struct {
        qdat alt;
        qdat other;
    } obs;

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
    std::unordered_set<std::string> qnames;
    std::array<int64_t, 4>          endpoints;
    bam1_upt                        mateb{bam_init1(), bam_destroy1};

    while ((plarr = bam_plp64_auto (buf.get(), &plp_tid, &plp_pos, &n_plp))
           != 0) {
        if (n_plp < 0 || plp_tid < 0 || plp_pos < 0)
            throw std::runtime_error ("pileup failed");

        if (plp_pos < v->pos || plp_pos >= v->pos + v->rlen) {
            continue;     // doesn't cover variant
        }

        for (size_t i = 0; i < static_cast<size_t> (n_plp); i++) {
            const auto pli = plarr + i;

            const std::string qname{bam_get_qname (pli->b)};
            // These are error states rather than skips
            // because if either occurs then something is fundamentally wrong
            // i.e. neither can occur during correct application of this function
            auto              raw_mc = bam_aux_get (pli->b, "MC");
            if (raw_mc == NULL)
                throw std::runtime_error (
                    std::format ("no MC tag for read {}", qname)
                );
            if (bam_aux_type (raw_mc) != 'Z')
                throw std::runtime_error (
                    std::format ("MC tag is not of type 'Z' for read {}. Record data corrupt; type 'Z' is mandated for MC tag by SAM format spec.", qname)
                );
            const std::string mc{bam_aux2Z (raw_mc)};
            if (bam_parse_cigar (mc.c_str(), NULL, mateb.get()) < 1) {
                throw std::runtime_error (
                    std::format ("unable to parse MC tag {} as cigar string for read {}", mc, qname)
                );
            }

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

            auto &bin = evaluate_support (pli, v, mtype) ? obs.alt
                                                         : obs.other;

            const auto l0 = pli->b->core.pos;

            bin.qpv.emplace_back (as_uint (pli->qpos));
            // don't double count templates,
            // shared between read pairs (by definition)
            if (qnames.find (qname) != qnames.end()) {     // qname already seen
                continue;
            }
            qnames.insert (qname);

            // TODO logging
            // std::cout << "retrieving template" << std::endl;

            //--- get template region ---//
            endpoints[0] = l0;
            endpoints[1] = l0
                           + bam_cigar2rlen (
                               static_cast<int> (pli->b->core.n_cigar),
                               bam_get_cigar (pli->b)
                           );
            const auto ml0 = pli->b->core.mpos;     // leftmost mate coord
            endpoints[2] = ml0;
            endpoints[3] = ml0
                           + bam_cigar2rlen (
                               static_cast<int> (mateb->core.n_cigar),
                               bam_get_cigar (mateb)
                           );

            const auto tco = std::minmax_element (
                endpoints.begin(),
                endpoints.end()
            );     // NOTE returns pair of *ptrs*

            bin.tev.emplace_back (as_uint (*tco.first), as_uint (*tco.second));
        }
    }

    return obs;
}
