#pragma once

#include <algorithm>
#include <cstdint>
#include <format>
#include <stdexcept>
#include <unordered_set>
#include <vector>

#include <htslib/kstring.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>

#include "hts_ptr_t.hpp"
#include "lib-stats/spatial.hpp"
#include "variant.hpp"
#include "util.hpp"

// nothing but C please
extern "C" {
struct pf_capture {
    htsFile   *fh;
    hts_itr_t *it;
    int        inc_flag;
    int        exc_flag;
};
inline int pileup_func (
    void   *data,
    bam1_t *b
) {
    pf_capture *d = (pf_capture *)(data);
    int         ret=-1;
    uint16_t    flag;
    // find the next good read
    while (1) {
        ret = sam_itr_next (d->fh, d->it, b);
        if (ret < 0) {
            break;     // EOF/err
        }
        flag = b->core.flag;
        if (!(flag & d->exc_flag) && ((flag & d->inc_flag) == d->inc_flag)) {
            // && b->core.qual >= d->min_mapq) {
            break;     // found good read
        };
    }
    return ret;
};
}
// end nothing but C


// Get pileup1_t covering a variant
using PileupColumn = std::vector<const bam_pileup1_t *>;
using AlnItr = hts_itr_t *;
using PileupItr = bam_plp_t;
std::tuple<PileupColumn, AlnItr, PileupItr>
inline get_pileup (
    htsFile   *aln_fh,
    hts_idx_t *aln_idx,
    bcf1_t    *var,
    int        sam_flag_include,
    int        sam_flag_exclude
) {
    auto aln_iter = sam_itr_queryi (
        aln_idx,
        var->rid,
        var->pos,
        var->pos + var->rlen
    );
    // TODO fix reported rid/positions to user expectations (and elsewhere)
    if (aln_iter == nullptr) {
        throw std::runtime_error (
            std::format (
                "could not create iterator for {}:{}-{}",
                var->rid,
                var->pos,
                var->pos + var->rlen
            )
        );
    }
    pf_capture  pfc{
        aln_fh,
        aln_iter,
        sam_flag_include,
        sam_flag_exclude
    };
    auto pileup_iter{bam_plp_init (pileup_func, &pfc)};
    bam_plp_set_maxcnt (
        pileup_iter,
        10000
    );     // TODO max depth placeholder

    int64_t              plp_pos = -1;
    int                  plp_tid = -1;
    int                  n_plp   = -1;
    const bam_pileup1_t *plarr;
    PileupColumn         pv;
    while (
        (plarr = bam_plp64_auto (pileup_iter, &plp_tid, &plp_pos, &n_plp))
        != 0
    ) {
        if (n_plp < 0 || plp_tid < 0 || plp_pos < 0)
            throw std::runtime_error ("pileup failed");

        if (plp_pos != var->pos) {
            continue;     // doesn't cover variant
        } else {
            pv.reserve((size_t)n_plp);
            for (size_t i = 0; i < static_cast<size_t> (n_plp); i++) {
                pv.push_back(plarr + i);
            }
        }
        break;
    }
    return {pv, aln_iter, pileup_iter};
}


inline PileupColumn partition_supporting (const PileupColumn &pc, const bcf1_t *var, int mutation_type) {
    PileupColumn out;
    for (const auto p : pc) {
        if (evaluate_support(p, var, mutation_type)) {
            out.push_back(p);
        }
    }
    return out;
}

struct PileupMetrics {
    std::vector<uint64_t> query_position;
    std::vector<double>   normalised_as;
    std::vector<spatial::line_seg> template_endpoints;
    size_t                nreads = 0;
};
inline PileupMetrics get_metrics (const PileupColumn &pc) {
    PileupMetrics out;

    std::unordered_set<const char*> qnames;
    std::array<int64_t, 4>          endpoints;
    bam1_upt                        mateb{bam_init1(), bam_destroy1};
    for (const auto p : pc) {
        const auto qname = bam_get_qname (p->b);

        // These are error states rather than skips
        // because if either occurs then something is fundamentally wrong
        // i.e. neither can occur during correct appcation of this function
        const auto raw_mc = bam_aux_get (p->b, "MC");
        if (raw_mc == NULL)
            throw std::runtime_error (
                std::format ("no MC tag for read {}", qname)
            );
        if (bam_aux_type (raw_mc) != 'Z')
            throw std::runtime_error (
                std::format (
                    "MC tag is not of type 'Z' for read {}. "
                    "Record data corrupt; type 'Z' is mandated "
                    "for MC tag by SAM format spec. Try samtools"
                    "fixmate?",
                    qname
                )
            );
        const auto mc{bam_aux2Z (raw_mc)};
        if (bam_parse_cigar (mc, NULL, mateb.get()) < 1) {
            throw std::runtime_error (
                std::format (
                    "unable to parse MC tag {} as cigar string "
                    "for read {}",
                    mc,
                    qname
                )
            );
        }

        auto raw_AS = bam_aux_get (p->b, "AS");
        if (raw_AS == NULL)
            throw std::runtime_error (
                std::format ("no AS tag for read {}", qname)
            );
        const auto raw_AS_tag_type = bam_aux_type (raw_AS);
        if (raw_AS_tag_type != 'i' && raw_AS_tag_type != 'C')
            throw std::runtime_error (
                std::format (
                    "AS tag is not of type 'i' for read {}. "
                    "Record data corrupt; type 'i' is mandated "
                    "for AS tag by SAM format spec.",
                    qname
                )
            );

        out.nreads++;

        out.normalised_as.push_back (
            static_cast<double> (bam_aux2i (raw_AS))
            / static_cast<double> (p->b->core.l_qseq)
        );     // length-normalised alignment score

        // get qpos
        if (!(p->is_del || p->is_refskip || p->qpos < 0)) {
            out.query_position.emplace_back (as_uint (p->qpos));
        }

        // don't double count templates
        if (qnames.find (qname)
            != qnames.end()) {     // qname already seen
            continue;
        }
        qnames.insert (qname);

        // check mate mapped to same reference
        // NOTE: this would need to be adjusted
        // if single end data is to be accepted
        if (p->b->core.tid != p->b->core.mtid) {
            continue;
        }

        //--- get template region ---//
        endpoints[0] = p->b->core.pos;
        endpoints[1] = endpoints[0]
                       + bam_cigar2rlen (
                           static_cast<int> (p->b->core.n_cigar),
                           bam_get_cigar (p->b)
                       );
        endpoints[2] = p->b->core.mpos;  // leftmost mate coordinate
        endpoints[3] = endpoints[2]
                       + bam_cigar2rlen (
                           static_cast<int> (mateb->core.n_cigar),
                           bam_get_cigar (mateb)
                       );

        const auto tco = std::minmax_element (
            endpoints.begin(),
            endpoints.end()
        );     // NOTE returns pair of *ptrs*

        out.template_endpoints.emplace_back (
            as_uint (*tco.first),
            as_uint (*tco.second)
        );
    }

    return out;
}


// returns metrics for supporting, and metrics for all
inline std::pair<PileupMetrics, PileupMetrics> pileup_partition_and_anaylse (
    htsFile   *aln_fh,
    hts_idx_t *aln_idx,
    bcf1_t    *var,
    int        mutation_type,
    int        sam_flag_include,
    int        sam_flag_exclude
) {
    auto [pileup_all,
          aln_iter,
          pileup_buf] = get_pileup (
                aln_fh,
                aln_idx,
                var,
                sam_flag_include,
                sam_flag_exclude
            );
    auto pileup_supporting = partition_supporting(pileup_all, var, mutation_type);

    // supporting, total
    std::pair<PileupMetrics, PileupMetrics> out;
    // we analyse both independently for overlap reasons
    out.first = get_metrics(pileup_supporting);
    out.second = get_metrics(pileup_all);

    hts_itr_destroy(aln_iter);
    bam_plp_destroy(pileup_buf);
    return out;
}

inline PileupMetrics pileup_analyse (
    htsFile   *aln_fh,
    hts_idx_t *aln_idx,
    bcf1_t    *var,
    int        sam_flag_include,
    int        sam_flag_exclude
) {
    auto [pileup_all,
          aln_iter,
          pileup_buf] = get_pileup (
                aln_fh,
                aln_idx,
                var,
                sam_flag_include,
                sam_flag_exclude
            );
    auto out = get_metrics(pileup_all);
    hts_itr_destroy(aln_iter);
    bam_plp_destroy(pileup_buf);
    return out;
}
