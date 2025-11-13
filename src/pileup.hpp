#pragma once

#include <algorithm>
#include <cstdint>
#include <format>
#include <unordered_set>
#include <vector>

#include <htslib/kstring.h>
#include <htslib/sam.h>

#include "bounds.hpp"
#include "region.hpp"
#include "variant.hpp"

// nothing but C please
extern "C" {
struct pf_capture {
    htsFile *fh;
    hts_itr_t *it;
    int exclude_flag;
    int include_flag;
    int min_mapq;
};
// can attach data via bam_pileup_cd
// + bam_plp_constructor/destructor
inline int pileup_func (void *data,
                        bam1_t *b) {
    pf_capture *d = static_cast<pf_capture *> (data);
    int ret;
    // find the next good read
    while (1) {
        ret = sam_itr_next (d->fh, d->it, b);
        if (ret < 0) {
            break; // EOF/err
        }
        // TODO hardcode check for is proper pair, mate mapped
        if (!(b->core.flag & d->exclude_flag) &&
            ((b->core.flag & d->include_flag) == d->include_flag) &&
            b->core.qual >= d->min_mapq) {
            break; // found good read
        };
    }
    return ret;
};
}
// end nothing but C

// get estimated template coordinates for
// all templates supporting a variant
inline std::vector<hts_region> get_templates (htsFile *aln_fh,
                                              hts_idx_t *aln_idx,
                                              Mut mut) {
    std::vector<hts_region> tmpl_reg;
    std::unordered_set<std::string> qnames;
    bam_plp_t buf = NULL;
    bam1_t *b     = NULL;

    safe_size_opts sso_plp_pos;
    sso_plp_pos.msg = "error translating htslib pileup position into "
                      "appropriate index for results array";

    // fetch a read overlapping the query region;
    // then do a pileup per base for the total region
    // covered by the retrieved read;
    // then count events on those pileups which overlap
    // the original query region.
    hts_itr_t *iter =
        sam_itr_queryi (aln_idx, mut.rid, mut.rpos0, mut.rpos0 + 1);

    pf_capture pfc{aln_fh, iter, 0, 0,
                   0}; // not using filtering params as of now
    buf = bam_plp_init (pileup_func,
                        &pfc); // initialize pileup
    bam_plp_set_maxcnt (buf, 10000);

    // could intialise these separately to the function call since they'll be
    // reused in the same way every time
    int64_t plp_pos = -1, rpos;
    std::pair<int64_t *, int64_t *> tco;
    int plp_tid = -1, n_plp = -1, tid;
    const bam_pileup1_t *pl;
    std::string mc, pqname;
    auto mateb = bam_init1();
    std::array<int64_t, 4> endpoints;

    // using smart pointers would be best
    auto cleanup = [&] {
        sam_itr_destroy (iter);
        bam_destroy1 (b);
        bam_destroy1 (mateb);
        bam_plp_destroy (buf);
    };

    while ((pl = bam_plp64_auto (buf, &plp_tid, &plp_pos, &n_plp)) !=
           0) {
        if (n_plp < 0 || plp_tid < 0 || plp_pos < 0) {
            throw std::runtime_error ("pileup failed");
        }
        if (mut.rpos0 != plp_pos || !mut.evaluate_support (pl)) {
            continue;
        }
        pqname = bam_get_qname (pl->b);
        if (qnames.find (pqname) !=
            qnames.end()) // qname already seen
            continue;
        qnames.insert (pqname);
        tid = pl->b->core.tid;
        if (tid != pl->b->core.mtid)
            continue;

        // These are error states rather than skips
        // because if either occurs then something is fundamentally wrong
        // i.e. neither can occur during correct application of this function
        mc = bam_aux2Z (bam_aux_get (pl->b, "MC"));
        if (mc.empty()) {
            cleanup();
            throw std::runtime_error (
                std::format ("no MC tag for qname {}", pqname));
        }
        if (bam_parse_cigar (mc.c_str(), NULL, mateb) < 1) {
            cleanup();
            throw std::runtime_error (
                std::format ("unable to parse MC tag {} for qname {}",
                             mc, pqname));
        }

        //--- get template region ---//
        rpos         = pl->b->core.pos; // leftmost coord
        endpoints[0] = rpos;
        endpoints[1] = rpos + bam_cigar2rlen (pl->b->core.n_cigar,
                                              bam_get_cigar (pl->b));
        rpos         = pl->b->core.mpos; // leftmost mate coord
        endpoints[2] = rpos;
        endpoints[3] = rpos + bam_cigar2rlen (mateb->core.n_cigar,
                                              bam_get_cigar (mateb));

        tco =
            std::minmax_element (endpoints.begin(), endpoints.end());

        tmpl_reg.push_back (
            hts_region::by_end (tid, *tco.first, *tco.second));
    }

    cleanup();

    return tmpl_reg;
}
