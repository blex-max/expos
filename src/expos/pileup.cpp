#include "pileup.hpp"

#include <algorithm>
#include <unordered_set>

#include <htslib/sam.h>


extern "C" {
struct PileupCapture {
    htsFile* fh;
    hts_itr_t* it;
    int inc_flag;
    int exc_flag;
};
int pileup_func (
    void* data,
    bam1_t* b
) {
    PileupCapture* d = (PileupCapture*)(data);
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
}
}


// non-owning view into bam_plp_t
using PileupColumn = std::vector<const bam_pileup1_t*>;

PileupColumn
get_pileup (
    hts_pos_t mutant_gpos,
    bam_plp_t piter  // ptr to bam_plp_s
)
{
    PileupColumn out;
    int64_t plp_pos = -1;
    int plp_tid = -1;
    int n_plp   = -1;
    const bam_pileup1_t* plarr;
    while (
        (plarr = bam_plp64_auto (piter, &plp_tid, &plp_pos, &n_plp))
        != 0
    ) {
        if (n_plp < 0 || plp_tid < 0 || plp_pos < 0)
            throw std::runtime_error ("pileup failed");

        if (plp_pos < mutant_gpos) {
            continue;     // doesn't cover variant
        }
        else if (plp_pos == mutant_gpos) {
            const auto nread = static_cast<size_t> (n_plp);
            out.reserve (nread);
            for (size_t i = 0; i < nread; ++i) {
                out.push_back(plarr + i);
            }
        }
        break;
    }
    return out;
}


double get_alignment_score_from_tag (const bam_pileup1_t* p) {
    const auto qname = bam_get_qname (p->b);

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

    // normalise
    return static_cast<double> (bam_aux2i (raw_AS))
           / static_cast<double> (p->b->core.l_qseq);
}

PileupMetrics
get_metrics_single_end
(const PileupColumn& pc)
{
    PileupMetrics out;

    for (const auto p : pc) {
        out.normalised_as.push_back (
            get_alignment_score_from_tag (p)
        );

        // get qpos
        if (!(p->is_del || p->is_refskip || p->qpos < 0)) {
            out.query_position.emplace_back (as_uint (p->qpos));
        }

        // get spanned region
        out.template_endpoints.emplace_back (
            as_uint (p->b->core.pos),
            as_uint (p->b->core.pos
                     + bam_cigar2rlen (
                         static_cast<int> (p->b->core.n_cigar),
                         bam_get_cigar (p->b)
                     ))
        );
    }
    return out;
}


PileupMetrics
get_metrics_paired_end
(const PileupColumn& pc)
{
    PileupMetrics out;

    std::unordered_set<const char*> qnames;
    std::array<int64_t, 4> endpoints;
    auto mateb = bam_init1();
    for (const auto p : pc) {
        const auto qname = bam_get_qname (p->b);

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
        if (bam_parse_cigar (mc, NULL, mateb) < 1) {
            throw std::runtime_error (
                std::format (
                    "unable to parse MC tag {} as cigar string "
                    "for read {}",
                    mc,
                    qname
                )
            );
        }

        out.normalised_as.push_back (
            get_alignment_score_from_tag(p)
        );

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

        if (p->b->core.tid != p->b->core.mtid
                || p->b->core.flag & BAM_FMUNMAP) {
            continue;
        }

        // get spanned region
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
    bam_destroy1 (mateb);
    return out;
}


PileupMetrics get_metrics
(const PileupColumn& pc, bool single_end=false)
{
    if (single_end) {
        return get_metrics_single_end (pc);
    }
    else {
        return get_metrics_paired_end (pc);
    }
}


PileupMetrics
pileup_analyse (
    htsFile* aln_fh,
    hts_idx_t* aln_idx,
    bcf1_t* var,
    int sam_flag_include,
    int sam_flag_exclude,
    const bcf_hdr_t* vcf_hdr,
    bool single_end
)
{
    auto aln_iter = sam_itr_queryi (
        aln_idx,
        var->rid,
        var->pos,
        var->pos + var->rlen
    );
    if (aln_iter == nullptr) {
        throw std::runtime_error (
            std::format (
                "could not create iterator for {}:{}-{}",
                bcf_hdr_id2name (vcf_hdr, var->rid),
                var->pos + 1,
                var->pos + 1 + var->rlen
            )
        );
    }
    PileupCapture  pfc{
        aln_fh,
        aln_iter,
        sam_flag_include,
        sam_flag_exclude
    };
    auto piter = bam_plp_init (pileup_func, &pfc);

    auto pile = get_pileup (var->pos, piter);
    auto out = get_metrics (pile, single_end);

    hts_itr_destroy (aln_iter);
    bam_plp_destroy (piter);

    return out;
}


PileupColumn filter
(const PileupColumn& pc, std::function<bool(const bam_pileup1_t*)> predicate_fn)
{
    PileupColumn out;
    for (const auto p : pc) {
        if (predicate_fn (p)) {
            out.push_back (p);
        }
    }
    return out;
}

GroupedMetrics
pileup_group_and_anaylse (
    htsFile* aln_fh,
    hts_idx_t* aln_idx,
    bcf1_t* var,
    int sam_flag_include,
    int sam_flag_exclude,
    std::function<bool (const bcf1_t*, const bam_pileup1_t*)> support_fn,
    const bcf_hdr_t* vcf_hdr,
    bool single_end
)
{
    auto aln_iter = sam_itr_queryi (
        aln_idx,
        var->rid,
        var->pos,
        var->pos + var->rlen
    );
    if (aln_iter == nullptr) {
        throw std::runtime_error (
            std::format (
                "could not create iterator for {}:{}-{}",
                bcf_hdr_id2name (vcf_hdr, var->rid),
                var->pos + 1,
                var->pos + 1 + var->rlen
            )
        );
    }
    PileupCapture  pfc{
        aln_fh,
        aln_iter,
        sam_flag_include,
        sam_flag_exclude
    };
    auto pileup_iter = bam_plp_init (pileup_func, &pfc);

    auto pile_all = get_pileup (var->pos, pileup_iter);

    auto pile_supporting =
        filter (pile_all, [&support_fn, var] (const auto p) { return support_fn (var, p); });

    GroupedMetrics out;
    // we analyse both independently since each
    // group may have different sets of
    // overlapping reads
    out.supporting = get_metrics (pile_supporting, single_end);
    out.all = get_metrics (pile_all, single_end);

    hts_itr_destroy (aln_iter);
    bam_plp_destroy (pileup_iter);

    return out;
}


PileupMetrics merge_pileup_metrics (PileupMetrics a, const PileupMetrics &b) {
    a.query_position.insert (end (a.query_position), begin (b.query_position), end (b.query_position));
    a.normalised_as.insert (end (a.normalised_as), begin (b.normalised_as), end (b.normalised_as));
    a.template_endpoints.insert (end (a.template_endpoints), begin (b.template_endpoints), end (b.template_endpoints));
    return a;
}


