// clang-format off
// GET TEMPLATES AND DISTRIBUTION OF TEMPLATES SUPPORTING A VARIANT
// for variant, report:
// - template start total range and MAD [Median Absolute Distribution]*
// - template size total range and MAD*
// - leftmost template start
// - rightmost template start
// - leftmost template end
// - rightmost template end
// optional reports to separate file:
// - template details per qname
// *If MAD ~= range, broad distribution, if range >> MAD, narrow distribution with outlier/s
// Using MAD because it is a robust spread statistic (i.e. not sensitive to outliers)
// and better able to deal with a low number of observations than IQR.
// NOTE could also report number of reads, bases, a pileup events style report (but different scope)
// NOTE (optionally) merge templates within +-n bp of eachother
// write (merged) templates with representative qnames + distribution/s
// into VCF
// ---
// The above in itself is probably an advancement in method for assessing variant
// legitimacy for low yield seq. If theres very little template start/size
// variation (compared to the average) then almost certainly it's an artefact.
// That is, you can filter on this data alone
// ---
// A subsequent pipeline for folding investigation might go:
// WITH BCFTOOLS CONSENSUS
// fetch template region from reference.
// attempt to assemble sample template from ref base,
// i.e. with mutations if any, by merging good reads
// WITH RNAlfold
// evaluate foldiness of recovered templates
// clang-format on

// ROUND 2
// MAD is in fact a proxy for distance between observations, from which we can
// draw conclusions about clustering
// in other words artefactual-variant relevant statistics are
// qpos clustering, via five number summary of the 1D gaps -- NOTE: 5 number summary usefulness in doubt
// NOTE: qpos clustering is equivalent to read endpoint clustering iff read lengths are ~all the same
// which they are for short read seq
// template clustering, via five number summary of the template endpoint chebyshev distance MST edges
// of these 5 number summaries, IQR is a good way of assessing clustering
// CRITICAL NOTE:
// -> maybe we can add a whole bunch of statistical power by comparing reference templates
// to supporting templates (and elsewhere), having appropriately adjusted for number of samples?
// -> my intution is, the distribution of non-supporting reads
// should match the distribution stats of supporting reads
// if nothing odd is going on
// -> shannon entropy of the reference/consensus would be useful


#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>

#include <cxxopts.hpp>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <limits>

#include "hts_ptr_t.hpp"
#include "pileup.hpp"
#include "stats.hpp"

// TODO allow exclusion of variants by filter flags
// TODO allow user defined samflags for include/exclude
// TODO support single ended?
// TODO options to output to VCF and/or TSV (VCF to stdout)
int main (
    int   argc,
    char *argv[]
) {
    namespace fs = std::filesystem;
    fs::path vcf_path;
    fs::path aln_path;

    cxxopts::Options options (
        "tempex",
        "variant fragment template extraction and statistics\n"
    );

    // clang-format off
    options.add_options() ("h,help", "Print usage")
        ("vcf", "VCF", cxxopts::value<fs::path>())
        ("aln", "Sample BAM", cxxopts::value<fs::path>());
        // ("-n,--name");  // key to match sample to VCF
    // clang-format on

    options.parse_positional ({"vcf", "aln"});
    options.positional_help ("<VCF> <ALN>");

    try {
        auto result = options.parse (argc, argv);

        if (result.count ("help")) {
            std::cout << options.help() << std::endl;
            return 0;
        }

        if (!result.count ("vcf") || !result.count ("aln"))
            throw std::runtime_error (
                "All positional arguments must be provided"
            );

        vcf_path = result["vcf"].as<fs::path>();
        aln_path = result["aln"].as<fs::path>();

        if (!fs::exists (vcf_path)) {
            throw std::runtime_error (
                "VCF file not found: " + vcf_path.string()
            );
        }

        std::cerr << "Using VCF: " << vcf_path << std::endl;

        if (!fs::exists (aln_path)) {
            throw std::runtime_error (
                "Alignment file not found: " + aln_path.string()
            );
        }

        std::cerr << "Using aln: " << aln_path << std::endl;
    } catch (const std::exception &e) {
        std::cerr << "Error parsing CLI options: " << e.what() << "\n";
        return 1;
    }

    auto _ain{hts_open (aln_path.c_str(), "r")};
    if (_ain == NULL) {
        std::cerr << std::format (
            "Could not open alignment file at {}",
            aln_path.string()
        ) << std::endl;
        return 1;
    }
    auto _aixin{sam_index_load (_ain, aln_path.c_str())};
    if (_aixin == NULL) {
        std::cerr << std::format (
            "Coud not open index for alignment "
            "file. Searched for {}.bai",
            aln_path.c_str()
        ) << std::endl;
    }

    auto _vin{hts_open (vcf_path.c_str(), "r")};
    if (_vin == NULL) {
        std::cerr << std::format (
            "Could not open VCF file at {}",
            vcf_path.string()
        ) << std::endl;
        return 1;
    }
    auto _vh{bcf_hdr_read (_vin)};
    if (_vh == NULL) {
        std::cerr << std::format (
            "Could not read header of VCF file at {}",
            vcf_path.string()
        );
    }

    htsFile_upt ap{std::move (_ain), hts_close};
    hts_idx_upt apit{std::move (_aixin), hts_idx_destroy};

    htsFile_upt vp{std::move (_vin), hts_close};
    bcf_hdr_upt vph{std::move (_vh), bcf_hdr_destroy};

    bcf1_upt b1{bcf_init(), bcf_destroy};

    // loop vars
    std::vector<endpoints1D<uint64_t>> tendpoints;
    std::optional<double>              start_mad, size_mad;
    while (bcf_read (vp.get(), vph.get(), b1.get()) == 0) {
        tendpoints.clear();
        start_mad.reset();
        size_mad.reset();
        // b1->errcode  // MUST CHECK BEFORE WRITE TO VCF
        auto vd = get_aln_data (ap.get(), apit.get(), b1.get());

        if (vd.qpv.empty() || vd.tev.empty())
            throw std::runtime_error ("no data?");     // TODO placeholder

        // consensus region of supporting templates
        uint64_t lmosttc = std::numeric_limits<uint64_t>::max();
        uint64_t rmosttc = 0ULL;
        for (const auto &te : vd.tev) {
            if (te.min < lmosttc)
                lmosttc = te.min;
            if (te.max > rmosttc)
                rmosttc = te.max;
        }

        std::sort (begin (vd.qpv), end (vd.qpv));     // necessary for stats tests
        const auto qpos_mad = mad (vd.qpv).value_or (
            NAN
        );     // my span50 is better than mad because data is not necessarily unimodal, but it's interesting to compare
        const auto qpos_distrib_spans =
            std::vector<double>{0.50, 0.90}
            | std::views::transform ([&] (double pt) {
                  return min_span_containing (vd.qpv, pt).value_or (NAN);
              })
            | std::ranges::to<std::vector<uint64_t>>();

        // TODO
        // get the ordered mst tree
        // from template endpoints
        // then do an equivalent "min span containing"
        // on the edge lengths

        // gap-based stats
        auto qgaps = gaps (vd.qpv);
        auto tcdv = mst_cheb_dists (vd.tev);     // template chebyshev mst distances

        std::sort (begin (qgaps), end (qgaps));
        std::sort (begin (tcdv), end (tcdv));

        const auto qgaps_tail_jump = tail_jump (qgaps).value_or (NAN);
        const auto tcd_tail_jump   = tail_jump (tcdv).value_or (NAN);
        const auto tcd_distrib     = summary_stats (tcdv);

        // for the qpos I've got a pretty good set of descriptive statistics now
        // but they HAVE to be compared to the null model of uniform distribution
        // and the distribution of non-supporting/reference calls to be of value
        // I'm hopeful that one can say something like, "at this variant
        // the supporting qpos distribution is very non-random compared to the null model
        // but note that so is the reference qpos distribution and so on".
        // Also, n.b. that the original plan was just to extract templates
        // and I've done that! I could fold some!
        // TODO null model
        // TODO ref data fetched from the get_aln_data func as well as supporting
        // TODO better descriptive stats for the template distribution (see other comment)
        // TODO eyeball comparison to hp2 -> subset a vcf to DVF and ADF marked
        // TODO some folding of templates!
        // TODO I may actually want to report the min/max template size
        // as it's relevant for folding potentially
        std::cout << std::format (
            "{}\t{}\t{}\t{}\t{}\t"
            "{}\t{}\t{}\t{}\t{}\t"
            "{}\t{}\t{}\t{}\t{}",
            std::to_string (b1->rid),
            std::to_string (b1->pos),
            std::to_string (qpos_mad),     // just for comparison to my stats for now
            std::to_string (qpos_distrib_spans[0]),     // what span contains 25,
            std::to_string (qpos_distrib_spans[1]),     // 90, and
            std::to_string (
                vd.qpv.back() - vd.qpv.front()
            ),     // 100% of supporting read coordinates
            std::to_string (qgaps_tail_jump),
            std::to_string (vd.qpv.size()),     // n supporting reads
            std::to_string (tcd_distrib.q25),     // placeholder
            std::to_string (tcd_distrib.q75),
            std::to_string (tcd_tail_jump),
            std::to_string (lmosttc),     // template consensus span
            std::to_string (rmosttc),
            std::to_string (rmosttc - lmosttc),     // span size
            std::to_string (vd.tev.size())     // n supporting templates (qname deduplicated)
        ) << std::endl;
    };

    std::cout << "complete" << std::endl;
    return 0;
}
