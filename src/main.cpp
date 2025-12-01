// clang-format off
// TODO UPDATE, NOT USING MAD AND SO ON NOW
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

// NOTE: qpos clustering is equivalent to read endpoint clustering iff read lengths are ~all the same
// which they are for short read seq
// template clustering, via five number summary of the template endpoint chebyshev distance MST edges
// CRITICAL NOTE:
// distribution of supporting data must not be meaningfully
// different to a random sampling of the total data
// if nothing odd is going on
// TODO eyeball comparison to hp2 -> subset a vcf to DVF and ADF marked
// TODO some folding of templates!

#include <algorithm>
#include <cstdint>
#include <format>
#include <iostream>

#include <cxxopts.hpp>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>

#include "hts_ptr_t.hpp"
#include "pileup.hpp"
#include "stats.hpp"


template <class T>
std::string opt_to_str (
    std::optional<T>               opt,
    std::string_view               sentinel,
    std::function<std::string (T)> conv = [] (const T &a) {
        return std::to_string (a);
    }
) {
    return opt ? conv (*opt) : std::string (sentinel);
}

std::string rdbl2 (const double &a) {
    return std::format ("{:.2f}", a);
}
std::string rdbl4 (const double &a) {
    return std::format ("{:.4f}", a);
}


// TODO allow exclusion of variants by filter flags
// TODO allow user defined samflags for include/exclude
// TODO support single ended? LATER
// TODO options for calculating subset of data
// TODO --encode-vcf <COLNAME1,2,3> && --ovcf && --otsv ("-") for stdout
// TODO options for more vcf data (e.g. REF,ALT) in TSV (if using expos as "genome browser by numbers")
// TODO use kolmogorov complexity of ref
// TODO formalise why span50 better than MAD
// TODO optionally use positional data from NORMAL/BULK/SOMATIC
// to as background for simulation (if it's the same protocol which I don't know - if possible great!)
// TODO compare to uniform distribution (less valuable than background but possibly useful if e.g. not enough reads otherwise)
int main (
    int   argc,
    char *argv[]
) {
    namespace fs = std::filesystem;
    fs::path vcf_path;
    fs::path aln_path;
    fs::path ref_path;

    cxxopts::Options options (
        "expos",
        "EXtract POSitional data and statistics from alignment at "
        "VCF-specified pileups\n"
    );

    // clang-format off
    options.add_options() ("h,help", "Print usage")
        ("vcf", "VCF", cxxopts::value<fs::path>())
        ("aln", "Sample BAM", cxxopts::value<fs::path>())
        ("r,ref",
         "Alignment Reference Fasta for optionally adding template kolmogorov to statistics", // TODO
         cxxopts::value<fs::path>());
    // clang-format on

    options.parse_positional ({"vcf", "aln"});
    options.positional_help ("<VCF> <ALN>");

    try {
        auto parsedargs = options.parse (argc, argv);

        if (parsedargs.count ("help")) {
            std::cout << options.help() << std::endl;
            return 0;
        }

        if (!parsedargs.count ("vcf") || !parsedargs.count ("aln"))
            throw std::runtime_error (
                "All positional arguments must be provided"
            );

        vcf_path = parsedargs["vcf"].as<fs::path>();
        aln_path = parsedargs["aln"].as<fs::path>();

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

        if (parsedargs.count ("ref")) {
            ref_path = parsedargs["ref"].as<fs::path>();
            if (!fs::exists (ref_path)) {
                throw std::runtime_error (
                    "Reference fasta not found: " + ref_path.string()
                );
            }
            std::cerr << "Using ref: " << ref_path << std::endl;
        }


    } catch (const std::exception &e) {
        std::cerr << "Error parsing CLI options: " << e.what()
                  << "\n";
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
    std::cout << "CHROM\tPOS\tQPOS_SPAN50\tQPOS_SPAN90\tQPOS_"
                 "SPAN\tQPOS_TAIL_JUMP\tSPAN50_EFF2BG\tSPAN50_"
                 "PVAL\tNALT\tNTOTAL\tTEMPL_RAD50\tTEMPL_"
                 "RAD90\tTEMPL_RAD100\tTEMPL_TAIL_JUMP\tRAD50_"
                 "EFF2BG\tRAD50_PVAL\tTEMPL_"
                 "LMOST\tTEMPL_RMOST\tTEMPL_SPAN\tNTEMPL"
              << "\n";
    while (bcf_read (vp.get(), vph.get(), b1.get()) == 0) {
        // b1->errcode  // MUST CHECK BEFORE WRITE TO VCF
        auto  vard = get_aln_data (ap.get(), apit.get(), b1.get());
        auto &altd = vard.alt;
        std::vector<uint64_t> qpos_popv;
        qpos_popv.insert (
            qpos_popv.end(),
            begin (vard.alt.qpv),
            end (vard.alt.qpv)
        );
        qpos_popv.insert (
            qpos_popv.end(),
            begin (vard.other.qpv),
            end (vard.other.qpv)
        );
        std::vector<endpoints1D<uint64_t>> te_popv;
        te_popv.insert (
            te_popv.end(),
            begin (vard.alt.tev),
            end (vard.alt.tev)
        );
        te_popv.insert (
            te_popv.end(),
            begin (vard.other.tev),
            end (vard.other.tev)
        );

        if (altd.qpv.empty() || altd.tev.empty())
            throw std::runtime_error (
                "no data?"
            );     // TODO placeholder

        // STATS ON QUERY POSITIONS
        // total span of query position
        const auto [qpmin, qpmax] = std::minmax_element (
            begin (altd.qpv),
            end (altd.qpv)
        );
        const auto
            qpos_distrib_spans = std::vector<double>{0.5, 0.9}
                                 | std::views::transform (
                                     [&altd] (double pt) {
                                         return min_span_containing (
                                             altd.qpv,
                                             pt
                                         );
                                     }
                                 )
                                 | std::ranges::to<std::vector>();
        // simulate against span50
        stat_eval_s span50sim;
        if (qpos_distrib_spans[0]) {
            span50sim = sim_to_bg<uint64_t, uint64_t> (
                *qpos_distrib_spans[0],
                altd.qpv.size(),
                qpos_popv,
                [] (const decltype (qpos_popv) &v) {
                    const auto ret = min_span_containing (v, 0.5);
                    if (!ret) {
                        throw std::runtime_error (
                            "bad call to min_span_containing, "
                            "program malformed"
                        );
                    }
                    return *ret;
                },
                [&qpos_distrib_spans] (const auto s) {
                    return s <= *qpos_distrib_spans[0];
                }
            );
        }
        // gap-based stats
        const auto qgaps           = gaps (altd.qpv);
        const auto qgaps_tail_jump = tail_jump (qgaps);

        // STATS ON TEMPLATE ENDPOINTS
        // consensus region of supporting templates
        uint64_t lmosttc = std::numeric_limits<uint64_t>::max();
        uint64_t rmosttc = 0ULL;
        for (const auto &te : altd.tev) {
            if (te.min < lmosttc)
                lmosttc = te.min;
            if (te.max > rmosttc)
                rmosttc = te.max;
        }
        const auto
            te_distrib_spans = std::vector<double>{0.5, 0.9, 1}
                               | std::views::transform ([&altd] (
                                                            double pt
                                                        ) {
                                     // does not require sorting
                                     return min_cheb_radius_containing (
                                         altd.tev,
                                         pt
                                     );
                                 })
                               | std::ranges::to<std::vector>();
        stat_eval_s rad50sim;
        if (te_distrib_spans[0]) {
            rad50sim = sim_to_bg<uint64_t, endpoints1D<uint64_t>> (
                *te_distrib_spans[0],
                altd.tev.size(),
                te_popv,
                [] (const decltype (te_popv) &v) {
                    const auto ret = min_cheb_radius_containing (
                        v,
                        0.5
                    );
                    if (!ret) {
                        throw std::runtime_error (
                            "bad call to min_cheb_radius_containing, "
                            "program malformed"
                        );
                    }
                    return *ret;
                },
                [&te_distrib_spans] (const auto s) {
                    return s <= *te_distrib_spans[0];
                }
            );
        }

        const auto tcdv = mst_cheb_dists (
            altd.tev
        );     // template chebyshev mst distances
        const auto tcd_tail_jump = tail_jump (tcdv);

        // report informative set of descriptive statistics
        // clang-format off
        std::cout << std::format (
            "{}\t{}\t{}\t{}\t{}\t"
            "{}\t{}\t{}\t{}\t{}\t"
            "{}\t{}\t{}\t{}\t{}"
            "\t{}\t{}\t{}\t{}\t{}",
            b1->rid,
            b1->pos,
            // TODO reorder all sizes to the end
            // report nsupporting templates + n ref templates
            opt_to_str (te_distrib_spans[0], "NA"),     // what span contains 50,
            opt_to_str (qpos_distrib_spans[1], "NA"),     // 90, and
            std::to_string (*qpmax - *qpmin),     // 100% of supporting read coordinates
            opt_to_str<double> (qgaps_tail_jump, "NA", rdbl2),
            opt_to_str<double> (span50sim.eff_sz, "NA", rdbl2),
            opt_to_str<double> (span50sim.pval, "NA", rdbl4),
            std::to_string (altd.qpv.size()),     // n supporting reads
            std::to_string (qpos_popv.size()),
            opt_to_str (te_distrib_spans[0], "NA"),     // what chebyshev radius contains 50,
            opt_to_str (te_distrib_spans[1], "NA"),     // 90, and
            opt_to_str (te_distrib_spans[2], "NA"),     // 100% of supporting template coordinates
            opt_to_str<double> (tcd_tail_jump, "NA", rdbl2),
            opt_to_str<double> (rad50sim.eff_sz, "NA", rdbl2),
            opt_to_str<double> (rad50sim.pval, "NA", rdbl4),
            std::to_string (lmosttc),     // template consensus span
            std::to_string (rmosttc),
            std::to_string (rmosttc - lmosttc),     // span size
            std::to_string (altd.tev.size())     // n supporting templates (qname deduplicated)
        ) << "\n";
        // clang-format on
    };

    std::cerr << "complete" << std::endl;
    return 0;
}
