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
#include <htslib/faidx.h>
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

    fai_upt rp{nullptr, fai_destroy};
    if (!ref_path.empty()) {
        auto _fin = fai_load (ref_path.string().c_str());
        if (_fin == NULL) {
            std::cerr << std::format (
                "Could not read reference fasta at {}",
                ref_path.string()
            );
        } else {
            rp.reset (_fin);
        }
    }

    htsFile_upt ap{std::move (_ain), hts_close};
    hts_idx_upt apit{std::move (_aixin), hts_idx_destroy};

    htsFile_upt vp{std::move (_vin), hts_close};
    bcf_hdr_upt vph{std::move (_vh), bcf_hdr_destroy};

    bcf1_upt b1{bcf_init(), bcf_destroy};


    // loop vars
    std::cout
        << "CHROM\tPOS\tQPOS_SPAN50\tQPOS_SPAN90\tQPOS_"
           "SPAN\tQPOS_TAIL_JUMP\tQPOS_SPAN50_EFF2BG\tQPOS_SPAN50_"
           "PVAL\tTEMPL_SPAN50\tTEMPL_"
           "SPAN90\tTEMPL_SPAN100\tTEMPL_TAIL_JUMP\t"
           "TEMPL_SPAN50_"
           "EFF2BG\tTEMPL_SPAN50_PVAL\tTEMPL_"
           "LMOST\tTEMPL_RMOST\tTEMPL_CONSENSUS_LEN\tCONSENSUS_"
           "CMPLX\tNALT_"
           "READS\tNTOTAL_READS\tNALT_TEMPL\tNTOTAL_TEMPL"
        << "\n";
    while (bcf_read (vp.get(), vph.get(), b1.get()) == 0) {
        // NOTE b1->errcode  // MUST CHECK BEFORE WRITE TO VCF
        if (b1->n_allele != 2) {
            std::cerr << std::format (
                "Variant {} has more than two (REF,ALT) alleles. "
                "Unnormalised variant calls are not supported. "
                "Skipping.",
                b1->d.id
            ) << std::endl;
            continue;
        }
        auto mtype = bcf_has_variant_type (
            b1.get(),
            1,
            VCF_DEL | VCF_INS | VCF_SNP | VCF_MNP
        );
        switch (mtype) {
            case (VCF_DEL):
            case (VCF_INS):
            case (VCF_SNP):
            case (VCF_MNP):
                break;
            default:
                std::cerr << std::format (
                    "Variant {} has an unsupported/complex mutation, "
                    "skipping.",
                    b1->d.id
                ) << std::endl;
                continue;
        }
        auto vard = get_aln_data (
            ap.get(),
            apit.get(),
            b1.get(),
            mtype
        );
        auto                 &altd = vard.alt;
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
        std::vector<line_seg<uint64_t>> te_popv;
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

        if (altd.qpv.empty() || altd.tev.empty()) {
            std::cerr << std::format (
                "no supporting reads found for variant {}, "
                "skipping.",
                b1->d.id
            ) << std::endl;
            continue;
        }

        // STATS ON QUERY POSITIONS
        // total span of query position
        const auto qpos_span50 = min_span_containing (altd.qpv, 0.5);
        const auto qpos_span90 = min_span_containing (altd.qpv, 0.9);
        const auto [qpmin, qpmax] = std::minmax_element (
            begin (altd.qpv),
            end (altd.qpv)
        );
        const auto qpos_span = *qpmax - *qpmin;

        // simulate against span50
        stat_eval_s span50sim;
        if (qpos_span50) {
            span50sim = sim_to_bg<uint64_t, uint64_t> (
                *qpos_span50,
                altd.qpv.size(),
                qpos_popv,
                [] (const decltype (qpos_popv) &v) {
                    const auto ret = min_span_containing (v, 0.5);
                    assert (ret);
                    return *ret;
                },
                [&qpos_span50] (const auto s) {
                    return s <= *qpos_span50;
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
            if (te.lmost < lmosttc)
                lmosttc = te.lmost;
            if (te.rmost > rmosttc)
                rmosttc = te.rmost;
        }

        const auto te_pwds = PairMatrix<uint64_t>::from_sample (
            altd.tev,
            ucheb<uint64_t>
        );
        const auto te_span50 = te_pwds.min_span_containing (0.5);
        const auto te_span90 = te_pwds.min_span_containing (0.9);
        const auto te_span   = te_pwds.min_span_containing (1);

        // simulate against span50
        stat_eval_s te_span50sim;
        if (te_span50) {
            te_span50sim = sim_to_bg<uint64_t, line_seg<uint64_t>> (
                *te_span50,
                altd.tev.size(),
                te_popv,
                [] (const decltype (te_popv) &v) {
                    const auto pwds = PairMatrix<uint64_t>::from_sample (
                        v,
                        ucheb<uint64_t>
                    );
                    const auto ret = pwds.min_span_containing (0.5);
                    assert (ret);
                    return *ret;
                },
                [&te_span50] (const auto s) { return s <= *te_span50; }
            );
        }

        const auto te_tail_jump = tail_jump (te_pwds.get1D());

        // const auto kolmc = nk_lz76()
        std::optional<double> kolmc;
        if (rp) {
            // TODO should really check if it's in bam header, and that this is the correct reference
            auto rid_name = bcf_hdr_id2name (vph.get(), b1->rid);
            if (rid_name == NULL) {
                throw std::runtime_error (
                    std::format (
                        "Could not find rid {} in VCF header - VCF "
                        "misformatted?",
                        b1->rid
                    )
                );
            }
            std::string refs = fai_autofetch (
                rp.get(),
                rid_name,
                lmosttc,
                rmosttc
            );
            // TODO warn if all N
            kolmc.emplace (nk_lz76 (refs));
        }

        // report informative set of descriptive statistics
        // clang-format off
        std::cout << std::format (
            "{}\t{}\t{}\t{}\t{}\t"
            "{}\t{}\t{}\t{}\t{}\t"
            "{}\t{}\t{}\t{}\t{}"
            "\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            b1->rid,
            b1->pos,
            opt_to_str (qpos_span50, "NA"),     // what span contains 50,
            opt_to_str (qpos_span90, "NA"),     // 90, and
            std::to_string (qpos_span),     // 100% of supporting read coordinates
            opt_to_str<double> (qgaps_tail_jump, "NA", rdbl2),
            opt_to_str<double> (span50sim.eff_sz, "NA", rdbl2),
            opt_to_str<double> (span50sim.pval, "NA", rdbl4),
            opt_to_str (te_span50, "NA"),     // what chebyshev radius contains 50,
            opt_to_str (te_span90, "NA"),     // 90, and
            opt_to_str (te_span, "NA"),     // 100% of supporting template coordinates
            opt_to_str<double> (te_tail_jump, "NA", rdbl2),
            opt_to_str<double> (te_span50sim.eff_sz, "NA", rdbl2),
            opt_to_str<double> (te_span50sim.pval, "NA", rdbl4),
            std::to_string (lmosttc),     // template consensus span
            std::to_string (rmosttc),
            std::to_string (rmosttc - lmosttc),     // span size
            opt_to_str<double> (kolmc, "NA", rdbl2),
            std::to_string (altd.qpv.size()),     // n supporting reads
            std::to_string (qpos_popv.size()),
            std::to_string (altd.tev.size()),     // n supporting templates (qname deduplicated)
            std::to_string (te_popv.size())
        ) << "\n";
        // clang-format on
    };

    std::cerr << "complete" << std::endl;
    return 0;
}
