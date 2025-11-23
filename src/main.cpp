// clang-format off
// GET TEMPLATES AND DISTRIBUTION OF TEMPLATES SUPPORTING A VARIANT
// for variant, report:
// - template start total range and MAD [Median Absolute Distribution]*
// - template size total range and MAD*
// - leftmost template start
// - rightmost template start
// - leftmost template end
// - rightmost template end
// across all variants tested, report to file:
// - distribution of template start range
// - distribution of template start MAD
// - distribution of template size range
// - distribution of template size MAD
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

#include <algorithm>
#include <cmath>
#include <iostream>

#include <cxxopts.hpp>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>

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

        std::cout << "Using VCF: " << vcf_path << std::endl;

        if (!fs::exists (aln_path)) {
            throw std::runtime_error (
                "Alignment file not found: " + aln_path.string()
            );
        }

        std::cout << "Using aln: " << aln_path << std::endl;
    } catch (const std::exception &e) {
        std::cerr << "Error parsing CLI options: " << e.what() << "\n";
        return 1;
    }

    auto _ain{hts_open (aln_path.c_str(), "r")};
    if (_ain == NULL) {
        std::cout << std::format (
            "Could not open alignment file at {}",
            aln_path.string()
        ) << std::endl;
        return 1;
    }
    auto _aixin{sam_index_load (_ain, aln_path.c_str())};
    if (_aixin == NULL) {
        std::cout << std::format (
            "Coud not open index for alignment "
            "file. Searched for {}.bai",
            aln_path.c_str()
        ) << std::endl;
    }

    auto _vin{hts_open (vcf_path.c_str(), "r")};
    if (_vin == NULL) {
        std::cout << std::format (
            "Could not open VCF file at {}",
            vcf_path.string()
        ) << std::endl;
        return 1;
    }
    auto _vh{bcf_hdr_read (_vin)};
    if (_vh == NULL) {
        std::cout << std::format (
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
    std::vector<read_template_s> var_tmpls;
    std::vector<int64_t> tstartv;     // for finding key endpoints,
                                      // and stats
    std::vector<size_t>  tsizev;      // key endpoints, and stats
    std::vector<int64_t> tendv;       // key endpoints
    std::array<int64_t, 4> key_endpoints;     // leftmost start, rmost start, lmost end, rmost end
    std::optional<double> start_mad, size_mad;
    while (bcf_read (vp.get(), vph.get(), b1.get()) == 0) {
        tstartv.clear();
        tendv.clear();
        tsizev.clear();
        start_mad.reset();
        size_mad.reset();
        // b1->errcode  // MUST CHECK BEFORE WRITE TO VCF
        var_tmpls = get_templates (ap.get(), apit.get(), b1.get());
        // if necessary, performance increase possible by calculating online during this loop
        if (var_tmpls.empty())
            throw std::runtime_error ("no templates?");
        for (const auto &t : var_tmpls) {
            tstartv.push_back (t.start);
            tsizev.push_back (t.len);
            tendv.push_back (t.end);
        }
        const auto [plmosts, prmosts] = std::minmax_element (
            begin (tstartv),
            end (tstartv)
        );
        const auto [plmoste, prmoste] = std::minmax_element (
            begin (tendv),
            end (tendv)
        );
        const auto [psize_min, psize_max] = std::minmax_element (
            begin (tsizev),
            end (tsizev)
        );
        key_endpoints = {*plmosts, *prmosts, *plmoste, *prmoste};
        start_mad     = mad (tstartv);
        size_mad      = mad (tsizev);
        // placeholder
        // NOTE multinomial distribution DOES matter because two tight clusters
        // could mean multiple duplicate groups from templates which both have
        // TODO clustering of template start
        // TODO clustering of qpos of variant (ADF)
        // (the same) issues.
        std::cout << std::format (
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            std::to_string (b1->rid),
            std::to_string (b1->pos),
            std::to_string (*prmosts - *plmosts),
            std::to_string (start_mad.value_or (NAN)),
            std::to_string (*psize_max - *psize_min),
            std::to_string (size_mad.value_or (NAN)),
            std::to_string (*plmosts),
            std::to_string (*prmoste),
            std::to_string (var_tmpls.size())     // n supporting templates
        ) << std::endl;
    };

    return 0;
}
