// clang-format off
// GET TEMPLATES AND DISTRIBUTION OF TEMPLATES SUPPORTING A VARIANT
// for variant, report:
// - template start total range and MAD [Median Absolute Distribution]*
// - template size total range and MAD*
// - leftmost template start
// - rightmost template start
// - leftmost template end
// - rightmost template end
// (optionally) merge templates within +-n bp of eachother
// write (merged) templates with representative qnames + distribution/s
// into VCF
// across all variants tested, report to file:
// - distribution of template start range
// - distribution of template size range
// optional reports to separate file:
// - template details per qname
// *If MAD ~= range, broad distribution, if range >> MAD, narrow distribution with outlier/s
// Using MAD because it is a robust spread statistic (i.e. not sensitive to outliers)
// and better able to deal with a low number of observations than IQR.
// ---
// The above in itself is probably an advancement in method for assessing variant
// legitimacy for low yield seq. If theres very little template variation (compared to the average)
// then almost certainly it's an artefact. That is, you can filter on this data alone
// ---
// WITH BCFTOOLS CONSENSUS
// fetch template region from reference.
// attempt to assemble sample template from ref base,
// i.e. with mutations if any, by merging good reads
// WITH RNAlfold
// evaluate foldiness of recovered templates
// clang-format on

#include <iostream>

#include <cxxopts.hpp>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>

#include "hts_ptr_t.hpp"
#include "pileup.hpp"

// TODO allow exclusion of variants by filter flags
// TODO allow user defined samflags for include/exclude
// TODO support single ended?
// TODO options to output to VCF and/or TSV (VCF to stdout)
int main (
    int argc,
    char *argv[]
) {
    namespace fs = std::filesystem;
    fs::path ref_path;
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

    options.parse_positional ({"vcf", "aln"});
    options.positional_help ("<VCF> <ALN>");
    options.show_positional_help();

    try {
        auto result = options.parse (argc, argv);

        if (result.count ("help")) {
            std::cout << options.help() << std::endl;
            return 0;
        }

        if (!result.count ("vcf") || !result.count ("aln"))
            throw std::runtime_error (
                "All positional arguments must be provided");

        vcf_path = result["vcf"].as<fs::path>();
        aln_path = result["aln"].as<fs::path>();

        if (!fs::exists (vcf_path)) {
            throw std::runtime_error ("VCF file not found: "
                                      + vcf_path.string());
        }

        std::cout << "Using VCF: " << vcf_path << std::endl;

        if (!fs::exists (aln_path)) {
            throw std::runtime_error ("Alignment file not found: "
                                      + aln_path.string());
        }

        std::cout << "Using aln: " << aln_path << std::endl;
    } catch (const std::exception &e) {
        std::cerr << "Error parsing CLI options: " << e.what()
                  << "\n";
        return 1;
    }
    // clang-format on

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
    htsFile_upt ap{std::move (_ain), hts_close};
    htsFile_upt vp{std::move (_vin), hts_close};
    hts_idx_upt apit{std::move (_aixin), hts_idx_destroy};

    auto _vh{bcf_hdr_read (vp.get())};
    bcf_hdr_upt vph{std::move (_vh), bcf_hdr_destroy};

    bcf1_upt b1{bcf_init(), bcf_destroy};

    template_vt var_tmpls;
    while (bcf_read (vp.get(), vph.get(), b1.get()) == 0) {
        // b1->errcode  // MUST CHECK BEFORE WRITE TO VCF
        var_tmpls = get_templates (ap.get(), apit.get(), b1.get());
        // TODO/NEXT
    };

    return 0;
}
