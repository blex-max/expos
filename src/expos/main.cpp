// clang-format off
// NOTE --
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

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <format>
#include <fstream>
#include <functional>
#include <htslib/faidx.h>
#include <iostream>

#include <cxxopts.hpp>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <random>
#include <stdexcept>
#include <unordered_map>

#include "hts_ptr_t.hpp"
#include "lib-stats/monte-carlo.hpp"
#include "pileup.hpp"
#include "lib-stats/string.hpp"

constexpr std::string PROG_NAME = "expos";
constexpr std::string VERSION   = "0.0.0";
struct field_s {
    std::string name;
    std::string desc;
    int         type, nrec;
};
// TODO sub in search radius
const std::unordered_map<std::string, field_s> FIELD_INF{
    {"QRK",
     {"QRK",
      "Array detailing Ripley's K for mutant query position, "
      "and monte-carlo simulation results:"
      "[0]calculated statistic;"
      "[1]log2 ratio effect size from comparisons to simulation against all reads;"
      "[2]two-sided P-value from comparisons to simulation against all reads;"
      "[3]log2 ratio effect size from comparisons to simulation against uniform distribution;"
      "[4]two-sided P-value from comparisons to simulation against uniform distribution",
      BCF_HT_REAL,
      5}},
    {"TRK",
     {"TRK",
      "Array detailing Ripley's K for endpoints of supporting templates, "
      "and monte-carlo simulation results:"
      "[0]calculated statistic;"
      "[1]log2 ratio effect size from comparisons to simulation against all reads;"
      "[2]two-sided P-value from comparisons to simluation against all reads.",
      BCF_HT_REAL,
      3}},
    {"RCMPLX",
     {"RCMPLX",
      "Array detailing complexity of the region of supporting templates: "
      "[0]Mean 100-base window complexity (Lempel-Ziv estimated entropy rate), scaled by x100"
      "[1]Mean length of the top 10% longest runs of repeated motifs of period 1-5"
      "[2]Longest run of a repeated motif of period 1-5",
      BCF_HT_REAL,
      3}},
    {"MLAS",
     {"MLAS",
      "Array of median read-length normalised alignment scores:"
      "[0]of reads supporting variant,"
      "[1]of all queried reads covering the variant location in the sample alignment",
      BCF_HT_REAL,
      2}}
};


// various helpers
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


// TODO consider multiple and adjustable t for clustering assessment -- if so encode t in field name like QRK-5
// TODO sim str run stats against reshuffled region
// TODO consider MFE/n as report for secondary structure propensity
// TODO add record of command to VCF!
// TODO fraction of supporting reads with soft clipping, eff sz, pval, and median number of clipped bases
// in reads with soft clipping (CLPM)
// TODO calculate ref statistics even if no supporting reads
// TODO options for more vcf data (e.g. REF,ALT) in TSV (if using expos as "genome browser by numbers")
// TODO assessment of cigar complexity/edit distance to ref of supporting reads
// compared to total population accounting for variant.
int main (
    int   argc,
    char *argv[]
) {
    namespace fs = std::filesystem;
    // TODO verify paths
    fs::path                 vcf_path;
    fs::path                 aln_path;
    fs::path                 norm_path;
    fs::path                 ref_path;
    fs::path                 otsv_path;
    std::vector<std::string> flt_inc;
    std::vector<std::string> flt_exc;
    uint32_t                 seed = 24601;
    size_t                   exp_read_len = 150;
    int                      flag_inc = 3;
    int                      flag_exc = 3852;
    bool                     no_gz = false;
    bool                     normal_only = false;
    // std::vector<std::string> wfields;

    // clang-format off
    cxxopts::Options options (
        "expos",
        "\n"
        "EXtract POSitional data and statistics from alignment at\n"
        "VCF variant sites, and encode them as INFO fields to VCF.\n"
        "Requires the presence of .(b/cr)ai indexes of the same name\n"
        "as the relevant alignment. Annotated VCF to stdout. See\n"
        "README or output VCF header for descriptions of fields\n"
        "added.\n"
    );

    options.add_options() ("h,help", "Print usage")
        // POSITIONAL
        ("vcf", "VCF", cxxopts::value<fs::path>())
        ("aln", "Sample BAM", cxxopts::value<fs::path>())

        // OPTS
        ("l,expected-read-len",
         "Sequencing read length. Default: 150",
         cxxopts::value<size_t>())
        ("r,ref",
         "Alignment Reference Fasta for optionally adding reference complexity to statistics.",
         cxxopts::value<fs::path>())
        ("n,normal",
         "Alignment for use as additional background data for simulation",
         cxxopts::value<fs::path>())
        ("normal-only",
         "Use only reads from the provided normal as background data, excluding non-supporting reads from the sample")

        ("i,include",
         "Only operate on VCF records with this value present in FILTER. e.g. -i PASS. May be passed multiple times.",
         cxxopts::value<std::vector<std::string>>()) // multiple allowed
        ("e,exclude",
         "Only operate on VCF records without this value present in FILTER. May be passed multiple times.",
         cxxopts::value<std::vector<std::string>>()) // multiple allowed
        // ("w,write",
        //  "Write specified field to output VCF. May be passed multiple times.",
        //  cxxopts::value<std::vector<std::string>>()->default_value("ALL"))

        ("t,tsv",
         "Write a tsv of extended statistics to file specified.",
         cxxopts::value<fs::path>())
        ("u,uncompressed", "output uncompressed VCF")
        ("seed",
        "Set random seed. Default: 24601",
        cxxopts::value<uint32_t>());
    // clang-format on

    options.parse_positional ({"vcf", "aln"});
    options.positional_help ("<VCF/BCF (- for stdin)> <ALN.(b/cr)am>");

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

        if (vcf_path.string() != "-" && !fs::exists (vcf_path)) {
            throw std::runtime_error (
                "VCF file not found: " + vcf_path.string()
            );
        }

        if (!fs::exists (aln_path)) {
            throw std::runtime_error (
                "Alignment file not found: " + aln_path.string()
            );
        }

        std::cerr << "Using VCF: " << vcf_path << std::endl;
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

        if (parsedargs.count ("expected-read-len")) {
            exp_read_len = parsedargs["expected-read-len"].as<size_t>();
        } else {
            std::cerr <<
            std::format("Read length not provided, assuming {}", exp_read_len)
            << std::endl;
        }

        if (parsedargs.count ("include")) {
            flt_inc = parsedargs["include"].as<std::vector<std::string>>();
        }
        if (parsedargs.count ("exclude")) {
            flt_exc = parsedargs["exclude"].as<std::vector<std::string>>();
        }
        if (parsedargs.count ("seed")) {
            seed = parsedargs["seed"].as<uint32_t>();
        }

        if (parsedargs.count ("tsv")) {
            otsv_path = parsedargs["tsv"].as<fs::path>();
        }

        if (parsedargs.count ("normal")) {
            norm_path = parsedargs["normal"].as<fs::path>();
            std::cerr << "Using normal: " << norm_path << std::endl;
        }
        if (parsedargs.count ("normal-only")) {
            if (norm_path.empty())
                throw std::runtime_error("a normal must be provided if normal-only is set.");
            std::cerr << "Using only normal data as background for simulation" << std::endl;
            normal_only = true;
        }

        if (parsedargs.count ("uncompressed")) {
            no_gz = true;
        }

    } catch (const std::exception &e) {
        std::cerr << "Error parsing CLI options: " << e.what()
                  << "\nTry expos --help.";
        return EXIT_FAILURE;
    }

    // inputs
    auto _ain{hts_open (aln_path.c_str(), "r")};
    if (_ain == NULL) {
        std::cerr << std::format (
            "Could not open alignment file at {}",
            aln_path.string()
        ) << std::endl;
        return EXIT_FAILURE;
    }
    auto _aixin{sam_index_load (_ain, aln_path.c_str())};
    if (_aixin == NULL) {
        std::cerr << std::format (
            "Coud not open index for alignment "
            "file. Searched for {}.bai",
            aln_path.c_str()
        ) << std::endl;
        return EXIT_FAILURE;
    }
    htsFile_upt alnfh{std::move (_ain), hts_close};
    hts_idx_upt aln_idx{std::move (_aixin), hts_idx_destroy};

    auto _vin{hts_open (vcf_path.c_str(), "r")};
    if (_vin == NULL) {
        std::cerr << std::format (
            "Could not open VCF file at {}",
            vcf_path.string()
        ) << std::endl;
        return EXIT_FAILURE;
    }
    auto _vh{bcf_hdr_read (_vin)};
    if (_vh == NULL) {
        std::cerr << std::format (
            "Could not read header of VCF file at {}",
            vcf_path.string()
        ) << std::endl;
        return EXIT_FAILURE;
    }
    htsFile_upt vcffh{std::move (_vin), hts_close};
    bcf_hdr_upt vcf_hdr{std::move (_vh), bcf_hdr_destroy};

    fai_upt reffh{nullptr, fai_destroy};
    if (!ref_path.empty()) {
        auto _fin = fai_load (ref_path.string().c_str());
        if (_fin == NULL) {
            std::cerr << std::format (
                "Could not read reference fasta at {}",
                ref_path.string()
            ) << std::endl;
            return EXIT_FAILURE;
        } else {
            reffh.reset (_fin);
        }
    }

    std::optional<std::pair<htsFile_upt, hts_idx_upt>> norm;
    if (!norm_path.empty()) {
        auto _nin{hts_open (norm_path.c_str(), "r")};
        if (_nin == NULL) {
            std::cerr << std::format (
                "Could not open alignment file at {}",
                norm_path.string()
            ) << std::endl;
            return EXIT_FAILURE;
        }
        auto _nixin{sam_index_load (_nin, norm_path.c_str())};
        if (_nixin == NULL) {
            std::cerr << std::format (
                "Coud not open index for alignment "
                "file. Searched for {}.bai",
                norm_path.c_str()
            ) << std::endl;
            return EXIT_FAILURE;
        }
        norm.emplace (
            htsFile_upt{std::move (_nin), hts_close},
            hts_idx_upt{std::move (_nixin), hts_idx_destroy}
        );
    }

    // outputs
    htsFile_upt ovcf{
        hts_open ("-", (no_gz ? "w" : "wz")),
        hts_close
    };     // stdout
    bcf_hdr_upt ohdr{bcf_hdr_dup (vcf_hdr.get()), bcf_hdr_destroy};

    // ADD LINES TO HDR
    constexpr auto make_info_line = [] (field_s i) {
        std::string t2s[4];
        t2s[BCF_HT_INT]  = "Integer";
        t2s[BCF_HT_REAL] = "Float";
        return std::format (
            "##INFO=<ID={},Number={},Type={},Description=\"{}\","
            "Source="
            "\"{}\",Version=\"{}\">",
            i.name,
            i.nrec,
            t2s[i.type],
            i.desc,
            PROG_NAME,
            VERSION
        );
    };
    for (const auto &l : FIELD_INF) {
        const auto &i = l.second;
        if (bcf_hdr_append (ohdr.get(), make_info_line (i).c_str())
            != 0) {
            throw std::runtime_error (
                std::format (
                    "failed to append to hdr for field {}",
                    i.name
                )
            );
        }
    }

    if (bcf_hdr_sync (ohdr.get()) < 0) {
        throw std::runtime_error ("failed to sync hdr");     // TODO
    }
    if (bcf_hdr_write (ovcf.get(), ohdr.get()) < 0) {
        throw std::runtime_error ("failed to write header");     // TODO
    };


    // optional tsv output
    std::optional<std::ofstream> otsv;
    if (!otsv_path.empty()) {
        otsv.emplace (otsv_path);
        // TODO -- pre-header comments explaining each field
        // NOTE -- TSV header should always stay the same.
        // Columns not calculated should simply be NA or other indicator
        *otsv << "CHROM\t"
                 "POS\t"
                 "MLAS\t"
                 "MLAS_EFFSZ\t"
                 "MLAS_PVAL\t"
                 "QPOS_RIPLEY\t"
                 "QPOS_RIPLEY_EFFSZ_TO_BACKGROUND\t"
                 "QPOS_RIPLEY_PVAL_TO_BACKGROUND\t"
                 "QPOS_RIPLEY_EFFSZ_TO_UNIFORM\t"
                 "QPOS_RIPLEY_PVAL_TO_UNIFORM\t"
                 "TEMPL_RIPLEY\t"
                 "TEMPL_RIPLEY_EFFSZ\t"
                 "TEMPL_RIPLEY_PVAL\t"
                 "CONSENSUS_CMPLXx100\t"
                 "REF_TOP10_STR_RUN_MEAN\t"
                 "REF_TOP_STR_RUN_LEN\t"
                 "LMOST_TEMPLATE_START\t"
                 "RMOST_TEMPLATE_END\t"
                 "NALT_READS\t"
                 "NTOTAL_READS"
              << "\n";
    }

    std::mt19937 rng{};
    rng.seed(seed);
    bool         firsti = true;
    bcf1_upt     b1{bcf_init(), bcf_destroy};
    while (bcf_read (vcffh.get(), vcf_hdr.get(), b1.get()) == 0) {
        if (firsti) {
            for (const auto &f : flt_inc) {
                if (bcf_has_filter (
                        vcf_hdr.get(),
                        b1.get(),
                        const_cast<char *> (f.c_str())
                    )
                    < 0) {
                    std::cerr << std::format (
                        "Unknown --include filter {} not present in "
                        "VCF",
                        f
                    ) << std::endl;     // unrecoverable
                    return EXIT_FAILURE;
                }
            }
            std::vector<std::string> tmp_ex;
            for (size_t i = 0; i < flt_exc.size(); ++i) {
                const auto &f = flt_exc[i];
                if (bcf_has_filter (
                        vcf_hdr.get(),
                        b1.get(),
                        const_cast<char *> (f.c_str())
                    )
                    < 0) {
                    std::cerr << std::format (
                        "Warning: Unknown --exclude filter {} not "
                        "present in VCF, "
                        "ignoring",
                        f
                    ) << std::endl;
                } else {
                    tmp_ex.push_back (flt_exc[i]);
                }
            }
            flt_exc = tmp_ex;
            firsti = false;
        }

        const auto iflt = has_filters (
            vcf_hdr.get(),
            b1.get(),
            flt_inc
        );
        if (std::any_of (begin (iflt), end (iflt), [] (const auto a) {
                return !a;
            })) {
            if (bcf_write (ovcf.get(), ohdr.get(), b1.get()) != 0) {
                std::cerr << std::format (
                    "failed to write record to output VCF"
                ) << std::endl;
                return EXIT_FAILURE;
            };
            continue;
        }
        const auto eflt = has_filters (vcf_hdr.get(), b1.get(), flt_exc);
        if (std::any_of (begin (eflt), end (eflt), [] (const auto a) {
                return a;
            })) {
            if (bcf_write (ovcf.get(), ohdr.get(), b1.get()) != 0) {
                std::cerr << std::format (
                    "failed to write record to output VCF"
                ) << std::endl;
                return EXIT_FAILURE;
            };
            continue;
        }

        // NOTE b1->errcode  // MUST CHECK BEFORE WRITE TO VCF
        if (b1->n_allele != 2) {
            std::cerr << std::format (
                "Variant {} {} {} has more than two (REF,ALT) "
                "alleles. "
                "Unnormalised variant calls are not supported. "
                "Skipping.",
                b1->rid,     // TODO convert rid to user facing
                b1->pos,     // ditto
                b1->d.id
            ) << std::endl;
            if (bcf_write (ovcf.get(), ohdr.get(), b1.get()) != 0) {
                std::cerr << std::format (
                    "failed to write record to output VCF"
                ) << std::endl;
                return EXIT_FAILURE;
            };
            continue;
        }
        auto mtype = bcf_has_variant_type (
            b1.get(),
            1,
            VCF_DEL | VCF_INS | VCF_SNP | VCF_MNP
        );
        // one and only one of
        switch (mtype) {
            case (VCF_DEL):
            case (VCF_INS):
            case (VCF_SNP):
            case (VCF_MNP):
                break;
            default:
                std::cerr << std::format (
                    "Variant {} has an unsupported/complex mutation, or could not be typed"
                    "skipping.",
                    b1->d.id
                ) << std::endl;
                if (bcf_write (ovcf.get(), ohdr.get(), b1.get())
                    != 0) {
                    std::cerr << std::format (
                        "failed to write record to output VCF"
                    ) << std::endl;
                    return EXIT_FAILURE;
                };
                continue;
        }

        auto [sample_supporting_pileup, sample_total_pileup] = pileup_partition_and_anaylse (
            alnfh.get(),
            aln_idx.get(),
            b1.get(),
            mtype,
            flag_inc,
            flag_exc
        );
        std::optional<PileupMetrics> normal_pileup;
        if (norm) {
            normal_pileup = pileup_analyse (
                norm->first.get(),
                norm->second.get(),
                b1.get(),
                flag_inc,
                flag_exc
            );
        }

        auto n_supporting_reads = sample_supporting_pileup.nreads;
        auto n_supporting_templates = sample_supporting_pileup.template_endpoints.size();

        if (n_supporting_reads == 0) {
            std::cerr << std::format (
                "no supporting reads found for variant {} {} {}, "
                "skipping.",
                b1->rid,     // TODO convert rid to user facing
                b1->pos,     // ditto
                b1->d.id
            ) << std::endl;
            if (bcf_write (ovcf.get(), ohdr.get(), b1.get()) != 0) {
                std::cerr << std::format (
                    "failed to write record to output VCF"
                ) << std::endl;
                return EXIT_FAILURE;
            };
            continue;
        }

        // TODO guard only relevant bits
        if (sample_supporting_pileup.query_position.size() == 0 || sample_supporting_pileup.template_endpoints.size() == 0) {
            std::cerr << std::format (
                "no supporting reads found with usable metrics for variant {} {} {}, "
                "skipping.",
                b1->rid,     // TODO convert rid to user facing
                b1->pos,     // ditto
                b1->d.id
            ) << std::endl;
            if (bcf_write (ovcf.get(), ohdr.get(), b1.get()) != 0) {
                std::cerr << std::format (
                    "failed to write record to output VCF"
                ) << std::endl;
                return EXIT_FAILURE;
            };
            continue;
        }
        if (normal_pileup && normal_pileup->nreads == 0) {
            std::cerr << std::format (
                "Warning: no reads covering variant location found in normal for variant {} {} {}",
                b1->rid,     // TODO convert rid to user facing
                b1->pos,     // ditto
                b1->d.id
            ) << std::endl;
        }

        // --- REF COMPLEXITY --- //
        // consensus region of supporting templates
        uint64_t lmosttc = std::numeric_limits<uint64_t>::max();
        uint64_t rmosttc = 0ULL;
        for (const auto &te : sample_supporting_pileup.template_endpoints) {
            if (te.lmost < lmosttc)
                lmosttc = te.lmost;
            if (te.rmost > rmosttc)
                rmosttc = te.rmost;
        }
        const auto span_length = rmosttc - lmosttc;

        // TODO should really check if it's in bam header not just the vcf,
        // and that this is the correct reference
        // NOTE not all needed to be fetched each loop
        auto rid_name = bcf_hdr_id2name (vcf_hdr.get(), b1->rid);
        std::optional<size_t> ref_entropy;
        std::optional<size_t> run_max;
        std::optional<double> run_mean;
        if (reffh) {
            if (rid_name == NULL) {
                std::cerr << std::format (
                    "Could not find reference ID {} in VCF header "
                    "- VCF "
                    "misformatted?",
                    b1->rid
                ) << std::endl;
            }
            std::string refs = fai_autofetch (
                reffh.get(),
                rid_name,
                lmosttc,
                rmosttc
            );
            // transform to ignore masks
            transform(begin(refs), end(refs), begin(refs), ::toupper);
            if (refs.find("N") == std::string::npos) {  // No Ns
                // lz complexity
                const size_t window_size = 100;
                double cmplx_sum = 0;
                size_t n_win = 0;
                for (;(n_win + window_size) <= refs.size(); ++n_win) {   // ++n_win == step of 1
                    cmplx_sum += string_stats::entropy_lz76(refs.substr(n_win, window_size));
                }
                const auto mean_window_entropy = static_cast<double> (cmplx_sum) / static_cast<double> (n_win);
                ref_entropy.emplace (
                    round (mean_window_entropy * 100)
                );     // x100 scaling factor

                auto str_runs = string_stats::periodic_rle(refs, 5);
                auto str_runs_n = str_runs.size();
                auto top10_i = static_cast<size_t> (floor(static_cast<double>(str_runs.size()) * 0.9));
                auto top10_n = str_runs_n - top10_i;
                std::nth_element(
                    begin(str_runs),
                    begin(str_runs) + static_cast<long> (top10_i), 
                    end(str_runs)
                );
                decltype(str_runs)::value_type str_run_max=0, top10_run_sum=0;
                for (size_t i = top10_i; i < str_runs_n; ++i) {
                    const auto e = str_runs[i];
                    if (e > str_run_max) {
                        str_run_max = e;
                    }
                    top10_run_sum += e;
                }
                run_max = str_run_max;
                run_mean = static_cast<double> (top10_run_sum) / static_cast<double> (top10_n);
            }
        }

        // --- CLUSTERING ANALYSIS - RIPLEY'S K --- //
        // get pairwise distances
        // simple 1D dist for query position
        constexpr auto dist_1D = [] (const auto &a,
                                 const auto &b) {
            return (a > b) ? (a - b) : (b - a);
        };
        const auto qpos_pwd = spatial::PairMatrix::from_sample (
            sample_supporting_pileup.query_position,
            dist_1D
        );     // empty if <2 samples

        // manhattan distance for template endpoints
        constexpr auto mannd = [] (const auto &a,
                                   const auto &b) {
            spatial::line_seg upper_pair{a.rmost, b.rmost};
            spatial::line_seg lower_pair{a.lmost, b.lmost};
            return upper_pair.diff() + lower_pair.diff();
        };
        const auto te_pwd = spatial::PairMatrix::from_sample (
            sample_supporting_pileup.template_endpoints,
            mannd
        );

        std::optional<double> qpos_rl;
        monte_carlo::stat_eval_s           qpos_rl_bgsim;  // compared to all reads
        monte_carlo::stat_eval_s           qpos_rl_unisim; // compared to expected distribution
        if (qpos_pwd) {
            const size_t search_radius = 5; // bases. Sensible values << read length.
            qpos_rl = 
                spatial::ripley_k(
                    *qpos_pwd,
                    search_radius,
                    static_cast<double>(qpos_pwd->dim())
                    / static_cast<double> (exp_read_len)
                );
            decltype(sample_supporting_pileup.query_position) qpos_popv;
            if (!normal_only) {
                if (normal_pileup) { // ADD NORMAL OBS
                    assert(normal_pileup);
                    qpos_popv.insert(
                        end(qpos_popv),
                        begin(sample_total_pileup.query_position),
                        end(sample_total_pileup.query_position)
                    );
                    qpos_popv.insert (
                        end(qpos_popv),
                        begin (normal_pileup->query_position),
                        end (normal_pileup->query_position)
                    );
                } else {  // total reads only
                    qpos_popv = sample_total_pileup.query_position;
                }
            } else {  // reads from normal only
                assert(normal_pileup);
                qpos_popv = normal_pileup->query_position;
            }

            auto stat_fn = [&dist_1D] (const auto &v) {
                const auto pwds = spatial::PairMatrix::from_sample (v, dist_1D);
                assert (pwds);
                return spatial::ripley_k(
                        *pwds,
                        search_radius,
                        static_cast<double>(v.size()) / 150.0  // read len
                    );
            };
            auto n_obs = sample_supporting_pileup.query_position.size();
            if (n_obs < 2) {
                qpos_rl_bgsim.err = "INSUFF_OBS";
            }
            else if (qpos_popv.size() < (n_obs * 2)) {
                // at a bare minimum, we want 2x more total samples than bg
                qpos_rl_bgsim.err = "INSUFF_BG";
            }  // TODO if less than e.g. 5x report low power?
            else {
                qpos_rl_bgsim = monte_carlo::sim_to_bg (
                    *qpos_rl,
                    [&qpos_popv, &rng, n_obs] () {
                        return monte_carlo::subsample_wo_replace(qpos_popv, n_obs, rng);
                    },
                    stat_fn,
                    monte_carlo::log2_effsz
                );
            }

            // simulate against uniform distribution
            // NOTE:
            // the same deviation from the null will always have the same
            // p value for a fixed sample size and read length, therefore
            // could precompute this and save wasted comuptation:
            // scale/normalise the statistic
            // S = (sample_size / read_length) * (Ripley)
            // Precompute one high-quality null CDF for S with a large sample size (>40)
            // Use that CDF for all datasets with moderate/large n.
            // for smaller n, use exactly stored CDF
            std::uniform_int_distribution<uint64_t> qpos_gen{0, exp_read_len};
            qpos_rl_unisim = monte_carlo::sim_to_bg(
                *qpos_rl,
                [&qpos_gen, &rng, n_obs] () {
                    std::vector<uint64_t> rand_qpos;
                    for (size_t i = 0; i < n_obs; ++i) {
                        rand_qpos.push_back(qpos_gen(rng));
                    }
                    return rand_qpos;
                },
                stat_fn,
                monte_carlo::log2_effsz
            );
        } else {
            qpos_rl_bgsim.err = "INSUFF_OBS";
            qpos_rl_unisim.err = "INSUFF_OBS";
        }

        std::optional<double> te_rl;
        monte_carlo::stat_eval_s           te_rl_sim;
        if (te_pwd) {
            const size_t search_radius = 6;
            te_rl = spatial::ripley_k(
                    *te_pwd,
                    search_radius,
                    static_cast<double>(te_pwd->dim()) / static_cast<double>(span_length)
            );
            decltype(sample_supporting_pileup.template_endpoints) te_popv;
            if (!normal_only) {
                if (normal_pileup) {
                    te_popv.insert(
                        end(te_popv),
                        begin(sample_total_pileup.template_endpoints),
                        end(sample_total_pileup.template_endpoints)
                    );
                    te_popv.insert (
                        te_popv.end(),
                        begin (normal_pileup->template_endpoints),
                        end (normal_pileup->template_endpoints)
                    );
                }
                else {
                    te_popv = sample_total_pileup.template_endpoints;
                }
            } else {
                te_popv = normal_pileup->template_endpoints;
            }

            if (n_supporting_templates < 2) {
                te_rl_sim.err = "INSUFF_OBS";
            }
            else if (te_popv.size() < n_supporting_templates * 2) {
                te_rl_sim.err = "INSUFF_BG";
            }
            else {
                te_rl_sim = monte_carlo::sim_to_bg (
                    *te_rl,
                    [&rng, &te_popv, n_supporting_templates] () {
                        return monte_carlo::subsample_wo_replace(te_popv, n_supporting_templates, rng);
                    },
                    [&mannd, span_length] (const auto &v) {
                        const auto pwds = spatial::PairMatrix::from_sample (
                            v,
                            mannd
                        );
                        assert (pwds);
                        return spatial::ripley_k(
                                *pwds,
                                search_radius,
                                static_cast<double>(pwds->dim()) / static_cast<double>(span_length));
                    },
                    monte_carlo::log2_effsz
                );
            }
        } else {
            te_rl_sim.err = "INSUFF_OBS";
        }

        // --- MEDIAN LENGTH-NORMALISED ALIGNMENT SCORE --- //
        const auto mlas_supporting = summary::percentile (
            sample_supporting_pileup.normalised_as,
            0.5
        );
        const auto mlas_total = summary::percentile (
            sample_total_pileup.normalised_as,
            0.5
        );

        // encode to vcf
        auto write_info = [&] (field_s i, const void *val) {
            if (bcf_update_info (
                    ohdr.get(),
                    b1.get(),
                    i.name.c_str(),
                    val,
                    i.nrec,
                    i.type
                )
                != 0) {
                std::cerr << std::format (
                    "failed to write {} as INFO field "
                    "in output VCF",
                    i.name
                ) << std::endl;
            }
        };
        // TODO should probably encode missingness into the vcf somehow... like an EXPOS_ERR info field
        if (qpos_rl) {
            float val[5]{
                static_cast<float> (*qpos_rl),
                static_cast<float> (qpos_rl_bgsim.eff_sz.value_or (0.0)),
                static_cast<float> (qpos_rl_bgsim.pval.value_or (1.0)),
                static_cast<float> (qpos_rl_unisim.eff_sz.value_or(1.0)),
                static_cast<float> (qpos_rl_unisim.pval.value_or (1.0)),
            };
            write_info (FIELD_INF.at ("QRK"), &val);
        }
        if (te_rl) {
            float val[3]{
                static_cast<float> (*te_rl),
                static_cast<float> (te_rl_sim.eff_sz.value_or (0.0)),
                static_cast<float> (te_rl_sim.pval.value_or (1.0))
            };
            write_info (FIELD_INF.at ("TRK"), &val);
        }
        if (ref_entropy) {
            float val[3] {
                static_cast<float> (*ref_entropy),
                static_cast<float> (*run_mean),
                static_cast<float> (*run_max)
            };
            write_info (FIELD_INF.at ("RCMPLX"), &val);
        }
        if (mlas_supporting) { // should still include total if no supporting, TODO
            float val[2]{
                static_cast<float> (*mlas_supporting),
                static_cast<float> (*mlas_total),
            };
            // TODO rounding
            write_info (FIELD_INF.at ("MLAS"), &val);
        }

        if (bcf_write (ovcf.get(), ohdr.get(), b1.get()) != 0) {
            std::cerr << std::format (
                "failed to write record to output VCF"
            ) << std::endl;
        };

        // report statistics to tsv
        // clang-format off
        if (otsv) {
            *otsv << std::format (
                "{}\t{}\t{}\t{}\t{}\t"
                "{}\t{}\t{}\t{}\t{}\t"
                "{}\t{}\t{}\t{}\t{}\t"
                "{}\t{}\t{}\t{}",
                rid_name,
                b1->pos + 1,
                opt_to_str<double> (mlas_supporting, "NA", rdbl2),
                opt_to_str<double> (mlas_total, "NA", rdbl2),
                opt_to_str<double> (qpos_rl, "NA", rdbl2),
                opt_to_str<double> (qpos_rl_bgsim.eff_sz, qpos_rl_bgsim.err, rdbl2),
                opt_to_str<double> (qpos_rl_bgsim.pval, qpos_rl_bgsim.err, rdbl4),
                opt_to_str<double> (qpos_rl_unisim.eff_sz, qpos_rl_unisim.err, rdbl2),
                opt_to_str<double> (qpos_rl_unisim.pval, qpos_rl_unisim.err, rdbl4),
                opt_to_str<double> (te_rl, "NA", rdbl2),
                opt_to_str<double> (te_rl_sim.eff_sz, te_rl_sim.err, rdbl2),
                opt_to_str<double> (te_rl_sim.pval, te_rl_sim.err, rdbl4),
                opt_to_str (ref_entropy, "NA"),
                opt_to_str (run_mean, "NA"),
                opt_to_str (run_max, "NA"),
                std::to_string(lmosttc),
                std::to_string(rmosttc),
                std::to_string (n_supporting_reads),
                std::to_string (sample_total_pileup.nreads)
            ) << "\n";
        }
        // clang-format on
    };

    std::cerr << "complete" << std::endl;
    return 0;
}
