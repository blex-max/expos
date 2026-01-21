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
// template clustering, via five number summary of the template endpoint chebyshev distance MST edges
// NOTE:
// distribution of supporting data must not be meaningfully
// different to a random sampling of the total data
// if nothing odd is going on

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
#include "pileup.hpp"
#include "stats.hpp"

constexpr std::string PROG_NAME = "expos";
constexpr std::string VERSION   = "0.0.0";
struct field_s {
    std::string name;
    std::string desc;
    int         type, nrec;
};
const std::unordered_map<std::string, field_s> FIELD_INF{
    {"MLAS",
     {"MLAS",
      "[0]Median read-Length-normalised Alignment Score (AS) of "
      "reads supporting variant;"
      "[1]delta (supporting - background) effect size and [2]two-sided P-value "
      "against "
      "background, from monte-carlo simulation",
      BCF_HT_REAL,
      3}},
    {"QM1NN",
     {"QM1NN",
      "[0]Median nearest neighbour distance of variant query "
      "position; [1]log2 ratio effect size and [2]two-sided P-value against "
      "background, from monte-carlo simulation",
      BCF_HT_REAL,
      3}},
    {"TM1NN",
     {"TM1NN",
      "[0]Median nearest neighbour distance of template endpoints "
      "from read pairs supporting variant; [1]log2 ratio effect size "
      "and [2]two-sided P-value against background, from "
      "monte-carlo simulation",
      BCF_HT_REAL,
      3}},
    {"RCMPLX",
     {"RCMPLX",
      "Complexity (Lempel-Ziv estimated entropy rate) of region spanned by supporting "
      "templates, scaled by x100",
      BCF_HT_INT,
      1}}     // NOTE -- kc is optional
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


// TODO add record of command to VCF!
// TODO fraction of supporting reads with soft clipping, eff sz, pval, and median number of clipped bases
// in reads with soft clipping (CLPM)
// TODO calculate average AS of all reads in sample region
// TODO calculate ref statistics even if no supporting reads
// TODO options for more vcf data (e.g. REF,ALT) in TSV (if using expos as "genome browser by numbers")
// NOTE uniform sim added for qpos!
// template endpoints is more complicated, they tend to show a right-skewed gaussian distribution
// around a target fragment size - TODO
// TODO allow option for only using normal as background - in progress
// 
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
    size_t                   read_len = 150;
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
        ("i,include",
         "Only operate on VCF records with this value present in FILTER. e.g. -i PASS. May be passed multiple times.",
         cxxopts::value<std::vector<std::string>>()) // multiple allowed
        ("e,exclude",
         "Only operate on VCF records without this value present in FILTER. May be passed multiple times.",
         cxxopts::value<std::vector<std::string>>()) // multiple allowed
        ("f,flag-include",
         "Only consider reads with these bits set in the SAM flag. Applies to both target and background alignment data. Default: 3",
         cxxopts::value<int>())
        ("F,flag-exclude",
         "Do not consider reads with these bits set in the SAM flag. Applies to both target and background alignment data. Default: 3852",
         cxxopts::value<int>())
        // ("w,write",
        //  "Write specified field to output VCF. May be passed multiple times.",
        //  cxxopts::value<std::vector<std::string>>()->default_value("ALL"))
        ("t,tsv",
         "Write a tsv of extended statistics to file specified.",
         cxxopts::value<fs::path>())
        ("n,normal",
         "Alignment for use as additional background data for simulation",
         cxxopts::value<fs::path>())
        ("normal-only",
         "Use only reads from the provided normal as background data, excluding non-supporting reads from the sample")
        ("r,ref",
         "Alignment Reference Fasta for optionally adding reference complexity to statistics.",
         cxxopts::value<fs::path>())
        ("seed",
        "Set random seed. Default: 24601",
        cxxopts::value<uint32_t>())
        ("u,uncompressed", "output uncompressed VCF");
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

        if (parsedargs.count ("include")) {
            flt_inc = parsedargs["include"].as<std::vector<std::string>>();
        }
        if (parsedargs.count ("exclude")) {
            flt_exc = parsedargs["exclude"].as<std::vector<std::string>>();
        }
        if (parsedargs.count ("flag-include")) {
            flag_inc = parsedargs["flag-include"].as<int>();
        }
        if (parsedargs.count ("flag-exclude")) {
            flag_exc = parsedargs["flag-exclude"].as<int>();
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
                 "QPOS_M1NN\t"
                 "QPOS_M1NN_EFFSZ_TO_BACKGROUND\t"
                 "QPOS_M1NN_PVAL_TO_BACKGROUND\t"
                 "QPOS_M1NN_EFFSZ_TO_UNIFORM\t"
                 "QPOS_M1NN_PVAL_TO_UNIFORM\t"
                 "TEMPL_M1NN\t"
                 "TEMPL_M1NN_EFFSZ\t"
                 "TEMPL_M1NN_PVAL\t"
                 "CONSENSUS_CMPLXx100\t"
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
        auto n_supporting_templates = sample_supporting_pileup.ntemplates;

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
        if (normal_pileup && normal_pileup->nreads == 0) {
            std::cerr << std::format (
                "Warning: no reads covering variant location foudn in normal for variant {} {} {}",
                b1->rid,     // TODO convert rid to user facing
                b1->pos,     // ditto
                b1->d.id
            ) << std::endl;
        }

        // --- CLUSTERING ANALYSIS (nearest neighbour) --- //
        // get pairwise distances
        // simple 1D dist for query position
        constexpr auto dist_1D = [] (const auto &a,
                                 const auto &b) {
            return (a > b) ? (a - b) : (b - a);
        };
        const auto qpos_pwd = PairMatrix::from_sample (
            sample_supporting_pileup.query_position,
            dist_1D
        );     // empty if <2 samples

        // manhattan distance for template endpoints
        constexpr auto mannd = [] (const auto &a,
                                   const auto &b) {
            line_seg upper_pair{a.rmost, b.rmost};
            line_seg lower_pair{a.lmost, b.lmost};
            return upper_pair.diff() + lower_pair.diff();
        };
        const auto te_pwd = PairMatrix::from_sample (
            sample_supporting_pileup.template_endpoints,
            mannd
        );

        // nearest neighbour monte carlo //
        std::optional<double> qpos_m1nn;
        stat_eval_s           qpos_m1nn_bgsim;  // compared to all reads
        stat_eval_s           qpos_m1nn_unisim; // compared to expected distribution
        if (qpos_pwd) {
            qpos_m1nn = medianNN (*qpos_pwd);
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
                const auto pwds = PairMatrix::from_sample (v, dist_1D);
                assert (pwds);
                const auto ret = medianNN (*pwds);
                return ret;
            };
            if (n_supporting_reads < 2) {
                qpos_m1nn_bgsim.err = "INSUFF_OBS";
            }
            else if (qpos_popv.size() < (n_supporting_reads * 2)) {
                // at a bare minimum, we want 2x more total samples than bg
                qpos_m1nn_bgsim.err = "INSUFF_BG";
            }  // TODO if less than e.g. 5/10 report low power?
            else {
                qpos_m1nn_bgsim = sim_to_bg (
                    *qpos_m1nn,
                    [&qpos_popv, &rng, n_supporting_reads] () {
                        return subsample_wo_replace(qpos_popv, n_supporting_reads, rng);
                    },
                    stat_fn,
                    log2_effsz
                );
            }

            // simulate against uniform distribution
            std::uniform_int_distribution<uint64_t> qpos_gen{0, read_len};
            qpos_m1nn_unisim = sim_to_bg(
                *qpos_m1nn,
                [&qpos_gen, &rng, n_supporting_reads] () {
                    std::vector<uint64_t> rand_qpos;
                    for (size_t i = 0; i < n_supporting_reads; ++i) {
                        rand_qpos.push_back(qpos_gen(rng));
                    }
                    return rand_qpos;
                },
                stat_fn,
                log2_effsz
            );
        } else {
            qpos_m1nn_bgsim.err = "INSUFF_OBS";
            qpos_m1nn_unisim.err = "INSUFF_OBS";
        }

        std::optional<double> te_m1nn;
        stat_eval_s           te_m1nn_sim;
        if (te_pwd) {
            te_m1nn = medianNN (*te_pwd);
            decltype(sample_supporting_pileup.template_endpoints) te_popv;
            if (!normal_only) {
                if (normal_pileup) {
                    te_popv.insert(
                        end(te_popv),
                        begin(sample_total_pileup.template_endpoints),
                        end(sample_supporting_pileup.template_endpoints)
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
                te_m1nn_sim.err = "INSUFF_OBS";
            }
            else if (te_popv.size() < n_supporting_templates * 2) {
                te_m1nn_sim.err = "INSUFF_BG";
            }
            else {
                te_m1nn_sim = sim_to_bg (
                    *te_m1nn,
                    [&rng, &te_popv, n_supporting_templates] () {
                        return subsample_wo_replace(te_popv, n_supporting_templates, rng);
                    },
                    [&mannd] (const auto &v) {
                        const auto pwds = PairMatrix::from_sample (
                            v,
                            mannd
                        );
                        assert (pwds);
                        const auto ret = medianNN (*pwds);
                        return ret;
                    },
                    [] (const auto &ev, const auto &simv) {
                        return log2 ((ev + 1) / (*mean (simv) + 1));
                    }
                );
            }
        } else {
            te_m1nn_sim.err = "INSUFF_OBS";
        }

        // --- MEDIAN LENGTH-NORMALISED ALIGNMENT SCORE --- //
        const auto  mlas = percentile (sample_supporting_pileup.normalised_as, 0.5);
        stat_eval_s mlas_sim;
        if (mlas) {
            decltype(sample_supporting_pileup.normalised_as) mlas_popv;
            if (!normal_only) {
                if (normal_pileup) {
                    mlas_popv.insert(
                        end(mlas_popv),
                        begin(sample_total_pileup.normalised_as),
                        end(sample_total_pileup.normalised_as)
                    );
                    mlas_popv.insert (
                        end(mlas_popv),
                        begin (normal_pileup->normalised_as),
                        end (normal_pileup->normalised_as)
                    );
                } else {
                    mlas_popv = sample_total_pileup.normalised_as;
                }
            } else {
                mlas_popv = normal_pileup->normalised_as;
            }

            if (n_supporting_reads < 2) {
                mlas_sim.err = "INSUFF_OBS";
            }
            else if (mlas_popv.size() < (2 * n_supporting_reads)) {
                mlas_sim.err = "INSUFF_BG";
            }
            else {
                mlas_sim = sim_to_bg (
                    *mlas,
                    [&rng, &mlas_popv, n_supporting_reads] () {
                        return subsample_wo_replace(mlas_popv, n_supporting_reads, rng);
                    },
                    [] (const std::vector<double> &v) {
                        const auto slas = percentile (v, 0.5);
                        assert (slas);
                        return *slas;
                    },
                    // effect size == raw delta
                    [] (const auto &ev, const auto &simv) {
                        return ev - *mean (simv);
                    }
                );
            }
        } else {
            mlas_sim.err = "INSUFF_OBS";
        }

        // consensus region of supporting templates
        // NOTE guarded earlier
        uint64_t lmosttc = std::numeric_limits<uint64_t>::max();
        uint64_t rmosttc = 0ULL;
        for (const auto &te : sample_supporting_pileup.template_endpoints) {
            if (te.lmost < lmosttc)
                lmosttc = te.lmost;
            if (te.rmost > rmosttc)
                rmosttc = te.rmost;
        }

        // TODO should really check if it's in bam header not just the vcf,
        // and that this is the correct reference
        // NOTE not all needed each loop
        auto rid_name = bcf_hdr_id2name (vcf_hdr.get(), b1->rid);
        std::optional<uint> kc;
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
            // TODO warn if all N
            kc.emplace (
                static_cast<uint> (round (entropy_lz76 (refs) * 100))
            );     // x100 scaling factor
        }

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
        if (mlas) {
            float val[3]{
                static_cast<float> (*mlas),
                static_cast<float> (mlas_sim.eff_sz.value_or (0.0)),
                static_cast<float> (mlas_sim.pval.value_or (0.0))
            };     // htslib requires conversion
            // TODO rounding
            write_info (FIELD_INF.at ("MLAS"), &val);
        }
        if (qpos_m1nn) {
            float val[3]{
                static_cast<float> (*qpos_m1nn),
                static_cast<float> (
                    qpos_m1nn_bgsim.eff_sz.value_or (0.0)
                ),
                static_cast<float> (qpos_m1nn_bgsim.pval.value_or (1.0))
            };
            write_info (FIELD_INF.at ("QM1NN"), &val);
        }
        if (te_m1nn) {
            float val[3]{
                static_cast<float> (*te_m1nn),
                static_cast<float> (te_m1nn_sim.eff_sz.value_or (0.0)),
                static_cast<float> (te_m1nn_sim.pval.value_or (1.0))
            };
            write_info (FIELD_INF.at ("TM1NN"), &val);
        }
        if (kc) {
            const auto val = *kc;
            write_info (FIELD_INF.at ("RCMPLX"), &val);
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
                "{}\t{}\t{}",
                rid_name,
                b1->pos + 1,
                opt_to_str<double> (mlas, "NA", rdbl2),
                opt_to_str<double> (mlas_sim.eff_sz, mlas_sim.err, rdbl2),
                opt_to_str<double> (mlas_sim.pval, mlas_sim.err, rdbl4),
                opt_to_str<double> (qpos_m1nn, "NA", rdbl2),
                opt_to_str<double> (qpos_m1nn_bgsim.eff_sz, qpos_m1nn_bgsim.err, rdbl2),
                opt_to_str<double> (qpos_m1nn_bgsim.pval, qpos_m1nn_bgsim.err, rdbl4),
                opt_to_str<double> (qpos_m1nn_unisim.eff_sz, qpos_m1nn_unisim.err, rdbl2),
                opt_to_str<double> (qpos_m1nn_unisim.pval, qpos_m1nn_unisim.err, rdbl4),
                opt_to_str<double> (te_m1nn, "NA", rdbl2),
                opt_to_str<double> (te_m1nn_sim.eff_sz, te_m1nn_sim.err, rdbl2),
                opt_to_str<double> (te_m1nn_sim.pval, te_m1nn_sim.err, rdbl4),
                opt_to_str (kc, "NA"),
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
