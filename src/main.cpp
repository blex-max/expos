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
#include <format>
#include <fstream>
#include <functional>
#include <htslib/faidx.h>
#include <iostream>

#include <cxxopts.hpp>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
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
      "[1]delta (supporting - background) effect size and [2]P-value "
      "against "
      "background, from monte-carlo simulation",
      BCF_HT_REAL,
      3}},
    {"QM1NN",
     {"QM1NN",
      "[0]Median nearest neighbour distance of variant query "
      "position; [1]log2 ratio effect size and [2]P-value against "
      "background, from monte-carlo simulation",
      BCF_HT_REAL,
      3}},
    {"TM1NN",
     {"TM1NN",
      "[0]Median nearest neighbour distance of template endpoints "
      "from read pairs supporting variant; [1]log2 ratio effect size "
      "and [2]P-value against background, from "
      "monte-carlo simulation",
      BCF_HT_REAL,
      3}},
    {"KC",
     {"KC",
      "Kolmogorov Complexity of region spanned by supporting "
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
// TODO MCLP (CLPM), and simulate
// TODO add consensus span region back to tsv
// TODO options for calculating subset of data only
// TODO options for more vcf data (e.g. REF,ALT) in TSV (if using expos as "genome browser by numbers")
// NOTE could compare to uniform distribution (less valuable than background but possibly useful if e.g. not enough reads otherwise)
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
    std::vector<std::string> flt_ex;
    bool                     no_gz = false;
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
        // ("w,write",
        //  "Write specified field to output VCF. May be passed multiple times.",
        //  cxxopts::value<std::vector<std::string>>()->default_value("ALL"))
        ("t,tsv",
         "Write a tsv of extended statistics to file specified.",
         cxxopts::value<fs::path>())
        ("n,normal",
         "Alignment for use as additional background data for simulation",
         cxxopts::value<fs::path>())
        ("r,ref",
         "Alignment Reference Fasta for optionally adding template kolmogorov complexity to statistics.",
         cxxopts::value<fs::path>())
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
            flt_ex = parsedargs["exclude"].as<std::vector<std::string>>();
        }

        if (parsedargs.count ("tsv")) {
            otsv_path = parsedargs["tsv"].as<fs::path>();
        }

        if (parsedargs.count ("normal")) {
            norm_path = parsedargs["normal"].as<fs::path>();
            std::cerr << "Using normal: " << norm_path << std::endl;
        }

        if (parsedargs.count ("uncompressed")) {
            no_gz = true;
        }

    } catch (const std::exception &e) {
        std::cerr << "Error parsing CLI options: " << e.what()
                  << "\n";
        return 1;
    }

    // inputs
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
    htsFile_upt alnfh{std::move (_ain), hts_close};
    hts_idx_upt aln_idx{std::move (_aixin), hts_idx_destroy};

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
    htsFile_upt vcffh{std::move (_vin), hts_close};
    bcf_hdr_upt vcf_hdr{std::move (_vh), bcf_hdr_destroy};

    fai_upt reffh{nullptr, fai_destroy};
    if (!ref_path.empty()) {
        auto _fin = fai_load (ref_path.string().c_str());
        if (_fin == NULL) {
            std::cerr << std::format (
                "Could not read reference fasta at {}",
                ref_path.string()
            );
        } else {
            reffh.reset (_fin);
        }
    }

    // inputs
    std::optional<std::pair<htsFile_upt, hts_idx_upt>> norm;
    if (!norm_path.empty()) {
        auto _nin{hts_open (aln_path.c_str(), "r")};
        if (_nin == NULL) {
            std::cerr << std::format (
                "Could not open alignment file at {}",
                aln_path.string()
            ) << std::endl;
            return 1;
        }
        auto _nixin{sam_index_load (_nin, aln_path.c_str())};
        if (_nixin == NULL) {
            std::cerr << std::format (
                "Coud not open index for alignment "
                "file. Searched for {}.bai",
                aln_path.c_str()
            ) << std::endl;
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
        *otsv << "CHROM\tPOS\t\tMLAS\tMLAS_EFFSZ\tMLAS_PVAL\tQPOS_"
                 "M1NN\tQPOS_M1NN_"
                 "EFFSZ\tQPOS_M1NN_PVAL\t"
                 "TEMPL_M1NN\tTEMPL_M1NN_EFFSZ\tTEMPL_M1NN_PVAL\t"
                 "CONSENSUS_CMPLXx100\tNALT_READS\t"
                 "NTOTAL_READS"
              << "\n";
    }


    bool     firsti = true;
    bcf1_upt b1{bcf_init(), bcf_destroy};
    while (bcf_read (vcffh.get(), vcf_hdr.get(), b1.get()) == 0) {
        if (firsti) {
            for (const auto &f : flt_inc) {
                if (bcf_has_filter (
                        vcf_hdr.get(),
                        b1.get(),
                        const_cast<char *> (f.c_str())
                    )
                    < 0)
                    throw std::runtime_error (
                        std::format ("Unknown --include filter {}", f)
                    );     // unrecoverable
            }
            std::vector<std::string> tmp_ex;
            for (size_t i = 0; i < flt_ex.size(); ++i) {
                const auto &f = flt_ex[i];
                if (bcf_has_filter (
                        vcf_hdr.get(),
                        b1.get(),
                        const_cast<char *> (f.c_str())
                    )
                    < 0) {
                    std::cerr << std::format (
                        "Warning: Unknown --exclude filter {}, "
                        "ignoring",
                        f
                    ) << std::endl;
                } else {
                    tmp_ex.push_back (flt_ex[i]);
                }
            }
            flt_ex = tmp_ex;
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
                throw std::runtime_error (
                    std::format ("failed to write record to VCF")
                );
            };
            continue;
        }
        const auto eflt = has_filters (vcf_hdr.get(), b1.get(), flt_ex);
        if (std::any_of (begin (eflt), end (eflt), [] (const auto a) {
                return a;
            })) {
            if (bcf_write (ovcf.get(), ohdr.get(), b1.get()) != 0) {
                throw std::runtime_error (
                    std::format ("failed to write record to VCF")
                );
            };
            continue;
        }

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
            alnfh.get(),
            aln_idx.get(),
            b1.get(),
            mtype,
            true
        );
        std::optional<aln_obs> normd;
        if (norm) {
            normd = get_aln_data (
                norm->first.get(),
                norm->second.get(),
                b1.get(),
                mtype,
                true
            );
        }

        auto &altd = vard.alt;
        if (altd.qp.empty() || altd.te.empty()) {
            std::cerr << std::format (
                "no supporting reads found for variant {}, "
                "skipping.",
                b1->d.id
            ) << std::endl;
            continue;
        }

        // --- CLUSTERING ANALYSIS (nearest neighbour) --- //
        using Tqposv = decltype (altd.qp);
        using Tqpos  = Tqposv::value_type;
        using Ttev   = decltype (altd.te);
        using Tte    = Ttev::value_type;

        // get pairwise distances
        constexpr auto d1d = [] (const Tqpos &a,
                                 const Tqpos &b) -> Tqpos {
            return (a > b) ? (a - b) : (b - a);
        };
        const auto qpos_pwd = PairMatrix<Tqpos>::from_sample (
            altd.qp,
            d1d
        );     // empty if <2 samples
        // manhattan distance
        constexpr auto mannd = [] (const Tte &a,
                                   const Tte &b) -> uint64_t {
            line_seg upper_pair{a.rmost, b.rmost};
            line_seg lower_pair{a.lmost, b.lmost};
            return upper_pair.diff() + lower_pair.diff();
        };
        const auto te_pwd = PairMatrix<uint64_t>::from_sample (
            altd.te,
            mannd
        );

        // nearest neighbour monte carlo //
        size_t                nread;     // extra stat
        std::optional<double> qpos_m1nn;
        stat_eval_s           qpos_m1nn_sim;
        if (qpos_pwd) {
            qpos_m1nn = medianNN (*qpos_pwd);
            Tqposv qpos_popv;
            qpos_popv.insert (
                qpos_popv.end(),
                begin (vard.alt.qp),
                end (vard.alt.qp)
            );
            qpos_popv.insert (
                qpos_popv.end(),
                begin (vard.other.qp),
                end (vard.other.qp)
            );
            nread = qpos_popv.size();     // don't include normal
            if (normd) {                  // ADD NORMAL OBS
                qpos_popv.insert (
                    qpos_popv.end(),
                    begin (normd->other.qp),
                    end (normd->other.qp)
                );
            }
            qpos_m1nn_sim = sim_to_bg (
                *qpos_m1nn,
                altd.qp.size(),
                qpos_popv,
                [&d1d] (const Tqposv &v) {
                    const auto pwds = PairMatrix<Tqpos>::from_sample (
                        v,
                        d1d
                    );
                    assert (pwds);
                    const auto ret = medianNN (*pwds);
                    return ret;
                },
                [] (const auto ev, const auto sim) {
                    return sim <= ev;
                },
                // +1 removes confusing values when 0,
                // log makes effect size symmetric around 0
                // log2 means -1 = half the size of background
                // +1 = double the size of background
                [] (const auto ev, const auto simv) {
                    return log2 ((ev + 1) / (*mean (simv) + 1));
                }
            );
        } else {
            qpos_m1nn_sim.err = "INSUFF_OBS";
        }

        std::optional<double> te_m1nn;
        stat_eval_s           te_m1nn_sim;
        if (te_pwd) {
            te_m1nn = medianNN (*te_pwd);
            Ttev te_popv;
            te_popv.insert (
                te_popv.end(),
                begin (vard.alt.te),
                end (vard.alt.te)
            );
            te_popv.insert (
                te_popv.end(),
                begin (vard.other.te),
                end (vard.other.te)
            );
            if (normd) {     // ADD NORMAL OBS
                te_popv.insert (
                    te_popv.end(),
                    begin (normd->other.te),
                    end (normd->other.te)
                );
            }
            te_m1nn_sim = sim_to_bg (
                *te_m1nn,
                altd.te.size(),
                te_popv,
                [&mannd] (const Ttev &v) {
                    const auto pwds = PairMatrix<uint64_t>::from_sample (
                        v,
                        mannd
                    );
                    assert (pwds);
                    const auto ret = medianNN (*pwds);
                    return ret;
                },
                [] (const auto ev, const auto sim) {
                    return sim <= ev;
                },
                [] (const auto ev, const auto simv) {
                    return log2 ((ev + 1) / (*mean (simv) + 1));
                }
            );
        } else {
            te_m1nn_sim.err = "INSUFF_OBS";
        }

        // --- MEDIAN LENGTH-NORMALISED ALIGNMENT SCORE --- //
        const auto  mlas = percentile_from_sorted (altd.las, 0.5);
        stat_eval_s mlas_sim;
        if (mlas) {
            std::vector<double> mlas_popv;
            mlas_popv.insert (
                mlas_popv.end(),
                begin (vard.alt.las),
                end (vard.alt.las)
            );
            mlas_popv.insert (
                mlas_popv.end(),
                begin (vard.other.las),
                end (vard.other.las)
            );
            if (normd) {     // ADD NORMAL OBS
                mlas_popv.insert (
                    mlas_popv.end(),
                    begin (normd->other.las),
                    end (normd->other.las)
                );
            }
            mlas_sim = sim_to_bg (
                *mlas,
                altd.las.size(),
                mlas_popv,
                [&mannd] (const std::vector<double> &v) {
                    const auto slas = percentile_from_sorted (v, 0.5);
                    assert (percentile_from_sorted);
                    return *slas;
                },
                [] (const auto ev, const auto sim) {
                    return sim <= ev;
                },
                // effect size == raw delta
                [] (const auto ev, const auto simv) {
                    return ev - *mean (simv);
                }
            );
        } else {
            mlas_sim.err = "INSUFF_OBS";
        }

        // consensus region of supporting templates
        uint64_t lmosttc = std::numeric_limits<uint64_t>::max();
        uint64_t rmosttc = 0ULL;
        for (const auto &te : altd.te) {
            if (te.lmost < lmosttc)
                lmosttc = te.lmost;
            if (te.rmost > rmosttc)
                rmosttc = te.rmost;
        }

        // TODO should really check if it's in bam header not just the vcf,
        // and that this is the correct reference
        auto rid_name = bcf_hdr_id2name (vcf_hdr.get(), b1->rid);
        std::optional<uint> kc;
        if (reffh) {
            if (rid_name == NULL) {
                throw std::runtime_error (
                    std::format (
                        "Could not find rid {} in VCF header "
                        "- VCF "
                        "misformatted?",
                        b1->rid
                    )
                );
            }
            std::string refs = fai_autofetch (
                reffh.get(),
                rid_name,
                lmosttc,
                rmosttc
            );
            // TODO warn if all N
            kc.emplace (
                static_cast<uint> (round (nk_lz76 (refs) * 100))
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
                throw std::runtime_error (
                    std::format (
                        "failed to write data {} as INFO field "
                        "in output VCF",
                        i.name
                    )
                );
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
                    qpos_m1nn_sim.eff_sz.value_or (0.0)
                ),
                static_cast<float> (qpos_m1nn_sim.pval.value_or (1.0))
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
            write_info (FIELD_INF.at ("KC"), &val);
        }

        if (bcf_write (ovcf.get(), ohdr.get(), b1.get()) != 0) {
            throw std::runtime_error (
                std::format ("failed to write record to VCF")
            );
        };

        // report statistics to tsv
        // clang-format off
        if (otsv) {
            *otsv << std::format (
                "{}\t{}\t{}\t{}\t{}\t"
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                rid_name,
                b1->pos + 1,
                opt_to_str<double> (mlas, "NA", rdbl2),
                opt_to_str<double> (mlas_sim.eff_sz, mlas_sim.err, rdbl2),
                opt_to_str<double> (mlas_sim.pval, mlas_sim.err, rdbl4),
                opt_to_str<double> (qpos_m1nn, "NA", rdbl2),
                opt_to_str<double> (qpos_m1nn_sim.eff_sz, qpos_m1nn_sim.err, rdbl2),
                opt_to_str<double> (qpos_m1nn_sim.pval, qpos_m1nn_sim.err, rdbl4),
                opt_to_str<double> (te_m1nn, "NA", rdbl2),
                opt_to_str<double> (te_m1nn_sim.eff_sz, te_m1nn_sim.err, rdbl2),
                opt_to_str<double> (te_m1nn_sim.pval, te_m1nn_sim.err, rdbl4),
                opt_to_str (kc, "NA"),
                std::to_string (altd.qp.size()),     // n supporting reads
                std::to_string (nread)
            ) << "\n";
        }
        // clang-format on
    };

    std::cerr << "complete" << std::endl;
    return 0;
}
