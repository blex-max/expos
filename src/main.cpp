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
#include <fstream>
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
    std::string type;
    std::string desc;
};
const std::unordered_map<std::string, field_s> FIELD_INF{
    {"MLAS",
     {"MLAS",
      "Float",
      "Median Length-normalised Alignment Score (AS) of "
      "reads supporting variant"}},
    {"QM1NN",
     {"QM1NN",
      "Float",
      "Median nearest neighbour distance of qpos on reads "
      "supporting variant"}},
    {"QES",
     {"QES",
      "Float",
      "Log2 ratio effect size of QM1NN against non-supporting reads, "
      "from monte carlo simulation"}},
    {"QPV",
     {"QPV",
      "Float",
      "P-value of QM1NN against non-supporting reads, from monte "
      "carlo simulation"}},
    {"TM1NN",
     {"TM1NN",
      "Float",
      "Median nearest neighbour distance of template endpoints as "
      "calcluated from read pairs supporting variant"}},
    {"TES",
     {"TES",
      "Float",
      "Log2 ratio effect size of TM1NN against non-supporting "
      "templates, from monte carlo simulation"}},
    {"TPV",
     {"TPV",
      "Float",
      "P-value of TM1NN against non-supporting templates, from monte "
      "carlo simulation"}},
    {"KC",
     {"KC",
      "Integer",
      "Kolmogorov Complexity of region spanned by supporting "
      "templates, scaled by x100"}}     // NOTE -- kc is optional
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

constexpr auto make_hdr_line (
    std::string name,
    std::string n,
    std::string t,
    std::string desc
) {
    return std::format (
        "##INFO=<ID={},Number={},Type={},Description=\"{}\",Source="
        "\"{}\",Version=\"{}\">",
        name,
        n,
        t,
        desc,
        PROG_NAME,
        VERSION
    );
}


// TODO allow user defined samflags for include/exclude
// TODO support single ended? LATER
// TODO options for calculating subset of data
// TODO options for more vcf data (e.g. REF,ALT) in TSV (if using expos as "genome browser by numbers")
// TODO optionally use positional data from NORMAL/BULK/SOMATIC
// to as background for simulation (if it's the same protocol which I don't know - if possible great!)
// TODO compare to uniform distribution (less valuable than background but possibly useful if e.g. not enough reads otherwise)
int main (
    int   argc,
    char *argv[]
) {
    namespace fs = std::filesystem;
    // TODO verify paths
    fs::path                 vcf_path;
    fs::path                 aln_path;
    fs::path                 ref_path;
    fs::path                 otsv_path;
    std::vector<std::string> flt_inc;
    std::vector<std::string> flt_ex;
    // std::vector<std::string> wfields;

    cxxopts::Options options (
        "expos",
        "EXtract POSitional data and statistics from alignment at "
        "VCF-specified pileups. Annotated VCF to stdout.\n"
    );

    // clang-format off
    options.add_options() ("h,help", "Print usage")
        // POSITIONAL
        ("vcf", "VCF", cxxopts::value<fs::path>())
        ("aln", "Sample BAM", cxxopts::value<fs::path>())

        // OPTS
        ("i,include",
         "Only operate on VCF records with this value present in FILTER. May be passed multiple times.",
         cxxopts::value<std::vector<std::string>>()) // multiple allowed
        ("e,exclude",
         "Only operate on VCF records without this value present in FILTER. May be passed multiple times.",
         cxxopts::value<std::vector<std::string>>()) // multiple allowed
        // ("w,write",
        //  "Write specified field to output VCF. May be passed multiple times.",
        //  cxxopts::value<std::vector<std::string>>()->default_value("ALL"))
        ("t,tsv",
         "Write a tsv of the output to file specified.",
         cxxopts::value<fs::path>())
        ("r,ref",
         "Alignment Reference Fasta for optionally adding template kolmogorov complexity to statistics.",
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

        if (parsedargs.count ("include")) {
            flt_inc = parsedargs["include"].as<std::vector<std::string>>();
        }
        if (parsedargs.count ("exclude")) {
            flt_ex = parsedargs["exclude"].as<std::vector<std::string>>();
        }

        if (parsedargs.count ("tsv")) {
            otsv_path = parsedargs["tsv"].as<fs::path>();
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


    // outputs
    htsFile_upt ovcf{hts_open ("-", "w"), hts_close};     // stdout
    bcf_hdr_upt ohdr{bcf_hdr_dup (vph.get()), bcf_hdr_destroy};

    // ADD LINES TO HDR
    // NOTE -- encodes a single field, regardless of number of samples (for now)
    for (const auto &l : FIELD_INF) {
        const auto &i = l.second;
        if (bcf_hdr_append (
                ohdr.get(),
                make_hdr_line (i.name, "1", i.type, i.desc).c_str()
            )
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
        *otsv << "CHROM\tPOS\t\tMLAS\tQPOS_M1NN\tQPOS_M1NN_"
                 "EFFSZ\tQPOS_M1NN_PVAL\t"
                 "TEMPL_M1NN\tTEMPL_M1NN_EFFSZ\tTEMPL_M1NN_PVAL\t"
                 "CONSENSUS_CMPLXx100\tNALT_READS\t"
                 "NTOTAL_READS"
              << "\n";
    }


    bool firsti = true;
    while (bcf_read (vp.get(), vph.get(), b1.get()) == 0) {
        if (firsti) {
            for (const auto &f : flt_inc) {
                if (bcf_has_filter (
                        vph.get(),
                        b1.get(),
                        const_cast<char *> (f.c_str())
                    )
                    < 0)
                    throw std::runtime_error (
                        "Unknown --include filter"
                    );     // unrecoverable
            }
            std::vector<std::string> tmp_ex;
            for (size_t i = 0; i < flt_ex.size(); ++i) {
                const auto &f = flt_ex[i];
                if (bcf_has_filter (
                        vph.get(),
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

        const auto iflt = has_filters (vph.get(), b1.get(), flt_inc);
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
        const auto eflt = has_filters (vph.get(), b1.get(), flt_ex);
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
            ap.get(),
            apit.get(),
            b1.get(),
            mtype
        );


        auto &altd   = vard.alt;
        using Tqposv = decltype (altd.qp);
        using Tqpos  = Tqposv::value_type;
        using Ttev   = decltype (altd.te);
        using Tte    = Ttev::value_type;

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

        if (altd.qp.empty() || altd.te.empty()) {
            std::cerr << std::format (
                "no supporting reads found for variant {}, "
                "skipping.",
                b1->d.id
            ) << std::endl;
            continue;
        }

        // --- MEDIAN LENGTH-NORMALISED ALIGNMENT SCORE --- //
        const auto mlas = percentile_from_sorted (altd.las, 0.5);


        // --- GET PAIRWISE DISTANCES --- //
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


        // --- NEAREST NEIGHBOUR MONTE CARLO --- //
        std::optional<double> qpos_ann;
        stat_eval_s           qpos_annsim;
        if (qpos_pwd) {
            qpos_ann    = medianNN (*qpos_pwd);
            qpos_annsim = sim_to_bg<double, Tqpos> (
                *qpos_ann,
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
                [&qpos_ann] (const auto s) { return s <= qpos_ann; }
            );
        } else {
            qpos_annsim.err = "INSUFF_OBS";
        }

        std::optional<double> te_ann = medianNN (*te_pwd);
        stat_eval_s           te_annsim;
        if (te_pwd) {
            te_ann    = medianNN (*te_pwd);
            te_annsim = sim_to_bg<double, Tte> (
                *te_ann,
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
                [&te_ann] (const auto s) { return s <= te_ann; }
            );
        } else {
            te_annsim.err = "INSUFF_OBS";
        }

        // TODO guard
        // consensus region of supporting templates
        uint64_t lmosttc = std::numeric_limits<uint64_t>::max();
        uint64_t rmosttc = 0ULL;
        for (const auto &te : altd.te) {
            if (te.lmost < lmosttc)
                lmosttc = te.lmost;
            if (te.rmost > rmosttc)
                rmosttc = te.rmost;
        }

        // TODO should really check if it's in bam header, and that this is the correct reference
        auto rid_name = bcf_hdr_id2name (vph.get(), b1->rid);
        std::optional<uint> kolmc;
        if (rp) {
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
                rp.get(),
                rid_name,
                lmosttc,
                rmosttc
            );
            // TODO warn if all N
            kolmc.emplace (
                static_cast<uint> (round (nk_lz76 (refs) * 100))
            );     // x100 scaling factor
        }

        // encode to vcf
        // TODO check0
        auto write_info =
            [&] (std::string nm, const void *val, int type) {
                if (bcf_update_info (
                        ohdr.get(),
                        b1.get(),
                        nm.c_str(),
                        val,
                        1,
                        type
                    )
                    != 0) {
                    throw std::runtime_error (
                        std::format (
                            "failed to write data {} as INFO field "
                            "in output VCF",
                            nm
                        )
                    );
                }
            };
        // TODO should probably encode missingness into the vcf somehow... like an EXPOS_ERR info field
        // TODO write other fields!
        // TODO extend the shorthand func further
        if (mlas) {
            const auto val = static_cast<float> (
                *mlas
            );     // htslib requires conversion
            // TODO rounding
            write_info ("MLAS", &val, BCF_HT_REAL);
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
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                rid_name,
                b1->pos + 1,
                opt_to_str<double> (mlas, "NA", rdbl2),
                opt_to_str<double> (qpos_ann, "NA", rdbl2),
                opt_to_str<double> (qpos_annsim.eff_sz, qpos_annsim.err, rdbl2),
                opt_to_str<double> (qpos_annsim.pval, qpos_annsim.err, rdbl4),
                opt_to_str<double> (te_ann, "NA", rdbl2),
                opt_to_str<double> (te_annsim.eff_sz, te_annsim.err, rdbl2),
                opt_to_str<double> (te_annsim.pval, te_annsim.err, rdbl4),
                opt_to_str (kolmc, "NA"),
                std::to_string (altd.qp.size()),     // n supporting reads
                std::to_string (qpos_popv.size())
            ) << "\n";
        }
        // clang-format on
    };

    std::cerr << "complete" << std::endl;
    return 0;
}
