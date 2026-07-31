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

#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <plog/Appenders/ColorConsoleAppender.h>
#include <plog/Formatters/TxtFormatter.h>
#include <plog/Init.h>
#include <plog/Log.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cxxopts.hpp>
#include <format>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>
#include <stdexcept>
#include <unordered_map>

#include "cli.hpp"
#include "expos/support.hpp"
#include "hts_ptr_t.hpp"
#include "lib-stats/monte-carlo.hpp"
#include "lib-stats/string.hpp"
#include "pileup.hpp"
#include "variant_record.hpp"

constexpr std::string PROG_NAME = "expos";
constexpr std::string VERSION = "0.0.0";
struct field_s {
  std::string name;
  std::string desc;
  int type, nrec;
};
// TODO sub in search radius
const std::unordered_map<std::string, field_s> FIELD_INF{
    {"QRK",
     {"QRK",
      "Array detailing Monte-Carlo simulation results for "
      "spatial clustering in mutant query position"
      "[1]standardised effect size (z-score) from comparisons "
      "to simulation against all reads;"
      "[2]one-sided P-value from comparisons to simulation "
      "against all reads."
      "Effect sizes greater than ~2.0 with significant p-value "
      "may indicate a spurious variant.",
      BCF_HT_REAL, 2}},
    {"TRK",
     {"TRK",
      "Array detailing Monte-Carlo simulation results for "
      "spatial clustering in endpoints of mutant templates:"
      "[1]standardised effect size (z-score) from comparisons "
      "to simulation against all reads;"
      "[2]one-sided P-value from comparisons to simulation "
      "against all reads."
      "Effect sizes greater than ~2.0 with significant p-value "
      "may indicate a spurious variant.",
      BCF_HT_REAL, 2}},
    {"RCMPLX",
     {"RCMPLX",
      "Mean 100-base window complexity (Lempel-Ziv estimated "
      "entropy rate) of "
      "the reference region spanned by supporting templates, "
      "scaled by x100.",
      BCF_HT_REAL, 1}},
    {"MLAS",
     {"MLAS",
      "Array of median read-length normalised alignment scores:"
      "[0]of reads supporting variant;"
      "[1]of all queried reads covering the variant location in "
      "the sample alignment.",
      BCF_HT_REAL, 2}}
};

// TODO optionally encode deviation from uniform to VCF
// TODO consider multiple and adjustable t for clustering assessment -- if so encode t in field name like QRK-5
// TODO consider MFE/n as report for secondary structure propensity
// TODO fraction of supporting reads with soft clipping, eff sz, pval, and median number of clipped bases
// in reads with soft clipping (CLPM)
// TODO options for more vcf data (e.g. REF,ALT) in TSV (if using expos as "genome browser by numbers")
// TODO consider assessment of cigar complexity/edit distance to ref of supporting reads
// compared to total population accounting for variant - (better than AS perhaps since AS is confounded with variant presence)
int main (int argc, char* argv[])
{
  std::string cmd_str;
  for (int i = 0; i < argc; ++i) {
    if (i != 0) {
      cmd_str += ' ';
    }
    cmd_str += argv[i];
  }
  auto now = std::chrono::system_clock::now();
  std::string datetime_str =
      std::format ("{:%Y-%m-%dT%H:%M:%S}", now);

  auto cli = cli::setup_cli();
  cli::ExposArgs ctx;

  try {
    auto parsedargs = cli.parse (argc, argv);
    if (parsedargs.contains ("help") != 0u) {
      std::cout << cli.help() << std::endl;
      return 0;
    }
    cli::parse_cli (ctx, parsedargs);
  }
  catch (const std::exception& e) {
    std::cerr << "Error parsing CLI options: " << e.what()
              << "\nTry expos --help.";
    return EXIT_FAILURE;
  }

  static plog::ColorConsoleAppender<plog::TxtFormatter>
      consoleAppender (plog::streamStdErr);
  plog::init (
      ctx.debug_mode ? plog::debug : plog::none, &consoleAppender
  );

  // TODO these uniq ptrs were not worthwhile
  // inputs
  auto *_ain{hts_open (ctx.aln_path.c_str(), "r")};
  if (_ain == NULL) {
    std::cerr << std::format (
                     "Could not open alignment file at {}",
                     ctx.aln_path.string()
                 )
              << std::endl;
    return EXIT_FAILURE;
  }
  auto *_aixin{sam_index_load (_ain, ctx.aln_path.c_str())};
  if (_aixin == NULL) {
    std::cerr << std::format (
                     "Coud not open index for alignment "
                     "file. Searched for {}.bai",
                     ctx.aln_path.c_str()
                 )
              << std::endl;
    return EXIT_FAILURE;
  }
  const htsFile_upt alnfh{_ain, hts_close};
  const hts_idx_upt aln_idx{_aixin, hts_idx_destroy};

  auto *_vin{hts_open (ctx.vcf_path.c_str(), "r")};
  if (_vin == NULL) {
    std::cerr << std::format (
                     "Could not open VCF file at {}",
                     ctx.vcf_path.string()
                 )
              << std::endl;
    return EXIT_FAILURE;
  }
  auto *_vh{bcf_hdr_read (_vin)};
  if (_vh == NULL) {
    std::cerr << std::format (
                     "Could not read header of VCF file at {}",
                     ctx.vcf_path.string()
                 )
              << std::endl;
    return EXIT_FAILURE;
  }
  const htsFile_upt vcffh{_vin, hts_close};
  const bcf_hdr_upt vcf_hdr{_vh, bcf_hdr_destroy};

  fai_upt reffh{nullptr, fai_destroy};
  if (!ctx.ref_path.empty()) {
    auto *_fin = fai_load (ctx.ref_path.string().c_str());
    if (_fin == NULL) {
      std::cerr << std::format (
                       "Could not read reference fasta at {}",
                       ctx.ref_path.string()
                   )
                << std::endl;
      return EXIT_FAILURE;
    }
          reffh.reset (_fin);
   
  }

  std::vector<std::pair<htsFile_upt, hts_idx_upt>> norms;
  for (const auto& norm_path : ctx.norm_paths) {
    auto *_nin{hts_open (norm_path.c_str(), "r")};
    if (_nin == NULL) {
      std::cerr << std::format (
                       "Could not open alignment file at {}",
                       norm_path.string()
                   )
                << std::endl;
      return EXIT_FAILURE;
    }
    auto *_nixin{sam_index_load (_nin, norm_path.c_str())};
    if (_nixin == NULL) {
      std::cerr << std::format (
                       "Coud not open index for alignment "
                       "file. Searched for {}.bai",
                       norm_path.c_str()
                   )
                << std::endl;
      return EXIT_FAILURE;
    }
    norms.emplace_back (
        htsFile_upt{_nin, hts_close},
        hts_idx_upt{_nixin, hts_idx_destroy}
    );
  }

  // outputs
  htsFile_upt ovcf{
      hts_open ("-", (ctx.no_gz ? "w" : "wz")), hts_close
  };  // stdout
  bcf_hdr_upt ohdr{bcf_hdr_dup (vcf_hdr.get()), bcf_hdr_destroy};

  // for adding to vcf header
  constexpr auto make_info_line = [] (field_s i) {
    std::string t2s[4];
    t2s[BCF_HT_INT] = "Integer";
    t2s[BCF_HT_REAL] = "Float";
    return std::format (
        "##INFO=<ID={},Number={},Type={},Description=\"{}\","
        "Source="
        "\"{}\",Version=\"{}\">",
        i.name, i.nrec, t2s[i.type], i.desc, PROG_NAME, VERSION
    );
  };
  for (const auto& l : FIELD_INF) {
    const auto& i = l.second;
    if (i.name == "RCMPLX" && !reffh) {
      continue;
    }
    if (bcf_hdr_append (
            ohdr.get(), make_info_line (i).c_str()
        ) != 0) {
      throw std::runtime_error (
          std::format (
              "failed to append to hdr for field {}", i.name
          )
      );
    }
  }

  if (bcf_hdr_append (
          ohdr.get(),
          std::format (
              "##source=\"expos v{} {}\"", VERSION, datetime_str
          )
              .c_str()
      ) != 0) {
    throw std::runtime_error ("failed to append source to hdr");
  }
  if (bcf_hdr_append (
          ohdr.get(),
          std::format ("##expos_cmd=\"{}\"", cmd_str).c_str()
      ) != 0) {
    throw std::runtime_error (
        "failed to append expos_cmd to hdr"
    );
  }

  if (bcf_hdr_sync (ohdr.get()) < 0) {
    throw std::runtime_error (
        "failed to sync VCF output header"
    );
  }
  if (bcf_hdr_write (ovcf.get(), ohdr.get()) < 0) {
    throw std::runtime_error (
        "failed to write VCF output header"
    );
  };

  // optional tsv output
  std::optional<std::ofstream> otsv;
  if (!ctx.otsv_path.empty()) {
    otsv.emplace (ctx.otsv_path);
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
             "REF_CMPLXx100\t"
             "REF_TOP_STR_RUN_MEAN\t"
             "REF_MAX_STR_RUN_LEN\t"
             "LMOST_TEMPLATE_START\t"
             "RMOST_TEMPLATE_END\t"
             "NALT_READS\t"
             "NTOTAL_READS"
          << "\n";
  }

  std::mt19937 rng{};
  rng.seed (ctx.seed);
  bool checked_read_len = false;
  constexpr double READ_LEN_TOL = 0.20;
  bcf1_upt b1{bcf_init(), bcf_destroy};

  for (const auto& f : ctx.flt_inc) {
    if (bcf_has_filter (
            vcf_hdr.get(), b1.get(),
            const_cast<char*> (f.c_str())
        ) == -1) {
      std::cerr << std::format (
                       "Unknown --include-records filter {}, "
                       "not present in "
                       "VCF",
                       f
                   )
                << std::endl;
      return EXIT_FAILURE;
    }
  }

  std::vector<std::string> tmp_ex;
  for (size_t i = 0; i < ctx.flt_exc.size(); ++i) {
    const auto& f = ctx.flt_exc[i];
    if (bcf_has_filter (
            vcf_hdr.get(), b1.get(),
            const_cast<char*> (f.c_str())
        ) == -1) {
      std::cerr << std::format (
                       "Warning: Unknown --exclude-records "
                       "filter {} not "
                       "present in VCF, "
                       "ignoring",
                       f
                   )
                << std::endl;
    }
    else {
      tmp_ex.push_back (ctx.flt_exc[i]);
    }
  }
  ctx.flt_exc = tmp_ex;

  auto write_record = [&]() {
    return bcf_write (ovcf.get(), ohdr.get(), b1.get());
  };
  while (bcf_read (vcffh.get(), vcf_hdr.get(), b1.get()) == 0) {
    if (!var::check_record_filters (
            b1.get(), vcf_hdr.get(), {ctx.flt_inc, ctx.flt_exc}
        )) {
      if (write_record() != 0) {
        goto ERR_WRITE;
      }
      else {
        continue;
      }
    }

    if (b1->n_allele != 2) {
      std::cerr
          << std::format (
                 "Variant {}:{} ({}) has more than two "
                 "(REF,ALT) "
                 "alleles. "
                 "Unnormalised variant calls are not supported. "
                 "Skipping.",
                 bcf_hdr_id2name (vcf_hdr.get(), b1->rid),
                 b1->pos + 1, b1->d.id
             )
          << std::endl;
      if (write_record() != 0) {
        goto ERR_WRITE;
      }
      else {
        continue;
      }
    }

    auto mtype = bcf_has_variant_type (
        b1.get(), 1, VCF_DEL | VCF_INS | VCF_SNP | VCF_MNP
    );
    std::function<bool (const bcf1_t*, const bam_pileup1_t*)>
        support_eval_fn;
    switch (mtype) {
      // one and only one of
      case (VCF_SNP):
        support_eval_fn = &eval_support_snp;
        break;
      case (VCF_MNP):
        support_eval_fn = &eval_support_mnp;
        break;
      case (VCF_DEL):
        support_eval_fn = &eval_support_del;
        break;
      case (VCF_INS):
        support_eval_fn = &eval_support_ins;
        break;
      default:
        std::cerr << std::format (
                         "Variant {} has an unsupported/complex "
                         "mutation, or could not be typed"
                         "skipping.",
                         b1->d.id
                     )
                  << std::endl;
        if (write_record() != 0) {
          goto ERR_WRITE;
        }
        else {
          continue;
        }
    }

    auto [sample_supporting_pileup, sample_total_pileup] =
        pileup_group_and_anaylse (
            alnfh.get(), aln_idx.get(), b1.get(), ctx.flag_inc,
            ctx.flag_exc, support_eval_fn, vcf_hdr.get()
        );
    std::optional<PileupMetrics> normal_pileup;
    for (auto& n : norms) {
      auto pm = pileup_analyse (
          n.first.get(), n.second.get(), b1.get(), ctx.flag_inc,
          ctx.flag_exc, vcf_hdr.get()
      );
      normal_pileup = normal_pileup
                          ? merge_pileup_metrics (
                                std::move (*normal_pileup), pm
                            )
                          : pm;
    }

    // guard incorrect read len
    if (!checked_read_len &&
        sample_total_pileup.read_lengths.size() >= 5) {
      auto lens = sample_total_pileup.read_lengths;
      std::nth_element (
          begin (lens), begin (lens) + lens.size() / 2,
          end (lens)
      );
      const double median_len = lens[lens.size() / 2];
      const double deviation =
          std::abs (
              median_len - static_cast<double> (ctx.exp_read_len)
          ) /
          static_cast<double> (ctx.exp_read_len);
      if (deviation > READ_LEN_TOL) {
        std::cerr
            << std::format (
                   "error: observed median read length {} "
                   "deviates from "
                   "expected {} by {:.0f}% (tolerance {}%). "
                   "Perhaps re-run with the correct -l value?",
                   static_cast<uint32_t> (median_len),
                   ctx.exp_read_len, deviation * 100.0,
                   READ_LEN_TOL * 100.0
               )
            << std::endl;
        return EXIT_FAILURE;
      }
      checked_read_len = true;
    }

    auto n_obs_qpos =
        sample_supporting_pileup.query_position.size();
    auto n_obs_templ =
        sample_supporting_pileup.template_endpoints.size();

    if (n_obs_qpos == 0 || n_obs_templ == 0) {
      std::cerr << std::format (
                       "no supporting reads found for variant "
                       "{}:{} ({}), "
                       "skipping.",
                       bcf_hdr_id2name (vcf_hdr.get(), b1->rid),
                       b1->pos + 1, b1->d.id
                   )
                << std::endl;
      // NOTE skips write to tsv
      if (write_record() != 0) {
        goto ERR_WRITE;
      }
      else {
        continue;
      }
    }

    // TODO this should be a skip I think,
    // with a cli flag override
    if (normal_pileup &&
        normal_pileup->query_position.empty()) {
      std::cerr
          << std::format (
                 "Warning: no reads covering variant location "
                 "found in normal/s for variant {}:{} ({})",
                 bcf_hdr_id2name (vcf_hdr.get(), b1->rid),
                 b1->pos + 1, b1->d.id
             )
          << std::endl;
    }

    // --- REF COMPLEXITY --- //
    // consensus region of supporting templates
    uint64_t lmosttc = std::numeric_limits<uint64_t>::max();
    uint64_t rmosttc = 0ULL;
    for (const auto& te :
         sample_supporting_pileup.template_endpoints) {
      lmosttc = std::min(te.lmost, lmosttc);
      rmosttc = std::max(te.rmost, rmosttc);
    }
    const auto span_length = rmosttc - lmosttc;

    // TODO should really check if it's in bam header not just the vcf,
    // and that this is the correct reference
    const auto *rid_name = bcf_hdr_id2name (vcf_hdr.get(), b1->rid);
    std::optional<size_t> ref_entropy;
    std::optional<double> run_mean;
    std::optional<size_t> run_max;
    if (reffh) {
      if (rid_name == NULL) {
        std::cerr << std::format (
                         "Could not find reference ID {} in VCF "
                         "header "
                         "- VCF "
                         "misformatted?",
                         b1->rid
                     )
                  << std::endl;
        std::cerr << "skipping reference complexity analysis";
      }
      else {
        std::string refs =
            fai_fetch (reffh.get(), rid_name, lmosttc, rmosttc);
        // transform to ignore masks
        transform (
            begin (refs), end (refs), begin (refs), ::toupper
        );
        if (refs.find ('N') == std::string::npos) {  // No Ns
          // lz complexity
          const size_t window_size = 100;
          double cmplx_sum = 0;
          size_t n_win = 0;
          for (; (n_win + window_size) <= refs.size();
               ++n_win) {  // ++n_win == step of 1
            cmplx_sum += string_stats::entropy_lz76 (
                refs.substr (n_win, window_size)
            );
          }
          const auto mean_window_entropy =
              cmplx_sum /
              static_cast<double> (n_win);
          ref_entropy.emplace (
              round (mean_window_entropy * 100)
          );  // x100 scaling factor, which may not be a great idea
        }

        // quantify homopolymer & STR presence
        if (ctx.assess_microhom) {
          auto str_runs = string_stats::periodic_rle (refs, 5);
          auto str_runs_n = str_runs.size();
          auto top10_i = static_cast<size_t> (
              floor (static_cast<double> (str_runs.size()) * 0.9)
          );
          auto top10_n = str_runs_n - top10_i;
          std::nth_element (
              begin (str_runs),
              begin (str_runs) + static_cast<long> (top10_i),
              end (str_runs)
          );
          decltype (str_runs)::value_type str_run_max = 0;
          decltype (str_runs)::value_type top10_run_sum = 0;
          for (size_t i = top10_i; i < str_runs_n; ++i) {
            const auto e = str_runs[i];
            str_run_max = std::max(e, str_run_max);
            top10_run_sum += e;
          }
          run_max = str_run_max;
          run_mean = static_cast<double> (top10_run_sum) /
                     static_cast<double> (top10_n);
        }
      }
    }

    // --- CLUSTERING ANALYSIS - RIPLEY'S K --- //
    // results for qpos
    std::optional<double> qpos_K;
    monte_carlo::stat_eval_s
        qpos_K_bgsim;  // compared to all reads
    monte_carlo::stat_eval_s
        qpos_K_unisim;  // compared to expected distribution

    // get pairwise distances
    // simple 1D dist for query position
    constexpr auto dist_1D = [] (const auto& a, const auto& b) {
      return (a > b) ? (a - b) : (b - a);
    };
    const auto qpos_pwd = spatial::PairMatrix::from_sample (
        sample_supporting_pileup.query_position,
        dist_1D
    );  // empty if <2 samples

    if (n_obs_qpos < 2 || !qpos_pwd) {
      qpos_K_bgsim.err = "INSUFF_OBS";
      if (ctx.uniform_sim) {
        qpos_K_unisim.err = "INSUFF_OBS";
      }
    }
    else {
      const size_t search_radius =
          5;  // bases. Sensible values << read length.

      // mutant sample value
      qpos_K = spatial::ripley_k (
          *qpos_pwd, search_radius,
          static_cast<double> (qpos_pwd->dim()) /
              static_cast<double> (ctx.exp_read_len)
      );

      // assemble total sample pool
      decltype (sample_supporting_pileup
                    .query_position) qpos_popv;
      if (!ctx.normal_only) {
        if (normal_pileup) {  // ADD NORMAL OBS
          assert (normal_pileup);
          qpos_popv.insert (
              end (qpos_popv),
              begin (sample_total_pileup.query_position),
              end (sample_total_pileup.query_position)
          );
          qpos_popv.insert (
              end (qpos_popv),
              begin (normal_pileup->query_position),
              end (normal_pileup->query_position)
          );
        }
        else {  // total reads only
          qpos_popv = sample_total_pileup.query_position;
        }
      }
      else {  // reads from normal only
        assert (normal_pileup);
        qpos_popv = normal_pileup->query_position;
      }

      // simulation function
      auto stat_fn = [&] (const auto& v) {
        const auto pwds =
            spatial::PairMatrix::from_sample (v, dist_1D);
        assert (pwds);
        return spatial::ripley_k (
            *pwds, search_radius,
            static_cast<double> (v.size()) /
                static_cast<double> (ctx.exp_read_len)
        );
      };
      // guard, then run simulation
      if (qpos_popv.size() < (n_obs_qpos * 2)) {
        // at a bare minimum, we want 2x more total samples than bg
        qpos_K_bgsim.err = "INSUFF_BG";
        std::cerr << std::format (
                         "Warning: insufficient total reads "
                         "(<2x supporting) for variant {}:{} "
                         "({}) for anaylsis",
                         bcf_hdr_id2name (
                             vcf_hdr.get(), b1->rid
                         ),
                         b1->pos + 1, b1->d.id
                     )
                  << std::endl;
      }
      else {
        qpos_K_bgsim = monte_carlo::sim_to_bg (
            *qpos_K,
            [&]() {
              return monte_carlo::subsample_wo_replace (
                  qpos_popv, n_obs_qpos, rng
              );
            },
            stat_fn, monte_carlo::ses
        );
      }

      // simulate against uniform distribution
      if (ctx.uniform_sim) {
        // NOTE:
        // the same deviation from the null will always have the same
        // p value for a fixed sample size and read length, therefore
        // could precompute this and save wasted comuptation:
        // scale/normalise the statistic
        // S = (sample_size / read_length) * (Ripley)
        // Precompute one high-quality null CDF for S with a large sample size (>40)
        // Use that CDF for all datasets with moderate/large n.
        // for smaller n, use exactly stored CDF
        std::uniform_int_distribution<uint64_t> qpos_gen{
            0, ctx.exp_read_len - 1
        };  // inclusive ends
        qpos_K_unisim = monte_carlo::sim_to_bg (
            *qpos_K,
            [&]() {
              std::vector<uint64_t> rand_qpos;
              rand_qpos.reserve (n_obs_qpos);
              for (size_t i = 0; i < n_obs_qpos; ++i) {
                rand_qpos.push_back (qpos_gen (rng));
              }
              return rand_qpos;
            },
            stat_fn, monte_carlo::ses
        );
      }
    }

    // --- TEMPLATE ENDPOINTS --- //
    // manhattan distance for template endpoints
    constexpr auto mannd = [] (const auto& a, const auto& b) {
      const spatial::line_seg upper_pair{a.rmost, b.rmost};
      const spatial::line_seg lower_pair{a.lmost, b.lmost};
      return upper_pair.diff() + lower_pair.diff();
    };
    const auto te_pwd = spatial::PairMatrix::from_sample (
        sample_supporting_pileup.template_endpoints, mannd
    );

    std::optional<double> te_K;
    monte_carlo::stat_eval_s te_K_bgsim;
    if (n_obs_templ < 2 || !te_pwd) {
      te_K_bgsim.err = "INSUFF_OBS";
    }
    else {
      const size_t search_radius = 6;
      te_K = spatial::ripley_k (
          *te_pwd, search_radius,
          static_cast<double> (te_pwd->dim()) /
              static_cast<double> (span_length)
      );

      // assemble total sample pool
      decltype (sample_supporting_pileup
                    .template_endpoints) te_popv;
      if (!ctx.normal_only) {
        if (normal_pileup) {
          te_popv.insert (
              end (te_popv),
              begin (sample_total_pileup.template_endpoints),
              end (sample_total_pileup.template_endpoints)
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
      }
      else {
        te_popv = normal_pileup->template_endpoints;
      }

      if (n_obs_templ < 2) {
        te_K_bgsim.err = "INSUFF_OBS";
        std::cerr << std::format (
                         "Warning: insufficient total templates "
                         "(<2x supporting) for variant {}:{} "
                         "({}) for analysis.",
                         bcf_hdr_id2name (
                             vcf_hdr.get(), b1->rid
                         ),
                         b1->pos + 1, b1->d.id
                     )
                  << std::endl;
      }
      else {
        te_K_bgsim = monte_carlo::sim_to_bg (
            *te_K,
            [&rng, &te_popv, n_obs_templ]() {
              return monte_carlo::subsample_wo_replace (
                  te_popv, n_obs_templ, rng
              );
            },
            [&mannd, span_length] (const auto& v) {
              const auto pwds =
                  spatial::PairMatrix::from_sample (v, mannd);
              assert (pwds);
              return spatial::ripley_k (
                  *pwds, search_radius,
                  static_cast<double> (pwds->dim()) /
                      static_cast<double> (span_length)
              );
            },
            monte_carlo::ses
        );
      }
    }

    // --- MEDIAN LENGTH-NORMALISED ALIGNMENT SCORE --- //
    const auto mlas_supporting = summary::percentile (
        sample_supporting_pileup.normalised_as, 0.5
    );
    const auto mlas_total = summary::percentile (
        sample_total_pileup.normalised_as, 0.5
    );

    // encode to vcf
    auto write_info = [&] (field_s i, const void* val) {
      if (bcf_update_info (
              ohdr.get(), b1.get(), i.name.c_str(), val, i.nrec,
              i.type
          ) != 0) {
        std::cerr << std::format (
                         "failed to write {} as INFO field "
                         "in output VCF",
                         i.name
                     )
                  << std::endl;
      }
    };
    if (qpos_K) {
      float val[2];
      if (qpos_K_bgsim.eff_sz) {
        val[0] = static_cast<float> (*qpos_K_bgsim.eff_sz);
      }
      else {
        bcf_float_set_missing (val[0]);
      }
      if (qpos_K_bgsim.pval) {
        val[1] = static_cast<float> (*qpos_K_bgsim.pval);
      }
      else {
        bcf_float_set_missing (val[1]);
      }
      write_info (FIELD_INF.at ("QRK"), &val);
    }
    if (te_K) {
      float val[2];
      if (te_K_bgsim.eff_sz) {
        val[0] = static_cast<float> (*te_K_bgsim.eff_sz);
      }
      else {
        bcf_float_set_missing (val[0]);
      }
      if (te_K_bgsim.pval) {
        val[1] = static_cast<float> (*te_K_bgsim.pval);
      }
      else {
        bcf_float_set_missing (val[1]);
      }
      write_info (FIELD_INF.at ("TRK"), &val);
    }
    if (ref_entropy) {
      float val = static_cast<float> (*ref_entropy);
      write_info (FIELD_INF.at ("RCMPLX"), &val);
    }
    if (mlas_supporting || mlas_total) {
      float val[2];
      if (mlas_supporting) {
        val[0] = static_cast<float> (*mlas_supporting);
      }
      else {
        bcf_float_set_missing (val[0]);
      }
      if (mlas_total) {
        val[1] = static_cast<float> (*mlas_total);
      }
      else {
        bcf_float_set_missing (val[1]);
      }
      write_info (FIELD_INF.at ("MLAS"), &val);
    }

    if (write_record() != 0) {
      goto ERR_WRITE;
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
                opt_to_str<double> (qpos_K, "NA", rdbl2),
                opt_to_str<double> (qpos_K_bgsim.eff_sz, qpos_K_bgsim.err, rdbl2),
                opt_to_str<double> (qpos_K_bgsim.pval, qpos_K_bgsim.err, rdbl4),
                opt_to_str<double> (qpos_K_unisim.eff_sz, qpos_K_unisim.err, rdbl2),
                opt_to_str<double> (qpos_K_unisim.pval, qpos_K_unisim.err, rdbl4),
                opt_to_str<double> (te_K, "NA", rdbl2),
                opt_to_str<double> (te_K_bgsim.eff_sz, te_K_bgsim.err, rdbl2),
                opt_to_str<double> (te_K_bgsim.pval, te_K_bgsim.err, rdbl4),
                opt_to_str (ref_entropy, "NA"),
                opt_to_str (run_mean, "NA"),
                opt_to_str (run_max, "NA"),
                std::to_string(lmosttc),
                std::to_string(rmosttc),
                std::to_string (sample_supporting_pileup.normalised_as.size()),
                std::to_string (sample_total_pileup.normalised_as.size())
            ) << "\n";
        }
    // clang-format on
  };

  std::cerr << "complete" << std::endl;
  return 0;

ERR_WRITE:
  std::cerr << std::format (
                   "failed to write record to output VCF"
               )
            << std::endl;
  return EXIT_FAILURE;
}
