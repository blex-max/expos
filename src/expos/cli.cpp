#include "cli.hpp"
#include "cxxopts.hpp"

#include <filesystem>
#include <iostream>
#include <format>

namespace cli {

cxo::Options setup_cli () {
    namespace fs = std::filesystem;

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
         "Alignment for use as additional background data for simulation. May be passed multiple times.",
         cxxopts::value<std::vector<fs::path>>())
        ("normal-only",
         "Use only reads from the provided normal as background data, excluding reads from the sample")

        ("i,include-records",
         "Only operate on VCF records with this value present in FILTER. e.g. -i PASS. May be passed multiple times.",
         cxxopts::value<std::vector<std::string>>()) // multiple allowed
        ("e,exclude-records",
         "Only operate on VCF records without this value present in FILTER. May be passed multiple times.",
         cxxopts::value<std::vector<std::string>>()) // multiple allowed
        ("I,include-reads",
         "SAM flag: only include reads with all of these bits set. Set 0 to disable. Default: 0",  // 3 requires reads to be paired + mapped in proper pair
         cxxopts::value<int>())
        ("E,exclude-reads",
         "SAM flag: exclude reads with any of these bits set. Default: 3844",  // 3852 also excludes unmapped mate
         cxxopts::value<int>())

        ("t,tsv",
         "Write a tsv of extended statistics to file specified.",
         cxxopts::value<fs::path>())
        ("u,uncompressed", "output uncompressed VCF")
        ("seed",
        "Set random seed. Default: 24601",
        cxxopts::value<uint32_t>())
        ("uniform",
        "additionally simulate against uniform null model for query position, and add result to --tsv output. For assessment of correlation with simulation against all-reads null.")
        ("assess-microhomology",
        "additionally assess STR and homopolymer content of reference regions, and add result to --tsv output. For assessment of correlation with drop in LZ.");

    options.add_options("") ("debug", "Enable debug logging", cxxopts::value<bool>()->default_value("false"));
    // clang-format on

    options.parse_positional ({"vcf", "aln"});
    options.positional_help ("<VCF/BCF (- for stdin)> <ALN.(b/cr)am>");

    return options;
}


void parse_cli
(ExposArgs& ctx, const cxo::ParseResult& input_args) {

  if (!input_args.count ("vcf") || !input_args.count ("aln"))
      throw std::runtime_error (
          "All positional arguments must be provided"
      );

  ctx.vcf_path = input_args["vcf"].as<fs::path>();
  ctx.aln_path = input_args["aln"].as<fs::path>();

  if (ctx.vcf_path.string() != "-" && !fs::exists (ctx.vcf_path)) {
      throw std::runtime_error (
          "VCF file not found: " + ctx.vcf_path.string()
      );
  }

  if (!fs::exists (ctx.aln_path)) {
      throw std::runtime_error (
          "Alignment file not found: " + ctx.aln_path.string()
      );
  }

  std::cerr << "Using VCF: " << ctx.vcf_path << std::endl;
  std::cerr << "Using aln: " << ctx.aln_path << std::endl;
  if (input_args.count ("ref")) {
      ctx.ref_path = input_args["ref"].as<fs::path>();
      if (!fs::exists (ctx.ref_path)) {
          throw std::runtime_error (
              "Reference fasta not found: " + ctx.ref_path.string()
          );
      }
      std::cerr << "Using ref: " << ctx.ref_path << std::endl;
  }

  if (input_args.count ("expected-read-len")) {
      ctx.exp_read_len = input_args["expected-read-len"].as<size_t>();
  } else {
      std::cerr <<
      std::format("Read length not provided, assuming {}", ctx.exp_read_len)
      << std::endl;
  }

  if (input_args.count ("include-records")) {
      ctx.flt_inc = input_args["include-records"].as<std::vector<std::string>>();
  }
  if (input_args.count ("exclude-records")) {
      ctx.flt_exc = input_args["exclude-records"].as<std::vector<std::string>>();
  }
  if (input_args.count ("include-reads")) {
      ctx.flag_inc = input_args["include-reads"].as<int>();
  }
  if (input_args.count ("exclude-reads")) {
      ctx.flag_exc = input_args["exclude-reads"].as<int>();
  }
  if (input_args.count ("seed")) {
      ctx.seed = input_args["seed"].as<uint32_t>();
  }

  if (input_args.count ("tsv")) {
      ctx.otsv_path = input_args["tsv"].as<fs::path>();
  }

  if (input_args.count ("normal")) {
      ctx.norm_paths = input_args["normal"].as<std::vector<fs::path>>();
      for (const auto &p : ctx.norm_paths) {
          if (!fs::exists (p))
              throw std::runtime_error ("Normal file not found: " + p.string());
          std::cerr << "Using normal: " << p << std::endl;
      }
  }
  if (input_args.count ("normal-only")) {
      if (ctx.norm_paths.empty())
          throw std::runtime_error("a normal must be provided if normal-only is set.");
      std::cerr << "Using only normal data as background for simulation" << std::endl;
      ctx.normal_only = true;
  }
  if (input_args.count ("uncompressed")) {
      ctx.no_gz = true;
  }
  if (input_args.count ("uniform")) {
      ctx.uniform_sim = true;
  }
  if (input_args.count ("assess-microhomology")) {
      ctx.assess_microhom = true;
  }
  if (input_args["debug"].as<bool>()) {
      ctx.debug_mode = true;
  }
}


}
