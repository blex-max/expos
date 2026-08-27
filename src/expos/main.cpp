#include <fmt/chrono.h>
#include <fmt/format.h>
#include <htslib/vcf.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <expected>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "argparse/argparse.hpp"
#include "expos/annotate.hpp"
#include "expos/compute_info_field.hpp"
#include "expos/encode_info_field.hpp"
#include "expos/pileup_fn.hpp"
#include "hts/hts_types.hpp"
#include "shared/err.hpp"

// Defined in CMakeLists.txt
#ifndef EXPOS_VERSION
#define EXPOS_VERSION "unknown"
#endif
inline constexpr std::string_view PROG_NAME = "expos";
static constexpr char VERSION[] = EXPOS_VERSION;

struct ExposArgs {
  std::string vcfPath;
  std::string alnPath;
  std::string refPath;
  uint32_t seed = DEFAULT_SEED;
  uint16_t flankSize = DEFAULT_FLANK;
  bool quiet = false;
  bool skipFiltered = false;
  std::vector<std::string> bgPaths;
  std::string invocation;  // CLI call
};

using ArgsOrErr = std::expected<ExposArgs, Err>;
static ArgsOrErr parse_args (int argc, char** argv)
{
  // Command line for the provenance header; handling single quotes
  std::string invocation;
  for (int i = 0; i < argc; ++i) {
    if (i != 0) {
      invocation += ' ';
    }
    invocation += argv[i];
  }
  std::replace (invocation.begin(), invocation.end(), '"', '\'');

  auto cli = argparse::ArgumentParser ("expos");
  cli.set_usage_max_line_width (80);
  cli.add_argument ("VCF").help (
      "input VCF/BCF of variants to annotate"
  );
  cli.add_argument ("REF").help (
      "indexed reference genome FASTA"
  );
  cli.add_argument ("ALN").help (
      "indexed alignment (BAM/CRAM) of the sample"
  );
  cli.add_argument ("--seed")
      .default_value (DEFAULT_SEED)
      .nargs (1)
      .metavar ("SEED")
      .scan<'u', uint32_t>()
      .help ("random seed for the Monte-Carlo simulation");
  cli.add_argument ("--flank")
      .default_value (DEFAULT_FLANK)
      .nargs (1)
      .metavar ("SIZE")
      .scan<'u', uint16_t>()
      .help (
          "Size of reference sequence flanks to retrieve from\n"
          "either side of variant, for use in calcuating\n"
          "reference complexity. It is suggested to set to\n"
          "approximately the average template size for the\n"
          "sequencing protocol used."
      );
  cli.add_argument ("-q", "--quiet")
      .flag()
      .help ("suppress per-record warnings to stderr");
  cli.add_argument ("--skip-filtered")
      .flag()
      .help (
          "Only analyse records where FILTER is PASS or . "
          "(unset)"
      );
  cli.add_argument ("-b", "--background-sample")
      .default_value (std::vector<std::string>{})
      .append()
      .nargs (1)
      .metavar ("PATH")
      .help (
          "additional indexed alignment file/s from which\n"
          "to merge data into Monte-Carlo background.\n"
          "Supporting reads are always taken from the\n"
          "primary ALN only."
      );

  try {
    cli.parse_args (argc, argv);
  }
  catch (const std::exception& ex) {
    std::ostringstream oss;
    oss << ex.what() << "\n" << cli;
    return std::unexpected (make_err (oss.str()));
  }

  return ExposArgs{
      .vcfPath = cli.get<std::string> ("VCF"),
      .alnPath = cli.get<std::string> ("ALN"),
      .refPath = cli.get<std::string> ("REF"),
      .seed = cli.get<uint32_t> ("--seed"),
      .flankSize = cli.get<uint16_t> ("--flank"),
      .quiet = cli.get<bool> ("--quiet"),
      .skipFiltered = cli.get<bool> ("--skip-filtered"),
      .bgPaths = cli.get<std::vector<std::string>> (
          "--background-sample"
      ),
      .invocation = std::move (invocation)
  };
}

static VcfOrErr create_output_vcf (
    const VcfFile& input, const std::string& invocation
)
{
  VcfFile out;
  out.fh = hts_open ("-", "w");
  if (out.fh == nullptr) {
    return std::unexpected (
        make_err ("Could not open stdout for writing output VCF")
    );
  }
  out.hdr = bcf_hdr_dup (input.hdr);
  if (out.hdr == nullptr) {
    return std::unexpected (
        make_err ("Could not copy header to output VCF")
    );
  }

  // Register INFO fields, write header
  for (const auto& stat : expos_field_registry()) {
    auto regRet = register_variant_stat_header (out.hdr, stat);
    if (!regRet) {
      return std::unexpected (regRet.error());
    }
  }
  if (auto regRet = register_expos_skip_header (out.hdr);
      !regRet) {
    return std::unexpected (regRet.error());
  }

  // Write provenance
  const auto now = std::chrono::floor<std::chrono::seconds> (
      std::chrono::system_clock::now()
  );
  const std::string sourceLine = fmt::format (
      "##source=\"{} v{} {:%Y-%m-%dT%H:%M:%SZ}\"", PROG_NAME,
      VERSION, now
  );
  if (bcf_hdr_append (out.hdr, sourceLine.c_str()) != 0) {
    return std::unexpected (
        make_err ("Could not add source header line")
    );
  }
  const std::string cmdLine =
      fmt::format ("##{}_command=\"{}\"", PROG_NAME, invocation);
  if (bcf_hdr_append (out.hdr, cmdLine.c_str()) != 0) {
    return std::unexpected (
        make_err ("Could not add command header line")
    );
  }

  if (bcf_hdr_write (out.fh, out.hdr) != 0) {
    return std::unexpected (
        make_err ("Could not write output VCF header")
    );
  }

  return out;
}

using CtxOrErr = std::expected<ExposCtx, Err>;
static CtxOrErr init_ctx (const ExposArgs& args)
{
  ExposCtx out;

  auto vcfRet = load_vcf (args.vcfPath);
  if (!vcfRet) {
    return std::unexpected (vcfRet.error());
  }
  out.vcfIn = std::move (*vcfRet);

  auto alnRet = load_aln (args.alnPath);
  if (!alnRet) {
    return std::unexpected (alnRet.error());
  }
  out.aln.handle = std::move (*alnRet);
  // takes ownership of PileupContext
  out.aln.plpIt = init_pileup_iterator (
      new PileupContext{out.aln.handle}, pileup_fn
  );

  out.backgrounds.reserve (args.bgPaths.size());
  for (const auto& bgPath : args.bgPaths) {
    auto bgRet = load_aln (bgPath);
    if (!bgRet) {
      return std::unexpected (bgRet.error());
    }
    // too many moves...
    // Could the initialiser be templated
    // on the derived pileupcontext type and
    // take arguments to construct in place?
    AlnBundle bundle;
    bundle.handle = std::move (*bgRet);
    bundle.plpIt = init_pileup_iterator (
        new PileupContext{bundle.handle}, pileup_fn
    );
    out.backgrounds.emplace_back (std::move (bundle));
  }

  auto refRet = load_fasta (args.refPath);
  if (!refRet) {
    return std::unexpected (refRet.error());
  }
  out.ref = std::move (*refRet);

  auto vcfOutRet =
      create_output_vcf (out.vcfIn, args.invocation);
  if (!vcfOutRet) {
    return std::unexpected (vcfOutRet.error());
  }
  out.vcfOut = std::move (*vcfOutRet);

  out.seed = args.seed;
  out.flankSize = args.flankSize;
  out.quiet = args.quiet;
  out.skipFiltered = args.skipFiltered;

  return out;
}

int main (int argc, char** argv)
{
  auto argRet = parse_args (argc, argv);
  if (!argRet) {
    std::cerr << argRet.error().msg << std::endl;
    return EXIT_FAILURE;
  }
  const auto args = std::move (*argRet);

  auto initRet = init_ctx (args);
  if (!initRet) {
    std::cerr << initRet.error().msg << std::endl;
    return EXIT_FAILURE;
  }
  auto ctx = std::move (*initRet);

  std::cerr << "analysing records" << std::endl;

  auto analyseRet = analyse_records (ctx);
  if (!analyseRet) {
    std::cerr << analyseRet.error().msg << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
