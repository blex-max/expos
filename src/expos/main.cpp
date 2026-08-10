#include <fmt/chrono.h>
#include <fmt/format.h>
#include <htslib/vcf.h>
#include <plog/Formatters/TxtFormatter.h>
#include <plog/Initializers/ConsoleInitializer.h>
#include <plog/Log.h>

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
  std::uint32_t seed = DEFAULT_SEED;
  bool uncompressed = false;
  bool debug = false;
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
      .scan<'u', std::uint32_t>()
      .help ("random seed for the Monte-Carlo simulation");
  cli.add_argument ("-u", "--uncompressed")
      .flag()
      .help (
          "write uncompressed VCF (default: bgzip-compressed)"
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
  cli.add_argument ("--additional-background-samples", "--bg")
      .default_value (std::vector<std::string>{})
      .nargs (argparse::nargs_pattern::at_least_one)
      .metavar ("PATH")
      .help (
          "additional indexed alignment file(s) whose reads are "
          "merged into "
          "the Monte-Carlo background. Supporting reads are "
          "always taken from the primary ALN only."
      );
  cli.add_argument ("--debug").flag().help (
      "enable debug logging to stderr"
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
      .seed = cli.get<std::uint32_t> ("--seed"),
      .uncompressed = cli.get<bool> ("--uncompressed"),
      .debug = cli.get<bool> ("--debug"),
      .quiet = cli.get<bool> ("--quiet"),
      .skipFiltered = cli.get<bool> ("--skip-filtered"),
      .bgPaths = cli.get<std::vector<std::string>> (
          "--additional-background-samples"
      ),
      .invocation = std::move (invocation)
  };
}

static VcfOrErr create_output_vcf (
    const VcfFile& input, bool uncompressed,
    const std::string& invocation
)
{
  VcfFile out;
  out.fh = hts_open ("-", uncompressed ? "w" : "wz");
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
static CtxOrErr init_ctx (const ExposArgs& filepaths)
{
  ExposCtx out;

  auto vcfRet = load_vcf (filepaths.vcfPath);
  if (!vcfRet) {
    return std::unexpected (vcfRet.error());
  }
  out.vcfIn = std::move (*vcfRet);

  auto alnRet = load_aln (filepaths.alnPath);
  if (!alnRet) {
    return std::unexpected (alnRet.error());
  }
  out.aln.handle = std::move (*alnRet);
  // takes ownership of PileupContext
  out.aln.plpIt = init_pileup_iterator (
      new PileupContext{out.aln.handle}, pileup_fn
  );

  out.backgrounds.reserve (filepaths.bgPaths.size());
  for (const auto& bgPath : filepaths.bgPaths) {
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

  auto refRet = load_fasta (filepaths.refPath);
  if (!refRet) {
    return std::unexpected (refRet.error());
  }
  out.ref = std::move (*refRet);

  auto vcfOutRet = create_output_vcf (
      out.vcfIn, filepaths.uncompressed, filepaths.invocation
  );
  if (!vcfOutRet) {
    return std::unexpected (vcfOutRet.error());
  }
  out.vcfOut = std::move (*vcfOutRet);

  out.seed = filepaths.seed;
  out.quiet = filepaths.quiet;
  out.skipFiltered = filepaths.skipFiltered;

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

  if (args.debug) {
    plog::init<plog::TxtFormatter> (
        plog::debug, plog::streamStdErr
    );
  }

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
