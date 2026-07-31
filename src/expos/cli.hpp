#pragma once

#include <cxxopts.hpp>
#include <filesystem>

namespace cli {

namespace cxo = cxxopts;
namespace fs = std::filesystem;

cxo::Options setup_cli();

struct ExposArgs {
  fs::path vcf_path;
  fs::path aln_path;
  std::vector<fs::path> norm_paths;
  fs::path ref_path;
  fs::path otsv_path;
  std::vector<std::string> flt_inc;
  std::vector<std::string> flt_exc;
  uint32_t seed = 24601;
  size_t exp_read_len = 150;
  int flag_inc = 0;
  int flag_exc = 3844;
  bool no_gz = false;
  bool normal_only = false;
  bool uniform_sim = false;
  bool assess_microhom = false;
  bool debug_mode = false;
};

void parse_cli (ExposArgs& ctx, const cxo::ParseResult& input_args);

}  // namespace cli
