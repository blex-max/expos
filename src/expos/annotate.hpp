#pragma once

#include <cstdint>
#include <vector>

#include "hts/hts_types.hpp"
#include "shared/err.hpp"

// Default Monte-Carlo seed (overridable with --seed).
inline constexpr std::uint32_t DEFAULT_SEED = 24601;

// State for one annotation run: I/O VCFs, sample + background alignments,
// reference, and run options. Built by init_ctx, consumed by analyse_records.
struct ExposCtx {
  VcfFile vcfIn;
  VcfFile vcfOut;
  AlnFile aln;
  std::vector<AlnFile> backgrounds;  // extra MC-background samples (--bg)
  FastaFile ref;
  std::uint32_t seed = DEFAULT_SEED;
  bool quiet = false;
  bool skipFiltered = false;
};

// Annotate every analysable record from ctx.vcfIn into ctx.vcfOut in place;
// unanalysable records pass through unchanged.
[[nodiscard]] VoidOrErr analyse_records (const ExposCtx& ctx);
