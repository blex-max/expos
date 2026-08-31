#pragma once

#include <cstdint>
#include <vector>

#include "hts/hts_types.hpp"
#include "shared/err.hpp"

// Default Monte-Carlo seed.
inline constexpr uint32_t DEFAULT_SEED = 24601;
// Default reference flank size for complexity analysis.
inline constexpr uint16_t DEFAULT_FLANK = 250;
// Default maximum fragment length to be considered valid,
// to avoid inclusion of ambiguously mapped fragments with
// huge fragment lengths.
inline constexpr uint16_t DEFAULT_MAX_FRAG_LEN = 2000;

struct AlnBundle {
  AlnFile handle;  // manages alignment lifetime
  ForwardPileupIterator
      plpIt;  // borrows ptrs from handle, responsible for pileup objects only

  AlnBundle() = default;

  // no copy
  AlnBundle (const AlnBundle&) = delete;
  AlnBundle& operator= (const AlnBundle&) = delete;

  // move ok: plpIt's context caches handle's htsFile/index *pointers*, not
  // the address of handle itself, and moving an AlnFile leaves those
  // pointees untouched. A context holding &handle would not survive this.
  AlnBundle (AlnBundle&&) = default;
  AlnBundle& operator= (AlnBundle&&) = default;
};

// State for one annotation run: I/O VCFs, sample + background alignments,
// reference, and run options. Built by init_ctx, consumed by analyse_records.
struct ExposCtx {
  VcfFile vcfIn;
  VcfFile vcfOut;
  AlnBundle aln;
  std::vector<AlnBundle>
      backgrounds;  // extra MC-background samples (-b/--background-sample)
  FastaFile ref;
  uint32_t seed = DEFAULT_SEED;
  bool quiet = false;
  bool skipFiltered = false;
  uint16_t flankSize = DEFAULT_FLANK;
  uint16_t maxFragLen = DEFAULT_MAX_FRAG_LEN;
};

// Annotate every analysable record from ctx.vcfIn into ctx.vcfOut in place;
// unanalysable records pass through unchanged.
[[nodiscard]] VoidOrErr analyse_records (ExposCtx& ctx);
