#pragma once

// Declarations for compute_info_field.cpp internals kept out of
// compute_info_field.hpp's public surface, exposed only for testing.
// Include from compute_info_field.cpp and its test suite only.

#include <cstdint>
#include <span>

#include "expos/pileup_features.hpp"

uint64_t count_pairs_within_1d (
    std::span<int32_t> pts, uint64_t radius
);

// Because each statistic is a pair-sum, the null's mean and variance under a
// uniform size-n draw without replacement are available in closed form and
// no simulation is needed for effect size.
struct KernelTerms {
  double t1;  // kernel summed over every background pair
  double t2;  // kernel value squaring each pair's value first
  double q;  // sum of the squared per-observation row sums
};

struct NullMoments {
  double mean;
  double var;
};

// derive exact mean and variance from kernel terms using inclusion
// probabilities at the draw size nObs.
NullMoments null_moments (
    uint64_t nObs, const KernelTerms& k, uint64_t nBackground
);

// Sorts pts in place. The kernel is an indicator, so every pair value is 0
// or 1 and squaring changes nothing: t2 == t1.
KernelTerms qpos_null_terms (
    std::span<int32_t> pts, uint64_t radius
);

KernelTerms jaccard_null_terms (
    std::span<const TemplateEndpoints> bg
);
