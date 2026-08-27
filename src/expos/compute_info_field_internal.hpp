#pragma once

// Declarations for compute_info_field.cpp internals kept out of
// compute_info_field.hpp's public surface, exposed only for testing.
// Include from compute_info_field.cpp and its test suite only.

#include <cstdint>
#include <span>

uint64_t count_pairs_within_1d (
    std::span<int32_t> pts, uint64_t radius
);
