#pragma once

// -- test-only declarations --- //

#include <cstddef>
#include <cstdint>
#include <span>

std::size_t count_pairs_within_1d (
    std::span<int32_t> pts, uint64_t radius
);
