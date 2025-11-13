// Copyright 2025 (c) Alex Byrne (alex@blex.bio), Luca Barbon; CASM
// Informatics, Wellcome Sanger Institute. All rights reserved. Use of
// this source code is governed by the MIT license that can be found
// in the LICENSE file.

#pragma once

#include <cstdint>
#include <limits>
#include <string>

struct safe_size_opts {
    size_t lower    = 0;
    size_t upper    = std::numeric_limits<size_t>::max();
    std::string msg = "";
};

inline size_t safe_size (int64_t i,
                         safe_size_opts opts = {}) {
    try {
        if (i < 0)
            throw std::out_of_range ("size would be negative");

        // cast fine since we know non negative now
        if (static_cast<uint64_t> (i) < opts.lower)
            throw std::out_of_range (
                "size would be below lower bound");

        if (static_cast<uint64_t> (i) > opts.upper)
            throw std::out_of_range ("size would exceed upper bound");

        // convert in peace
        return static_cast<size_t> (i);
    } catch (std::exception &e) {
        throw std::runtime_error (opts.msg + ": " + e.what());
    }
}
