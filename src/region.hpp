// Copyright 2025 (c) Alex Byrne (alex@blex.bio), Luca Barbon; CASM
// Informatics, Wellcome Sanger Institute. All rights reserved. Use of
// this source code is governed by the MIT license that can be found
// in the LICENSE file.

#pragma once

#include <cstdint>
#include <htslib/hts.h>
#include <limits>
#include <string>

#include "bounds.hpp"


struct hts_region {
    int32_t rid;
    int64_t start;
    int64_t end;
    size_t rlen;

    // default invalid constructor
    hts_region () noexcept
        : rid (-1),
          start (-1),
          end (-1),
          rlen (0) {}

    static hts_region by_end (int32_t rid_,
                              int64_t gstart_,
                              int64_t gend_) {
        hts_region r;
        r.rid = rid_;
        r.start = gstart_;
        r.end = gend_;
        if (!r.valid_rid() || !r.valid_span())
            throw std::invalid_argument (
                "hts_region::by_end invalid parameters\nrid " +
                std::to_string (r.rid) + "\ngstart " +
                std::to_string (r.start) + "\ngend " +
                std::to_string (r.end));
        safe_size_opts sso;
        sso.msg = "hts_region::by_end - span too large " +
            std::to_string (r.start) + " " + std::to_string (r.end);
        r.rlen = safe_size (r.end - r.start, sso);

        return r;
    }

    static hts_region by_len (int32_t rid_,
                              int64_t gstart_,
                              size_t rlen_) {
        hts_region r;
        r.rid = rid_;
        r.start = gstart_;
        r.rlen = rlen_;
        if (!r.valid_rid() || !r.valid_rlen())
            throw std::invalid_argument (
                "hts_region::by_end invalid parameters\nrid " +
                std::to_string (r.rid) + "\ngstart " +
                std::to_string (r.start) + "\rlen" +
                std::to_string (r.rlen));
        r.end = r.start + static_cast<int64_t> (r.rlen);
        return r;
    }

    bool valid_rid () const noexcept { return rid >= 0; }
    bool valid_span () const noexcept { return end > start; }
    bool valid_rlen () const noexcept {
        return (std::numeric_limits<int64_t>::max() -
                    static_cast<int64_t> (rlen) >
                start);
    }
};

