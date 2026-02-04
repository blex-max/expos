#pragma once

#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <math.h>
#include <optional>
#include <vector>


namespace spatial {

// object described by two
// coordinates on the same axis
struct line_seg {
    uint64_t lmost;
    uint64_t rmost;

    line_seg() = delete;
    line_seg (
        uint64_t a,
        uint64_t b
    ) {
        lmost = a > b ? b : a;
        rmost = a > b ? a : b;
    }

    uint64_t diff () const { return rmost - lmost; }
};


// 2D symmetric square matrix via vector
// rows are contiguous in vector
class PairMatrix {
  private:
    std::vector<uint64_t> mat;
    const size_t          dim_;

    PairMatrix (
        std::vector<uint64_t> v,
        size_t                dim
    )
        : mat (v),
          dim_ (dim) {
        assert (dim > 1);
        assert (v.size() == (dim * dim));
    }

  public:
    PairMatrix() = delete;

    auto dim () const noexcept { return dim_; }
    auto size () const noexcept { return mat.size(); }

    auto get (
        size_t i,
        size_t j
    ) const {
        assert (i < dim_);
        assert (j < dim_);
        assert (!mat.empty());
        assert (mat.size() == (dim() * dim()));
        return mat[(i * dim_) + j];
    }

    const auto &get1D () const noexcept { return mat; }
    auto        copy1D () const noexcept { return mat; }

    // assumes symmetric distance
    // clang-format off
    template <typename U, typename F>
        requires std::invocable<F &, const U &, const U &>
                 && std::same_as<std::invoke_result_t<F &,const U &,const U &>, uint64_t>
    // clang-format on
    static std::optional<PairMatrix> from_sample (
        const std::vector<U> &obs,
        F                   &&sym_pairfn
    ) {
        assert (!obs.empty());
        const auto dim = obs.size();
        if (dim < 2)
            return std::nullopt;
        const auto            nel = dim * dim;
        std::vector<uint64_t> in (nel);     // nel-long vector
        for (size_t i = 0; i < dim; ++i) {
            for (size_t j = 0; j < (i + 1); ++j) {
                const auto val    = sym_pairfn (obs[i], obs[j]);
                in[(i * dim) + j] = val;
                if (j != i)
                    in[(j * dim) + i] = val;
            }
        }
        return PairMatrix{in, dim};
    }
};


// no edge correction - fine if only comparing to monte carlo the same way
inline double ripley_k (const PairMatrix &pwd, uint64_t t, double point_intensity) {
    const auto n = pwd.dim();  // square matrix
    const auto intensity_factor = 1 / point_intensity;

    size_t lt_sum = 0;
    for (size_t row = 0; row < n; ++row) {
        size_t lt_count = 0;
        for (size_t col = 0; col < n; ++col) {
            if (row == col) continue;
            if (pwd.get(row, col) < t) {
                ++lt_count;
            }
        }
        lt_sum += lt_count;
    }

    const auto mean_lt = static_cast<double>(lt_sum) / static_cast<double>(n);
    return intensity_factor * mean_lt;
}


}  // end namespace
