#pragma once

#include <cassert>
#include <optional>
#include <random>
#include <vector>

#include "lib-stats/summary.hpp"

namespace monte_carlo {

struct stat_eval_s {
    std::optional<double> eff_sz = std::nullopt;
    std::optional<double> pval   = std::nullopt;
    std::string           err    = "";
};
// are event supporting observations
// meaningfully different than total observations
// for a given statistic?
// compare via simulation: if we take random samples
// of the total observations of the same size
// as the true event supporting sample,
// how often do we get a lt/gt value for the
// statistic in question (pvalue), and
// how large is the effect size.
template <
    typename StatT,
    typename StatFn,
    typename DrawFn,
    typename EffFn,
    typename ObsT = typename std::invoke_result_t<DrawFn&>::value_type
>
requires
    std::same_as<
        std::invoke_result_t<DrawFn&>,
        std::vector<ObsT>
    > &&
    std::same_as<
        std::invoke_result_t<StatFn&, const std::vector<ObsT> &>,
        StatT
    > &&
    std::same_as<
        std::invoke_result_t<EffFn&, StatT, const std::vector<StatT> &>,
        double
    >
inline stat_eval_s sim_to_bg (
    StatT    ev_stat,
    DrawFn &&drawfn,
    StatFn &&statfn,
    EffFn  &&efffn,
    size_t   nsim=2500
) {
    stat_eval_s res;
    assert (nsim > 0);
    std::vector<ObsT> sim_obs;

    std::vector<StatT> draw_results;
    draw_results.reserve(nsim);

    size_t sim_count_le = 0;  // count of simulated stats that are <= observed stat
    for (size_t i = 0; i < nsim; ++i) {
        sim_obs.clear();
        sim_obs = drawfn();

        const auto draw_stat = statfn (sim_obs);
        if (draw_stat <= ev_stat) {
            ++sim_count_le;
        }
        draw_results.push_back (draw_stat);
    }

    // report effect size
    // if eff_sz is large then we can
    // get away with a low number of samples
    // if not it's just noise
    res.eff_sz = efffn (ev_stat, draw_results);

    // TODO "power analysis"

    // two sided p val
    res.pval = 2
               * static_cast<double> (std::min(sim_count_le + 1, nsim - sim_count_le + 1))
               / static_cast<double> (nsim + 1);

    // one sided
    // res.pval = static_cast<double> (sim_count_le + 1)
    //            / static_cast<double> (nsim + 1);
    return res;
}


// log2-fold change
// requires strictly postive values!
// 
// +1 removes confusing values when 0,
// log makes effect size symmetric around 0
// log2 means -1 = half the size of background
// +1 = double the size of background
inline auto log2_effsz (const double &ev, const std::vector<double> &simv) {
    return log2 ((ev + 1) / (*summary::mean (simv) + 1));
}


// get a random sample without replacement of
// size n from input obs vector.
template <typename T>
inline std::vector<T> subsample_wo_replace (
    const std::vector<T>& obs,
    size_t n,
    std::mt19937 &rng
) {
    size_t nobs = obs.size();
    assert (n < nobs);

    // partial Fisher–Yates shuffle an index vector so that
    // idx[0..n] is a uniform sample without replacement.
    std::vector<size_t> all_idx(nobs);
    std::iota(begin(all_idx), end(all_idx), 0);

    // shuffle
    for (size_t i = 0; i < n; ++i) {
        std::uniform_int_distribution<size_t> dist(i, nobs - 1);
        size_t j = dist(rng);
        std::swap(all_idx[i], all_idx[j]);
    }

    // get obs
    std::vector<T> out;
    for (size_t k = 0; k < n; ++k) {
        out.push_back(obs[all_idx[k]]);
    }
    return out;
}

}  // end namespace

