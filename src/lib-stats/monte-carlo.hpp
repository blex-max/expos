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
        std::invoke_result_t<StatFn&, const std::vector<ObsT>&>,
        StatT
    > &&
    std::same_as<
        std::invoke_result_t<EffFn&, StatT, const std::vector<StatT>&>,
        std::optional<double>
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

    size_t sim_count_gt = 0;
    for (size_t i = 0; i < nsim; ++i) {
        sim_obs.clear();
        sim_obs = drawfn();

        const auto draw_stat = statfn (sim_obs);
        if (draw_stat >= ev_stat) {
            ++sim_count_gt;
        }
        draw_results.push_back (draw_stat);
    }

    // NOTE: effsz and pval calcuations have changed!
    // these changes have not been propagated through
    // the rest of the codebase or documentation

    // why does this take a fn
    res.eff_sz = efffn (ev_stat, draw_results);

    // TODO "power analysis"

    // one sided
    res.pval = static_cast<double> (sim_count_gt + 1)
               / static_cast<double> (nsim + 1);

    return res;
}


// standardised effect size
// z-score relative to the simulated null distribution
inline std::optional<double> ses
(const double& testv, const std::vector<double>& simv)
{
    if (simv.empty()) {
        return {};
    }

    // mean of simulated results
    const auto sim_mean = *summary::mean (simv);
    // calc pop stddev
    const auto pn = simv.size();
    double idev = 0.0;
    double devsum = 0.0;
    for (size_t i = 0; i < pn; ++i) {
        idev = simv[i] - sim_mean;
        devsum += idev * idev;
    }
    const auto sim_sd =
        std::sqrt (devsum / static_cast<double> (pn));

    if (sim_sd == 0.0) {  // unlikely but possible
        return {};  // otherwise div by zero
    }

    return (testv - sim_mean) / sim_sd;
}


// get a random sample without replacement of
// size n from input obs vector using a partial
// Fisher-Yates shuffle. Faster than Floyd
// sampling when n is not << than Nobs -
// the expected domain for low yield sequencing.
template <typename T>
inline std::vector<T> subsample_wo_replace (
    const std::vector<T>& obs,
    size_t n,
    std::mt19937 &rng
) {
    std::vector<T> out;
    size_t nobs = obs.size();
    assert (n < nobs);

    // partial Fisher–Yates shuffle an index vector so that
    // idx[0..n] is a uniform sample without replacement.
    std::vector<size_t> all_idx(nobs);
    std::iota(begin(all_idx), end(all_idx), 0);
    for (size_t i = 0; i < n; ++i) {
        std::uniform_int_distribution<size_t> dist(i, nobs - 1);
        size_t j = dist(rng);
        std::swap(all_idx[i], all_idx[j]);
        out.push_back(obs[all_idx[i]]);
    }

    return out;
}

}  // end namespace

