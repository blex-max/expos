#include "variant_stats.hpp"

#include <htslib/vcf.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <string_view>

#include "shared/stats.hpp"

// --- helpers --- //

static uint64_t manhattan (
    const TemplateEndpoints& a, const TemplateEndpoints& b
)
{
  const int64_t d1 =
      a.first > b.first ? a.first - b.first : b.first - a.first;
  const int64_t d2 =
      a.second > b.second ? a.second - b.second : b.second - a.second;
  return static_cast<uint64_t> (d1) + static_cast<uint64_t> (d2);
}

// Split the MC result into two subfields. The effect size is missing (with a
// zero-variance reason) when the null had no spread; the p-value is always
// present.
static std::vector<StatValue> mc_to_fields (const McResult& mc)
{
  return {
      mc.effectSize ? stat_value (*mc.effectSize)
                    : stat_missing (REASON_ZERO_VARIANCE),
      stat_value (mc.pValue)
  };
}

// Shared Monte-Carlo driver for a Ripley's-K-vs-background statistic. The
// intensity normalisation cancels in both the z-score and the p-value (all
// draws share the observed sample size), so only a within-radius pair count
// is needed here — no intensity term.
template <typename Point, typename CountWithin>
static McResult ripley_k_vs_background (
    const std::vector<Point>& observed,
    const std::vector<Point>& background, uint64_t radius,
    CountWithin countWithin, std::mt19937& rng
)
{
  const double observedStat =
      static_cast<double> (countWithin (observed, radius));
  auto draw = [&] {
    return subsample_wo_replace (background, observed.size(), rng);
  };
  auto stat = [&] (const std::vector<Point>& sample) {
    return static_cast<double> (countWithin (sample, radius));
  };
  return monte_carlo_pvalue (observedStat, draw, stat);
}

// A Ripley statistic needs >= 2 supporting observations and at least twice as
// many background reads. Returns the skip reason if a guard fails.
static std::optional<std::string_view> ripley_guard_reason (
    std::size_t nObs, std::size_t nBackground
)
{
  if (nObs < 2) {
    return REASON_INSUFFICIENT_SUPPORT;
  }
  if (nBackground < 2 * nObs) {
    return REASON_INSUFFICIENT_BACKGROUND;
  }
  return std::nullopt;
}

// --- statistics --- //

static std::vector<StatValue> compute_qrk (const VariantStatInputs& in)
{
  constexpr uint64_t k_qposRadius = 5;
  const auto& obs = in.supporting.qPos;
  const auto& bg = in.all.qPos;
  if (const auto reason = ripley_guard_reason (obs.size(), bg.size())) {
    return {stat_missing (*reason), stat_missing (*reason)};
  }
  // Query position is an offset within the read, so a heterogeneous read
  // population confounds this statistic specifically (unlike TRK).
  if (!in.readLenHomogeneous) {
    return {stat_missing (REASON_HETEROGENEOUS_READ_LENGTH),
            stat_missing (REASON_HETEROGENEOUS_READ_LENGTH)};
  }
  const auto mc = ripley_k_vs_background (
      obs, bg, k_qposRadius,
      [] (const std::vector<int32_t>& sample, uint64_t radius) {
        return count_pairs_within_1d (sample, radius);
      },
      in.rng
  );
  return mc_to_fields (mc);
}

static std::vector<StatValue> compute_trk (const VariantStatInputs& in)
{
  constexpr uint64_t k_templRadius = 5;
  const auto& obs = in.supporting.endpoints;
  const auto& bg = in.all.endpoints;
  if (const auto reason = ripley_guard_reason (obs.size(), bg.size())) {
    return {stat_missing (*reason), stat_missing (*reason)};
  }
  // Template endpoints track fragment boundaries, not read length, so this
  // statistic is not gated on read-length homogeneity.
  const auto mc = ripley_k_vs_background (
      obs, bg, k_templRadius,
      [] (const std::vector<TemplateEndpoints>& sample, uint64_t radius) {
        return count_pairs_within (sample, radius, manhattan);
      },
      in.rng
  );
  return mc_to_fields (mc);
}

static std::vector<StatValue> compute_mlas (const VariantStatInputs& in)
{
  constexpr double k_median = 0.5;
  return {
      stat_or (
          percentile (in.supporting.normalisedAs, k_median), REASON_NO_SUPPORT
      ),
      stat_or (
          percentile (in.all.normalisedAs, k_median), REASON_NO_BACKGROUND
      )
  };
}

static std::vector<StatValue> compute_rcmplx (const VariantStatInputs& in)
{
  constexpr std::size_t k_windowSize = 100;
  constexpr double k_scale = 100.0;
  const std::string_view ref = in.refSlice;
  // too-short spans have no full window; masked/ambiguous bases make the
  // complexity meaningless.
  if (ref.size() < k_windowSize) {
    return {stat_missing (REASON_REFERENCE_TOO_SHORT)};
  }
  if (ref.find_first_of ("Nn") != std::string_view::npos) {
    return {stat_missing (REASON_REFERENCE_HAS_N)};
  }
  double entropySum = 0.0;
  std::size_t nWin = 0;
  for (; nWin + k_windowSize <= ref.size(); ++nWin) {  // step of 1
    entropySum += entropy_lz76 (ref.substr (nWin, k_windowSize));
  }
  const double meanWindowEntropy =
      entropySum / static_cast<double> (nWin);
  return {stat_value (std::round (meanWindowEntropy * k_scale))};
}

// --- INFO header definitions --- //

constexpr std::string_view QRK_HEADER =
    "##INFO=<ID=QRK,Number=2,Type=Float,Description=\"Monte-Carlo "
    "results for spatial clustering of mutant query positions: [1] "
    "standardised effect size (z-score) versus simulation against all "
    "reads; [2] one-sided p-value. Effect sizes greater than ~2.0 with a "
    "significant p-value may indicate a spurious variant.\">";
constexpr std::string_view TRK_HEADER =
    "##INFO=<ID=TRK,Number=2,Type=Float,Description=\"Monte-Carlo "
    "results for spatial clustering of mutant template endpoints: [1] "
    "standardised effect size (z-score) versus simulation against all "
    "reads; [2] one-sided p-value. Effect sizes greater than ~2.0 with a "
    "significant p-value may indicate a spurious variant.\">";
constexpr std::string_view MLAS_HEADER =
    "##INFO=<ID=MLAS,Number=2,Type=Float,Description=\"Median "
    "read-length-normalised alignment scores: [1] of reads supporting the "
    "variant; [2] of all reads covering the variant site.\">";
constexpr std::string_view RCMPLX_HEADER =
    "##INFO=<ID=RCMPLX,Number=1,Type=Float,Description=\"Mean 100-base "
    "window complexity (Lempel-Ziv 76 entropy rate) of the reference "
    "region spanned by supporting templates, scaled by 100.\">";
constexpr std::string_view EXPOS_SKIP_HEADER =
    "##INFO=<ID=EXPOS_SKIP,Number=.,Type=String,Description=\"Why expos "
    "produced no value, as '<scope>:<reason>' tokens. Scope is 'record' "
    "(whole record skipped) or a statistic ID. Reasons: multiallelic, "
    "complex, insufficient_support, insufficient_background, "
    "heterogeneous_read_length, zero_variance, no_support, no_background, "
    "reference_too_short, reference_has_n.\">";

// --- registry --- //

// The full statistics package, fixed at compile time (no runtime selection).
constexpr std::array<VariantStat, 4> VARIANT_STATS = {{
    {{"QRK", QRK_HEADER, 2}, &compute_qrk},
    {{"TRK", TRK_HEADER, 2}, &compute_trk},
    {{"MLAS", MLAS_HEADER, 2}, &compute_mlas},
    {{"RCMPLX", RCMPLX_HEADER, 1}, &compute_rcmplx},
}};

std::span<const VariantStat> variant_stats ()
{
  return VARIANT_STATS;
}

// --- VCF output --- //

VoidOrErr register_variant_stat_header (
    bcf_hdr_t* hdr, const VariantStat& stat
)
{
  const std::string line{stat.field.headerLine};
  if (bcf_hdr_append (hdr, line.c_str()) != 0) {
    return std::unexpected (make_err (
        "failed to add INFO header line for " +
        std::string (stat.field.id)
    ));
  }
  return {};
}

VoidOrErr register_expos_skip_header (bcf_hdr_t* hdr)
{
  const std::string line{EXPOS_SKIP_HEADER};
  if (bcf_hdr_append (hdr, line.c_str()) != 0) {
    return std::unexpected (
        make_err ("failed to add EXPOS_SKIP INFO header line")
    );
  }
  return {};
}

VoidOrErr encode_variant_stat (
    bcf_hdr_t* hdr, bcf1_t* rec, const VariantStatField& field,
    const std::vector<StatValue>& values
)
{
  assert (static_cast<int> (values.size()) == field.nValues);

  // vcf works at float precision
  std::vector<float> buf (values.size());
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (values[i].value) {
      buf[i] = static_cast<float> (*values[i].value);
    }
    else {
      bcf_float_set_missing (buf[i]);
    }
  }

  const std::string id{field.id};
  if (bcf_update_info_float (
          hdr, rec, id.c_str(), buf.data(),
          static_cast<int> (buf.size())
      ) != 0) {
    return std::unexpected (
        make_err ("failed to write INFO field " + std::string (field.id))
    );
  }
  return {};
}

std::vector<std::string> stat_skip_tokens (
    const VariantStatField& field, const std::vector<StatValue>& values
)
{
  std::vector<std::string> tokens;
  for (const auto& v : values) {
    if (!v.reason) {
      continue;
    }
    std::string token =
        std::string (field.id) + ":" + std::string (*v.reason);
    if (std::find (tokens.begin(), tokens.end(), token) == tokens.end()) {
      tokens.push_back (std::move (token));
    }
  }
  return tokens;
}

VoidOrErr set_expos_skip (
    bcf_hdr_t* hdr, bcf1_t* rec, const std::vector<std::string>& tokens
)
{
  if (tokens.empty()) {
    return {};
  }
  std::string joined;
  for (std::size_t i = 0; i < tokens.size(); ++i) {
    if (i != 0) {
      joined += ',';
    }
    joined += tokens[i];
  }
  if (bcf_update_info_string (hdr, rec, "EXPOS_SKIP", joined.c_str()) != 0) {
    return std::unexpected (
        make_err ("failed to write EXPOS_SKIP INFO field")
    );
  }
  return {};
}
