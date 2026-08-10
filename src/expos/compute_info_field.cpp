#include "compute_info_field.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string_view>
#include <vector>

#include "expos/guards.hpp"
#include "expos/variant_stats.hpp"
#include "shared/stats.hpp"

// --- helpers --- //

// Builders keeping the compute_* functions readable.
static StatValue stat_value (double v)
{
  return {v, std::nullopt};
}

static StatValue stat_missing (StatSkipReason reason)
{
  return {std::nullopt, reason};
}

static StatValue stat_or (
    std::optional<double> v, StatSkipReason reasonIfMissing
)
{
  return v ? stat_value (*v) : stat_missing (reasonIfMissing);
}

static ValuesOrSkip mc_to_fields (const McResult& mc)
{
  if (!mc.effectSize) {
    return std::unexpected (StatSkipReason::zeroVariance);
  }
  return std::vector<StatValue>{
      stat_value (*mc.effectSize), stat_value (mc.pValue)
  };
}

void set_qrk_guard (
    StatContext& sCtx, const PileupFeatures& primaryAll
)
{
  if (!sufficient_reads (primaryAll.readLen.size())) {
    // fail-closed, but say so: too few reads to reach any verdict on
    // the spread is not the same claim as reads of uneven length.
    sCtx.readLenSuppression =
        StatSkipReason::readLengthUnverified;
    return;
  }
  sCtx.readLenSuppression =
      read_lens_within_tol (primaryAll.readLen)
          ? std::nullopt
          : std::optional{
                StatSkipReason::heterogeneousReadLength
            };
}

// --- statistics --- //

static ValuesOrSkip compute_qrk (
    const VariantStatInputs& in, const StatContext& ctx
)
{
  constexpr uint64_t QPOS_RADIUS = 5;
  const auto& obs = in.supporting.qPos;
  const auto& bg = in.all.qPos;

  // Size before quality: a locus with too few reads must report that it
  // was too thin, not a verdict on read-length spread that the guard
  // never had the data to reach.
  if (const auto reason = size_guard (obs.size(), bg.size())) {
    return std::unexpected (*reason);
  }

  // Query position is an offset within the read, so a heterogeneous read
  // population would confound
  if (ctx.readLenSuppression) {
    return std::unexpected (*ctx.readLenSuppression);
  }

  // count_pairs_within_1d sorts, so the observed statistic works on a copy
  // and the record's features stay const. One copy per record, not per draw.
  std::vector<int32_t> obsSorted = obs;
  const double observedStat = static_cast<double> (
      count_pairs_within_1d (obsSorted, QPOS_RADIUS)
  );

  auto draw = [&] {
    return subsample_wo_replace (
        bg, obs.size(), ctx.mc.rng, ctx.mc.qPos
    );
  };
  auto stat = [&] (std::span<int32_t> sample) {
    return static_cast<double> (
        count_pairs_within_1d (sample, QPOS_RADIUS)
    );
  };
  const auto mc = run_monte_carlo (observedStat, draw, stat);
  return mc_to_fields (mc);
}

static ValuesOrSkip compute_template_jaccard (
    const VariantStatInputs& in, const StatContext& ctx
)
{
  const auto& obs = in.supporting.endpoints;
  const auto& bg = in.all.endpoints;

  if (const auto reason = size_guard (obs.size(), bg.size())) {
    return std::unexpected (*reason);
  }

  // Template endpoints track fragment boundaries, not read length,
  // so not gated on read-length homogeneity.

  const auto obsJaccardSum = pairwise_jaccard_sum (obs);

  const auto drawFn = [&] {
    return subsample_wo_replace (
        bg, obs.size(), ctx.mc.rng, ctx.mc.endpoints
    );
  };

  const auto mc = run_monte_carlo (
      obsJaccardSum, drawFn, pairwise_jaccard_sum
  );

  return mc_to_fields (mc);
}

// The only statistic with independently-missing subfields: two unrelated
// medians, either of which can be absent while the other is reportable.
static ValuesOrSkip compute_mlas (
    const VariantStatInputs& in, const StatContext&
)
{
  return std::vector<StatValue>{
      stat_or (
          percentile (
              in.supporting.normalisedAs, PERCENTILE_MEDIAN
          ),
          StatSkipReason::noSupport
      ),
      stat_or (
          percentile (in.all.normalisedAs, PERCENTILE_MEDIAN),
          StatSkipReason::noBackground
      )
  };
}

static ValuesOrSkip compute_rcmplx (
    const VariantStatInputs& in, const StatContext&
)
{
  constexpr size_t WIN_SZ = 100;
  constexpr size_t WIN_STEP = 10;

  const std::string_view ref = in.refSlice;
  // too-short spans have no full window; masked/ambiguous bases make the
  // complexity meaningless.
  if (ref.size() < WIN_SZ) {
    return std::unexpected (StatSkipReason::referenceTooShort);
  }
  if (ref.find_first_of ("Nn") != std::string_view::npos) {
    return std::unexpected (StatSkipReason::referenceHasN);
  }

  double entropySum = 0.0;
  size_t nWin = 0;
  for (size_t winStart = 0; winStart + WIN_SZ <= ref.size();
       winStart += WIN_STEP) {
    entropySum += entropy_lz76 (ref.substr (winStart, WIN_SZ));
    ++nWin;
  }
  const double meanWindowEntropy =
      entropySum / static_cast<double> (nWin);
  return std::vector<StatValue>{stat_value (meanWindowEntropy)};
}

// --- INFO header definitions --- //

constexpr std::string_view QRK_HEADER =
    "##INFO=<ID=QRK,Number=2,Type=Float,Description=\"Monte-"
    "Carlo "
    "results for spatial clustering of mutant query positions: "
    "[0] "
    "standardised effect size (z-score) versus simulation "
    "against all "
    "reads; [1] one-sided p-value. Effect sizes greater than "
    "~3.0 with a "
    "significant p-value may indicate a spurious variant.\">";
constexpr std::string_view TJAC_HEADER =
    "##INFO=<ID=TJAC,Number=2,Type=Float,Description=\"Monte-"
    "Carlo "
    "results for graded pairwise overlap of mutant templates: "
    "[0] "
    "standardised effect size (z-score) versus simulation "
    "against all "
    "reads; [1] one-sided p-value. Effect sizes greater than "
    "~3.0 with a "
    "significant p-value may indicate a spurious variant.\">";
constexpr std::string_view MLAS_HEADER =
    "##INFO=<ID=MLAS,Number=2,Type=Float,Description=\"Median "
    "read-length-normalised alignment scores: [0] of reads "
    "supporting the "
    "variant; [1] of all reads covering the variant site.\">";
constexpr std::string_view RCMPLX_HEADER =
    "##INFO=<ID=RCMPLX,Number=1,Type=Float,Description=\"Mean "
    "100-base "
    "window complexity (Lempel-Ziv 76 entropy rate) of the "
    "reference "
    "within 400 bases either side of the REF allele.\">";

// --- registry --- //
constexpr std::array<VariantStat, 4> VARIANT_STATS = {{
    {{"QRK", QRK_HEADER, 2}, &compute_qrk},
    {{"TJAC", TJAC_HEADER, 2}, &compute_template_jaccard},
    {{"MLAS", MLAS_HEADER, 2}, &compute_mlas},
    {{"RCMPLX", RCMPLX_HEADER, 1}, &compute_rcmplx},
}};

std::span<const VariantStat> expos_field_registry()
{
  return VARIANT_STATS;
}
