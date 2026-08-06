#include "compute_field.hpp"

#include "expos/variant_stats.hpp"
#include "shared/stats.hpp"

// --- helpers --- //

static ValuesOrSkip mc_to_fields (const McResult& mc)
{
  if (!mc.effectSize) {
    return std::unexpected (REASON_ZERO_VARIANCE);
  }
  return std::vector<StatValue>{
      stat_value (*mc.effectSize), stat_value (mc.pValue)
  };
}

// A statistic needs >= 2 supporting observations and at least twice as
// many background reads (very, very liberal). Returns the skip reason if a guard fails.
static std::optional<std::string_view> size_guard (
    std::size_t nObs, std::size_t nBackground
)
{
  if (nObs < 2) {
    return REASON_INSUFFICIENT_SUPPORT;
  }
  if (nBackground < (2 * nObs)) {
    return REASON_INSUFFICIENT_BACKGROUND;
  }
  return std::nullopt;
}

// --- statistics --- //

static ValuesOrSkip compute_qrk (const VariantStatInputs& in)
{
  // Query position is an offset within the read, so a heterogeneous read
  // population would confound
  if (!in.readLenHomogeneous) {
    return std::unexpected (REASON_HETEROGENEOUS_READ_LENGTH);
  }

  constexpr uint64_t k_qposRadius = 5;
  const auto& obs = in.supporting.qPos;
  const auto& bg = in.all.qPos;
  if (const auto reason = size_guard (obs.size(), bg.size())) {
    return std::unexpected (*reason);
  }

  const double observedStat = static_cast<double> (
      count_pairs_within_1d (obs, k_qposRadius)
  );
  auto draw = [&] {
    return subsample_wo_replace (bg, obs.size(), in.rng);
  };
  auto stat = [&] (decltype (bg) sample) {
    return static_cast<double> (
        count_pairs_within_1d (sample, k_qposRadius)
    );
  };
  const auto mc = run_monte_carlo (observedStat, draw, stat);
  return mc_to_fields (mc);
}

static ValuesOrSkip compute_template_jaccard (
    const VariantStatInputs& in
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

  const auto drawFn = [&]() {
    return subsample_wo_replace (bg, obs.size(), in.rng);
  };

  const auto mc = run_monte_carlo (
      obsJaccardSum, drawFn, pairwise_jaccard_sum
  );

  return mc_to_fields (mc);
}

// The only statistic with independently-missing subfields: two unrelated
// medians, either of which can be absent while the other is reportable.
static ValuesOrSkip compute_mlas (const VariantStatInputs& in)
{
  return std::vector<StatValue>{
      stat_or (
          percentile (
              in.supporting.normalisedAs, PERCENTILE_MEDIAN
          ),
          REASON_NO_SUPPORT
      ),
      stat_or (
          percentile (in.all.normalisedAs, PERCENTILE_MEDIAN),
          REASON_NO_BACKGROUND
      )
  };
}

static ValuesOrSkip compute_rcmplx (const VariantStatInputs& in)
{
  constexpr std::size_t k_windowSize = 100;
  const std::string_view ref = in.refSlice;
  // too-short spans have no full window; masked/ambiguous bases make the
  // complexity meaningless.
  if (ref.size() < k_windowSize) {
    return std::unexpected (REASON_REFERENCE_TOO_SHORT);
  }
  if (ref.find_first_of ("Nn") != std::string_view::npos) {
    return std::unexpected (REASON_REFERENCE_HAS_N);
  }
  double entropySum = 0.0;
  std::size_t nWin = 0;
  for (; nWin + k_windowSize <= ref.size();
       ++nWin) {  // step of 1
    entropySum += entropy_lz76 (ref.substr (nWin, k_windowSize));
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
    "region spanned by supporting templates.\">";
constexpr std::string_view EXPOS_SKIP_HEADER =
    "##INFO=<ID=EXPOS_SKIP,Number=.,Type=String,Description="
    "\"Why expos "
    "produced no value, as '<scope>:<reason>' tokens. Scope is "
    "'record' "
    "(whole record skipped) or a statistic ID.\">";

// --- registry --- //
constexpr std::array<VariantStat, 4> VARIANT_STATS = {{
    {{"QRK", QRK_HEADER, 2}, &compute_qrk},
    {{"TJAC", TJAC_HEADER, 2}, &compute_template_jaccard},
    {{"MLAS", MLAS_HEADER, 2}, &compute_mlas},
    {{"RCMPLX", RCMPLX_HEADER, 1}, &compute_rcmplx},
}};

std::span<const VariantStat> variant_stats()
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

std::vector<StatValue> stat_all_missing (
    const VariantStatField& field, std::string_view reason
)
{
  return std::vector<StatValue> (
      static_cast<std::size_t> (field.nValues),
      stat_missing (reason)
  );
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
    return std::unexpected (make_err (
        "failed to write INFO field " + std::string (field.id)
    ));
  }
  return {};
}

std::vector<std::string> stat_skip_tokens (
    const VariantStatField& field,
    const std::vector<StatValue>& values
)
{
  std::vector<std::string> tokens;
  for (const auto& v : values) {
    if (!v.reason) {
      continue;
    }
    std::string token =
        std::string (field.id) + ":" + std::string (*v.reason);
    if (std::find (tokens.begin(), tokens.end(), token) ==
        tokens.end()) {
      tokens.push_back (std::move (token));
    }
  }
  return tokens;
}

VoidOrErr set_expos_skip (
    bcf_hdr_t* hdr, bcf1_t* rec,
    const std::vector<std::string>& tokens
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
  if (bcf_update_info_string (
          hdr, rec, "EXPOS_SKIP", joined.c_str()
      ) != 0) {
    return std::unexpected (
        make_err ("failed to write EXPOS_SKIP INFO field")
    );
  }
  return {};
}
