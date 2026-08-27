#include "encode_info_field.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "shared/err.hpp"

constexpr std::string_view EXPOS_SKIP_HEADER =
    "##INFO=<ID=EXPOS_SKIP,Number=.,Type=String,Description="
    "\"Why expos "
    "produced no value, as '<scope>:<reason>' tokens. Scope is "
    "'record' "
    "(whole record skipped) or a statistic ID.\">";

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
    const VariantStatField& field, StatSkipReason reason
)
{
  return std::vector<StatValue> (
      static_cast<uint8_t> (field.nValues),
      StatValue{std::nullopt, reason}
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
  for (uint64_t i = 0; i < values.size(); ++i) {
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

std::vector<Skip> stat_skips (
    const VariantStatField& field,
    const std::vector<StatValue>& values
)
{
  std::vector<Skip> skips;
  for (const auto& v : values) {
    if (!v.reason) {
      continue;
    }
    const Skip skip = make_stat_skip (field.id, *v.reason);
    if (std::find (skips.begin(), skips.end(), skip) ==
        skips.end()) {
      skips.push_back (skip);
    }
  }
  return skips;
}

VoidOrErr set_expos_skip (
    bcf_hdr_t* hdr, bcf1_t* rec, std::span<const Skip> skips
)
{
  if (skips.empty()) {
    return {};
  }
  std::string joined;
  for (uint64_t i = 0; i < skips.size(); ++i) {
    if (i != 0) {
      joined += ',';
    }
    joined += to_info_format (skips[i]);
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
