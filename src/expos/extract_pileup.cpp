#include "extract_pileup.hpp"

#include <fmt/format.h>
#include <htslib/sam.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <expected>
#include <optional>
#include <string>
#include <unordered_set>

#include "vcf_record.hpp"

static std::expected<TemplateEndpoints, Err> extract_endpoints (
    const bam1_t* b1, bam1_t* bCigBuf
)
{
  std::array<int64_t, 4> endpointArr;

  const auto* qname = bam_get_qname (b1);

  auto* const raw_mc = bam_aux_get (b1, "MC");
  if (raw_mc == NULL) {
    return std::unexpected (make_err (
        fmt::format (
            "no MC tag for read {}. Try samtools fixmate?", qname
        )
    ));
  }
  if (bam_aux_type (raw_mc) != 'Z') {
    return std::unexpected (make_err (
        fmt::format (
            "MC tag is not of type 'Z' for read {}. "
            "Record data corrupt; type 'Z' is mandated "
            "for MC tag by SAM format spec. Try samtools"
            "fixmate?",
            qname
        )
    ));
  }
  auto* const mc{bam_aux2Z (raw_mc)};
  if (bam_parse_cigar (mc, NULL, bCigBuf) < 1) {
    return std::unexpected (make_err ((fmt::format (
        "unable to parse MC tag {} as cigar string "
        "for read {}. Try samtools fixmate?",
        mc, qname
    ))));
  }

    // get template region
  endpointArr[0] = b1->core.pos;
  endpointArr[1] =
      endpointArr[0] +
      bam_cigar2rlen (
          static_cast<int> (b1->core.n_cigar), bam_get_cigar (b1)
      );
  endpointArr[2] = b1->core.mpos;  // leftmost mate coordinate
  endpointArr[3] = endpointArr[2] +
                   bam_cigar2rlen (
                       static_cast<int> (bCigBuf->core.n_cigar),
                       bam_get_cigar (bCigBuf)
                   );

  const auto endpoints = std::minmax_element (
      endpointArr.begin(),
      endpointArr.end()
  );  // NOTE: returns pair of ptrs

  return std::pair{*endpoints.first, *endpoints.second};
}

// Read-length-normalised alignment score (AS tag / l_qseq).
static std::expected<double, Err> extract_normalised_as (
    const bam1_t* b1
)
{
  const auto* qname = bam_get_qname (b1);
  auto* const raw_as = bam_aux_get (b1, "AS");
  if (raw_as == NULL) {
    return std::unexpected (
        make_err (fmt::format ("no AS tag for read {}", qname))
    );
  }
  const auto asType = bam_aux_type (raw_as);
  if (asType != 'i' && asType != 'C') {
    return std::unexpected (make_err (fmt::format (
        "AS tag is not an integer type for read {}; SAM spec "
        "mandates type 'i'",
        qname
    )));
  }
  return static_cast<double> (bam_aux2i (raw_as)) /
         static_cast<double> (b1->core.l_qseq);
}

VoidOrErr extract_features (
    const PreparedPileup& plp, PileupFeatures& out
)
{
  std::unordered_set<std::string> qnamesSeen;
  auto* ru_bCigBuf = bam_init1();

  for (size_t i = 0; i < plp.nPlp; ++i) {
    const auto& p1 = plp.plpArr[i];
    const auto* b1 = p1.b;

    const auto mateMapped =
        (b1->core.tid == b1->core.mtid &&
         ((b1->core.flag & BAM_FMUNMAP) == 0));

    if (mateMapped) {
      const std::string qname = bam_get_qname (b1);
      if (!qnamesSeen.contains (qname)) {
        qnamesSeen.insert (qname);
        const auto epRet = extract_endpoints (b1, ru_bCigBuf);
        if (!epRet) {
          bam_destroy1 (ru_bCigBuf);
          return std::unexpected (epRet.error());
        }
        out.endpoints.emplace_back (*epRet);
      }
    }

    if (!p1.is_del && !p1.is_refskip && p1.qpos >= 0) {
      out.qPos.emplace_back (p1.qpos);
    }
    out.readLen.emplace_back (b1->core.l_qseq);

    const auto nasRet = extract_normalised_as (b1);
    if (!nasRet) {
      bam_destroy1 (ru_bCigBuf);
      return std::unexpected (nasRet.error());
    }
    out.normalisedAs.emplace_back (*nasRet);
  }

  bam_destroy1 (ru_bCigBuf);

  return {};
}

VoidOrErr extract_partition_features (
    const PreparedPileup& plp, const VcfRec& rec,
    PileupFeatures& outSupport, PileupFeatures& outAll
)
{
  std::unordered_set<std::string> qnamesSeenSupporting;
  std::unordered_set<std::string> qnamesSeenAll;
  auto* ru_bCigBuf = bam_init1();
  std::optional<TemplateEndpoints> ru_endpoints;

  for (size_t i = 0; i < plp.nPlp; ++i) {
    ru_endpoints.reset();
    const auto& p1 = plp.plpArr[i];
    const auto* b1 = p1.b;

    const auto mateMapped =
        (b1->core.tid == b1->core.mtid &&
         ((b1->core.flag & BAM_FMUNMAP) == 0));

    const std::string qname = bam_get_qname (b1);
    auto qnameNewSupporting = false;
    auto qnameNewAll = false;

    if (mateMapped) {
      qnameNewSupporting = !qnamesSeenSupporting.contains (qname);
      qnameNewAll = !qnamesSeenAll.contains (qname);

      if (qnameNewSupporting || qnameNewAll) {
        const auto epRet = extract_endpoints (b1, ru_bCigBuf);
        if (!epRet) {
          bam_destroy1 (ru_bCigBuf);
          return std::unexpected (epRet.error());
        }
        ru_endpoints = *epRet;
      }
    }

    const auto qPosAln =
        !p1.is_del && !p1.is_refskip && p1.qpos >= 0;

    const auto nasRet = extract_normalised_as (b1);
    if (!nasRet) {
      bam_destroy1 (ru_bCigBuf);
      return std::unexpected (nasRet.error());
    }
    const double normalisedAs = *nasRet;

    if (qPosAln) {
      outAll.qPos.emplace_back (p1.qpos);
    }
    if (qnameNewAll && ru_endpoints) {
      outAll.endpoints.emplace_back (*ru_endpoints);
      qnamesSeenAll.insert (qname);
    }
    outAll.readLen.emplace_back (b1->core.l_qseq);
    outAll.normalisedAs.emplace_back (normalisedAs);

    if (rec.eval_read_support (rec, p1)) {
      if (qPosAln) {
        outSupport.qPos.emplace_back (p1.qpos);
      }
      if (qnameNewSupporting && ru_endpoints) {
        outSupport.endpoints.emplace_back (*ru_endpoints);
        qnamesSeenSupporting.insert (qname);
      }
      outSupport.readLen.emplace_back (b1->core.l_qseq);
      outSupport.normalisedAs.emplace_back (normalisedAs);
    }
  }

  bam_destroy1 (ru_bCigBuf);

  return {};
}
