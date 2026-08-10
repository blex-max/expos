#include "extract_pileup.hpp"

#include <fmt/format.h>
#include <htslib/sam.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <cstdint>
#include <expected>
#include <optional>
#include <string>
#include <unordered_set>

#include "hts/hts_types.hpp"
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
  auto* const mc{bam_aux2Z (raw_mc)};
  if (mc == NULL && errno == EINVAL) {
    return std::unexpected (make_err (
        fmt::format (
            "MC tag is not a string type for read {}. Try "
            "samtools fixmate?",
            qname
        )
    ));
  }
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
  auto* const rawAs = bam_aux_get (b1, "AS");
  if (rawAs == NULL) {
    return std::unexpected (
        make_err (fmt::format ("no AS tag for read {}", qname))
    );
  }
  const auto intAs = bam_aux2i (rawAs);
  if (intAs == 0 && errno == EINVAL) {
    return std::unexpected (make_err (
        fmt::format (
            "AS tag is not an integer type for read {}", qname
        )
    ));
  }
  return static_cast<double> (intAs) /
         static_cast<double> (b1->core.l_qseq);
}

VoidOrErr extract_features (PileupView plp, PileupFeatures& _out)
{
  std::unordered_set<std::string> qnamesSeen;
  auto* ru_bCigBuf = bam_init1();

  for (const auto& p1 : plp) {
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
        _out.endpoints.emplace_back (*epRet);
      }
    }

    if (!p1.is_del && !p1.is_refskip && p1.qpos >= 0) {
      _out.qPos.emplace_back (p1.qpos);
    }
    _out.readLen.emplace_back (b1->core.l_qseq);

    const auto nasRet = extract_normalised_as (b1);
    if (!nasRet) {
      bam_destroy1 (ru_bCigBuf);
      return std::unexpected (nasRet.error());
    }
    _out.normalisedAs.emplace_back (*nasRet);
  }

  bam_destroy1 (ru_bCigBuf);

  return {};
}

VoidOrErr extract_partition_features (
    PileupView plp, const VcfRec& rec,
    PileupFeatures& _outSupport, PileupFeatures& _outAll
)
{
  std::unordered_set<std::string> qnamesSeenSupporting;
  std::unordered_set<std::string> qnamesSeenAll;
  auto* ru_bCigBuf = bam_init1();
  std::optional<TemplateEndpoints> ru_endpoints;

  for (const auto& p1 : plp) {
    ru_endpoints.reset();
    const auto* b1 = p1.b;

    const auto mateMapped =
        (b1->core.tid == b1->core.mtid &&
         ((b1->core.flag & BAM_FMUNMAP) == 0));

    const std::string qname = bam_get_qname (b1);
    auto qnameNewSupporting = false;
    auto qnameNewAll = false;

    if (mateMapped) {
      qnameNewSupporting =
          !qnamesSeenSupporting.contains (qname);
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
      _outAll.qPos.emplace_back (p1.qpos);
    }
    if (qnameNewAll && ru_endpoints) {
      _outAll.endpoints.emplace_back (*ru_endpoints);
      qnamesSeenAll.insert (qname);
    }
    _outAll.readLen.emplace_back (b1->core.l_qseq);
    _outAll.normalisedAs.emplace_back (normalisedAs);

    if (rec.eval_read_support (rec, p1)) {
      if (qPosAln) {
        _outSupport.qPos.emplace_back (p1.qpos);
      }
      if (qnameNewSupporting && ru_endpoints) {
        _outSupport.endpoints.emplace_back (*ru_endpoints);
        qnamesSeenSupporting.insert (qname);
      }
      _outSupport.readLen.emplace_back (b1->core.l_qseq);
      _outSupport.normalisedAs.emplace_back (normalisedAs);
    }
  }

  bam_destroy1 (ru_bCigBuf);

  return {};
}
