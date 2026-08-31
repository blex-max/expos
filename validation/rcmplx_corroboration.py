#!/usr/bin/env python3
"""Does a low RCMPLX co-occur with degraded alignment score at the locus?

Reads the AS tag directly from the BAM for reads overlapping each variant.

Usage: rcmplx_corroboration.py ANNOTATED.vcf sample.bam [--cutoff 20] [--region CHR:START-END]
"""
import argparse
import os
import re
import statistics
import subprocess
import sys
import tempfile
from collections import defaultdict

EXCLUDE = 3840  # BAM_FSECONDARY|FQCFAIL|FDUP|FSUPPLEMENTARY, per expos/pileup_fn.hpp
CONSUME_REF = set("MDN=X")
CONSUME_QRY = set("MIS=X")
ALIGNED = set("M=X")  # ref and query consumed together: a real base call


def info_float(info, key, idx=0):
    m = re.search(r"(?:^|;)" + key + r"=([^;]+)", info)
    if not m:
        return None
    try:
        return float(m.group(1).split(",")[idx])
    except ValueError:
        return None


def scan_read(pos0, cigar, seq, targets, on_hit):
    """targets: {ref_pos0: variant idx} for this read's chromosome."""
    rp, qp, num = pos0, 0, ""
    for ch in cigar:
        if ch.isdigit():
            num += ch
            continue
        ln, num = int(num), ""
        if ch in CONSUME_REF:
            aligned = ch in ALIGNED
            for k in range(ln):
                idx = targets.get(rp + k)
                if idx is not None:
                    on_hit(idx, seq[qp + k] if aligned else None)
            rp += ln
        if ch in CONSUME_QRY:
            qp += ln


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("vcf")
    ap.add_argument("bam")
    ap.add_argument("--cutoff", type=float, default=20.0)
    ap.add_argument("--bins", type=int, default=5, help="number of equal-count RCMPLX bins to report")
    ap.add_argument("--region", default="21:9000000-12000000",
                     help="externally-defined difficult region, e.g. acrocentric chr21")
    a = ap.parse_args()

    chrom, span = a.region.split(":")
    rstart, rend = (int(x) for x in span.split("-"))

    records, n = [], 0
    for line in open(a.vcf):
        if line.startswith("#"):
            continue
        n += 1
        f = line.split("\t")
        rcmplx = info_float(f[7], "RCMPLX")
        if rcmplx is None:
            continue
        ref, alt = f[3].strip(), f[4].split(",")[0].strip()
        if len(ref) != 1 or len(alt) != 1:
            continue  # SNVs only -- read support for other types isn't checked here
        records.append(dict(chrom=f[0], pos=int(f[1]), alt=alt.upper(), rcmplx=rcmplx))

    N = len(records)
    print(f"records: {n:,}   RCMPLX computed (SNVs): {N:,} ({100.0*N/n:.1f}%)")
    if not N:
        return 1

    by_chrom = defaultdict(dict)
    for i, r in enumerate(records):
        by_chrom[r["chrom"]][r["pos"] - 1] = i

    with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as bed:
        for r in sorted(records, key=lambda r: (r["chrom"], r["pos"])):
            bed.write(f"{r['chrom']}\t{r['pos']-1}\t{r['pos']}\n")
        bed_path = bed.name

    support_as, all_as = defaultdict(list), defaultdict(list)
    try:
        proc = subprocess.run(
            ["samtools", "view", "-F", str(EXCLUDE), "-M", "-L", bed_path, a.bam],
            capture_output=True, text=True, check=True,
        )
    finally:
        os.unlink(bed_path)

    for line in proc.stdout.splitlines():
        f = line.split("\t")
        rname, pos1, cigar, seq = f[2], int(f[3]), f[5], f[9]
        if cigar == "*" or seq == "*":
            continue
        targets = by_chrom.get(rname)
        if not targets:
            continue
        as_tag = next((t for t in f[11:] if t.startswith("AS:i:")), None)
        if as_tag is None:
            continue
        normalised_as = int(as_tag[5:]) / len(seq)

        def on_hit(idx, base):
            all_as[idx].append(normalised_as)
            if base is not None and base.upper() == records[idx]["alt"]:
                support_as[idx].append(normalised_as)

        scan_read(pos1 - 1, cigar, seq, targets, on_hit)

    rows = [
        dict(
            **r,
            as_support=statistics.median(support_as[i]) if support_as[i] else None,
            as_all=statistics.median(all_as[i]) if all_as[i] else None,
        )
        for i, r in enumerate(records)
    ]

    flagged = [r for r in rows if r["rcmplx"] < a.cutoff]
    unflagged = [r for r in rows if r["rcmplx"] >= a.cutoff]
    print(f"RCMPLX < {a.cutoff}: {len(flagged):,} ({100.0*len(flagged)/N:.1f}%)  "
          f"vs >= {a.cutoff}: {len(unflagged):,}")

    def in_region(r):
        return r["chrom"] == chrom and rstart <= r["pos"] < rend

    def summarize(field, label):
        fv = [r[field] for r in flagged if r[field] is not None]
        uv = [r[field] for r in unflagged if r[field] is not None]
        if not fv or not uv:
            print(f"  {label:17s}: insufficient data")
            return
        print(f"  {label:17s}: flagged median={statistics.median(fv):.4f} (n={len(fv)})   "
              f"unflagged median={statistics.median(uv):.4f} (n={len(uv)})")

    print("\nAlignment score at the locus, read directly from the BAM's AS tag")
    print("(higher is better -- a real effect means flagged is worse on both):")
    summarize("as_support", "AS[supporting]")
    summarize("as_all", "AS[all covering]")

    freg = 100.0 * sum(1 for r in flagged if in_region(r)) / len(flagged) if flagged else float("nan")
    ureg = 100.0 * sum(1 for r in unflagged if in_region(r)) / len(unflagged) if unflagged else float("nan")
    breg = 100.0 * sum(1 for r in rows if in_region(r)) / N
    print(f"\nRegion membership ({a.region}):")
    print(f"  baseline (all)   {breg:5.1f}%")
    print(f"  flagged          {freg:5.1f}%")
    print(f"  unflagged        {ureg:5.1f}%")

    print(f"\nRCMPLX decile breakdown (n, region%, median AS[supporting]/AS[all]):")
    srt = sorted(rows, key=lambda r: r["rcmplx"])
    NBINS = a.bins
    qsize = len(srt) // NBINS
    for q in range(NBINS):
        chunk = srt[q * qsize: (q + 1) * qsize] if q < NBINS - 1 else srt[q * qsize:]
        lo, hi = chunk[0]["rcmplx"], chunk[-1]["rcmplx"]
        reg = 100.0 * sum(1 for r in chunk if in_region(r)) / len(chunk)
        sup = [r["as_support"] for r in chunk if r["as_support"] is not None]
        allv = [r["as_all"] for r in chunk if r["as_all"] is not None]
        print(f"  Q{q+1} [{lo:.3f}-{hi:.3f}] n={len(chunk):4d} region={reg:5.1f}%  "
              f"AS[support]={statistics.median(sup) if sup else float('nan'):.3f} "
              f"AS[all]={statistics.median(allv) if allv else float('nan'):.3f}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
