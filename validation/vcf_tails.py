#!/usr/bin/env python3
"""Behaviour of QRK and TJAC on a real annotated VCF.

Usage:
    vcf_tails.py ANNOTATED.vcf [--region CHR:START-END] [--stat QRK]
"""
import argparse
import re
import statistics
import sys
from collections import Counter

BINS = [(-1e9, 0), (0, 1.95), (1.95, 3.0), (3.0, 4.0), (4.0, 8.0), (8.0, 1e9)]
TAILS = [(1.95, 5.0), (2.7, None), (3.0, 1.0)]  # cutoff, expected null %


def info_float(info, key, idx=0):
    m = re.search(r"(?:^|;)" + key + r"=([^;]+)", info)
    if not m:
        return None
    try:
        return float(m.group(1).split(",")[idx])
    except ValueError:
        return None  # missing subfield '.'


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("vcf")
    ap.add_argument("--stat", default="QRK", help="INFO ID to analyse")
    ap.add_argument("--region", help="difficult region, e.g. 21:9000000-12000000")
    a = ap.parse_args()

    chrom = start = end = None
    if a.region:
        chrom, span = a.region.split(":")
        start, end = (int(x) for x in span.split("-"))

    rows, skips, n = [], Counter(), 0
    for line in open(a.vcf):
        if line.startswith("#"):
            continue
        n += 1
        f = line.split("\t")
        s = re.search(r"(?:^|;)EXPOS_SKIP=([^;]+)", f[7])
        if s:
            skips.update(s.group(1).split(","))
        z = info_float(f[7], a.stat, 0)
        if z is None:
            continue
        rows.append((f[0], int(f[1]), z, info_float(f[7], a.stat, 1),
                     info_float(f[7], "RCMPLX", 0)))

    N = len(rows)
    print(f"records: {n:,}   {a.stat} computed: {N:,} ({100.0 * N / n:.1f}%)")
    for tok, c in skips.most_common(8):
        print(f"   {tok:42s} {c:7,}  ({100.0 * c / n:.1f}%)")
    if not N:
        return 1

    print(f"\nmedian z {statistics.median(r[2] for r in rows):.2f}; upper tail:")
    for cut, expect in TAILS:
        c = sum(1 for r in rows if r[2] > cut)
        exp = f"{expect:.1f}%" if expect else "-"
        print(f"   z > {cut:<5} {100.0 * c / N:5.2f}%  ({c:,})   null {exp}")

    def in_region(r):
        return chrom and r[0] == chrom and start <= r[1] < end

    print(f"\n{'z bin':>14} {'n':>7} {'in region':>10} {'RCMPLX':>7}")
    if chrom:
        base = 100.0 * sum(1 for r in rows if in_region(r)) / N
        print(f"{'(baseline)':>14} {N:>7} {base:>9.1f}%")
    for lo, hi in BINS:
        sub = [r for r in rows if lo <= r[2] < hi]
        if not sub:
            continue
        lbl = f"{lo:g} - {hi:g}".replace("-1e+09", "min").replace("1e+09", "max")
        reg = (f"{100.0 * sum(1 for r in sub if in_region(r)) / len(sub):>9.1f}%"
               if chrom else f"{'-':>10}")
        rcx = [r[4] for r in sub if r[4] is not None]
        print(f"{lbl:>14} {len(sub):>7} {reg} "
              f"{statistics.median(rcx) if rcx else float('nan'):>7.3f}")

    if chrom:
        print(f"\n{'subset':>22} {'n':>7} " + " ".join(f"z>{c:<5}" for c, _ in TAILS))
        for lbl, sel in [("all loci", lambda r: True),
                         ("excluding region", lambda r: not in_region(r)),
                         ("region only", in_region)]:
            sub = [r for r in rows if sel(r)]
            if not sub:
                continue
            pcts = " ".join(
                f"{100.0 * sum(1 for r in sub if r[2] > c) / len(sub):5.1f}%"
                for c, _ in TAILS)
            print(f"{lbl:>22} {len(sub):>7} {pcts}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
