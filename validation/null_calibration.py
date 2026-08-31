#!/usr/bin/env python3
"""Independently confirm the suggestions made in the documentation of a 3.0 cutoff
on effect size for QRK and TJAC.

Usage:
    null_calibration.py BAM REGION [REGION ...]
"""
import argparse
import subprocess
import sys
from collections import defaultdict

import numpy as np

RADIUS = 5  # must match QPOS_RADIUS in expos/compute_info_field.cpp
EXCLUDE = 3840  # BAM_FSECONDARY|FQCFAIL|FDUP|FSUPPLEMENTARY, per expos/pileup_fn.hpp
N_VALUES = [2, 3, 4, 5, 8, 15, 30, 60]
MIN_BACKGROUND = 10  # must match MIN_BACKGROUND in expos/guards.hpp
CUTOFFS = [1.95, 2.0, 2.33, 2.7, 3.0, 3.5]

CONSUME_REF = set("MDN=X")
CONSUME_QRY = set("MIS=X")


def rlen(cigar):
    total, num = 0, ""
    for ch in cigar:
        if ch.isdigit():
            num += ch
        else:
            if ch in CONSUME_REF:
                total += int(num)
            num = ""
    return total


def scan(bam, region):
    """(query-position columns, template-interval columns) for one region."""
    qcols, tcols = defaultdict(list), defaultdict(dict)
    proc = subprocess.run(
        ["samtools", "view", "-F", str(EXCLUDE), bam, region],
        capture_output=True, text=True, check=True,
    )
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        qname, flag, rname, pos, cigar = f[0], int(f[1]), f[2], int(f[3]) - 1, f[5]
        if cigar == "*":
            continue
        rp, qp, num = pos, 0, ""
        for ch in cigar:
            if ch.isdigit():
                num += ch
                continue
            ln, num = int(num), ""
            if ch in CONSUME_REF and ch in CONSUME_QRY:
                for k in range(ln):
                    qcols[rp + k].append(qp + k)
            if ch in CONSUME_REF:
                rp += ln
            if ch in CONSUME_QRY:
                qp += ln
        # template interval, one per qname, matching extract_pileup.cpp
        if flag & 0x8 or f[6] not in ("=", rname):
            continue
        mc = next((t[5:] for t in f[11:] if t.startswith("MC:Z:")), None)
        if mc is None:
            continue
        mpos = int(f[7]) - 1
        c = (pos, pos + rlen(cigar), mpos, mpos + rlen(mc))
        for k in range(rlen(cigar)):
            tcols[pos + k].setdefault(qname, (min(c), max(c)))
    return qcols, {k: list(v.values()) for k, v in tcols.items()}


def draw_idx(N, n, ndraw, rng):
    return np.argpartition(rng.random((ndraw, N)), n - 1, axis=1)[:, :n]


def qrk_stat(pop, n, ndraw, rng):
    """Unordered query-position pairs closer than RADIUS (strictly <)."""
    s = np.sort(pop[draw_idx(len(pop), n, ndraw, rng)], axis=1)
    tot = np.zeros(ndraw, dtype=np.int64)
    for k in range(1, n):
        tot += (s[:, k:] - s[:, :-k] < RADIUS).sum(axis=1)
    return tot.astype(np.float64)


def tjac_stat(pop, n, ndraw, rng):
    """Sum of pairwise Jaccard overlap over template intervals."""
    idx = draw_idx(len(pop), n, ndraw, rng)
    a, b = pop[idx, 0].astype(np.float64), pop[idx, 1].astype(np.float64)
    ln = b - a
    tot = np.zeros(ndraw)
    for i in range(n):
        for j in range(i + 1, n):
            ov = np.minimum(b[:, i], b[:, j]) - np.maximum(a[:, i], a[:, j])
            np.maximum(ov, 0.0, out=ov)
            tot += ov / (ln[:, i] + ln[:, j] - ov)
    return tot


def calibrate(name, colsets, statfn, ndraw, ncols, rng, is_count=False):
    # var/mean is only interpretable for a count; TJAC sums a continuous
    # per-pair quantity, so the ratio carries units and means nothing.
    disp = f"{'var/mean':>9} " if is_count else ""
    print(f"\n{'=' * 79}\n{name}\n{'=' * 79}")
    print(f"{'n':>3} {'cols':>5} {disp}{'z@95':>6} {'z@99':>6} "
          f"{'z@99 by col':>13}   FPR of (z > cutoff AND p < 0.05)"
          f"        z>3.0")
    print(f"{'':>3} {'':>5} {'':>9}"[:11 + len(disp)] + f"{'':>6} {'':>6} {'':>13}   "
          + " ".join(f"{c:>5}" for c in CUTOFFS) + "        alone")
    for n in N_VALUES:
        zs, ps, percol, vms = [], [], [], []
        for cols in colsets:
            keys = sorted(cols)
            step = max(1, len(keys) // ncols)
            for k in keys[::step][:ncols]:
                pop = np.asarray(cols[k])
                # size_guard: nBg >= 2*nObs and nBg >= MIN_BACKGROUND. The
                # floor only bites below n=5; above that 2*n subsumes it.
                if len(pop) < 2 * n or len(pop) < MIN_BACKGROUND:
                    continue
                c = statfn(pop, n, ndraw, rng)
                mu, sd = c.mean(), c.std()
                if sd == 0:  # zero_variance; the code reports no effect size
                    continue
                zs.append((c - mu) / sd)
                # Each draw in turn plays the observed value; the rest are its
                # null. p = (#{others >= obs} + 1)/ndraw = #{all >= obs}/ndraw.
                srt = np.sort(c)
                ps.append((len(c) - np.searchsorted(srt, c, "left")) / len(c))
                percol.append(np.percentile(zs[-1], 99))
                vms.append(sd * sd / mu if mu > 0 else np.nan)
        if not zs:
            print(f"{n:>3} {0:>5}   (no columns deep enough)")
            continue
        z, p = np.concatenate(zs), np.concatenate(ps)
        sig = p < 0.05
        fpr = [100.0 * ((z > c) & sig).mean() for c in CUTOFFS]
        # The conjunction is what the docs recommend, but users do read the
        # effect size on its own. At low support the two diverge sharply, so
        # the unconditional rate at the published cutoff is reported too.
        zonly = 100.0 * (z > 3.0).mean()
        print(f"{n:>3} {len(percol):>5} "
              + (f"{np.nanmedian(vms):>9.2f} " if is_count else "")
              + f"{np.percentile(z, 95):>6.2f} {np.percentile(z, 99):>6.2f} "
              f"{f'{min(percol):.2f}-{max(percol):.2f}':>13}   "
              + " ".join(f"{v:>4.1f}%" for v in fpr)
              + f"        {zonly:>4.1f}%")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("bam")
    ap.add_argument("regions", nargs="+")
    ap.add_argument("--draws-qrk", type=int, default=20000)
    ap.add_argument("--draws-tjac", type=int, default=4000)
    ap.add_argument("--cols-qrk", type=int, default=50)
    ap.add_argument("--cols-tjac", type=int, default=25)
    ap.add_argument("--seed", type=int, default=20260806)
    a = ap.parse_args()

    rng = np.random.default_rng(a.seed)
    qsets, tsets = [], []
    for region in a.regions:
        q, t = scan(a.bam, region)
        if not q:
            print(f"{region}: no usable reads", file=sys.stderr)
            continue
        qsets.append(q)
        tsets.append(t)
        print(f"{region:>26}: {len(q)} columns, read depth median "
              f"{int(np.median([len(v) for v in q.values()]))}, template depth "
              f"median {int(np.median([len(v) for v in t.values()]))}")
    if not qsets:
        return 1

    print("\nvar/mean is 1.0 for a Poisson count; a sum of Bernoulli(p) pair")
    print("indicators gives 1-p, and sampling a large fraction of a finite")
    print("column pulls it lower still. z@99 is the 1-in-100 point of the null")
    print("effect size -- 2.33 if it were normal.")
    calibrate("QRK -- query-position pair count, radius 5", qsets,
              qrk_stat, a.draws_qrk, a.cols_qrk, rng, is_count=True)
    calibrate("TJAC -- template pairwise Jaccard sum", tsets,
              tjac_stat, a.draws_tjac, a.cols_tjac, rng)
    return 0


if __name__ == "__main__":
    sys.exit(main())
