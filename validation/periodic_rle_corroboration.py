#!/usr/bin/env python3
"""Do low RCMPLX windows contain bounded-period tandem repeats?

Usage: periodic_rle_corroboration.py ANNOTATED.vcf ref.fa [--flank 250]
"""
import argparse
import re
import statistics
import subprocess
import sys


def info_float(info, key, idx=0):
    m = re.search(r"(?:^|;)" + key + r"=([^;]+)", info)
    if not m:
        return None
    try:
        return float(m.group(1).split(",")[idx])
    except ValueError:
        return None


def periodic_rle(s, max_k):
    """Greedily segment s into maximal exact tandem-repeat runs, period 1..max_k."""
    seg_lens = []
    n = len(s)
    i = 0
    while i < n:
        best_run_bp = 1
        k_limit = min(max_k, n - i)
        for k in range(1, k_limit + 1):
            unit = s[i:i + k]
            run_period = 1
            while i + (run_period + 1) * k <= n:
                cand = s[i + run_period * k: i + (run_period + 1) * k]
                if cand == unit:
                    run_period += 1
                else:
                    break
            if run_period >= 2:
                run_bp = run_period * k
                if run_bp > best_run_bp:
                    best_run_bp = run_bp
        seg_lens.append(best_run_bp)
        i += best_run_bp
    return seg_lens


def top10_stats(seg_lens):
    n = len(seg_lens)
    if n == 0:
        return 0, 0.0
    top = sorted(seg_lens)[int(n * 0.9):]
    return max(top), sum(top) / len(top)


def spearman(xs, ys):
    def ranks(vs):
        order = sorted(range(len(vs)), key=lambda i: vs[i])
        r = [0.0] * len(vs)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and vs[order[j + 1]] == vs[order[i]]:
                j += 1
            avg_rank = (i + j) / 2.0
            for m in range(i, j + 1):
                r[order[m]] = avg_rank
            i = j + 1
        return r

    rx, ry = ranks(xs), ranks(ys)
    n = len(xs)
    mx, my = sum(rx) / n, sum(ry) / n
    cov = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    vx = sum((a - mx) ** 2 for a in rx)
    vy = sum((b - my) ** 2 for b in ry)
    if vx == 0 or vy == 0:
        return float("nan")
    return cov / (vx * vy) ** 0.5


def fetch_batch(ref_path, regions, batch_size=500):
    seqs = {}
    for i in range(0, len(regions), batch_size):
        chunk = regions[i:i + batch_size]
        proc = subprocess.run(
            ["samtools", "faidx", ref_path, *chunk],
            capture_output=True, text=True, check=True,
        )
        header, seq = None, []
        for line in proc.stdout.splitlines():
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(seq)
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if header is not None:
            seqs[header] = "".join(seq)
    return seqs


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("vcf")
    ap.add_argument("ref")
    ap.add_argument("--flank", type=int, default=250,
                     help="must match the --flank expos was run with")
    ap.add_argument("--max-k", type=int, default=5)
    ap.add_argument("--bins", type=int, default=10)
    a = ap.parse_args()

    records = []
    for line in open(a.vcf):
        if line.startswith("#"):
            continue
        f = line.split("\t")
        rcmplx = info_float(f[7], "RCMPLX", 0)
        if rcmplx is None:
            continue
        chrom, pos = f[0], int(f[1])
        start = max(1, pos - a.flank)
        end = pos - 1 + a.flank
        region = f"{chrom}:{start}-{end}"
        records.append(dict(rcmplx=rcmplx, region=region))

    print(f"RCMPLX-computed records: {len(records):,}   fetching reference...")
    seqs = fetch_batch(a.ref, [r["region"] for r in records])

    rows = []
    for r in records:
        seq = seqs.get(r["region"])
        if not seq or "N" in seq.upper():
            continue
        seg_lens = periodic_rle(seq.upper(), a.max_k)
        run_max, run_top10_mean = top10_stats(seg_lens)
        rows.append(dict(**r, run_max=run_max, run_top10_mean=run_top10_mean))

    N = len(rows)
    print(f"periodic_rle computed: {N:,}\n")

    rc10 = spearman([r["rcmplx"] for r in rows], [r["run_top10_mean"] for r in rows])
    rc = spearman([r["rcmplx"] for r in rows], [r["run_max"] for r in rows])
    print(f"Spearman(RCMPLX, periodic top10-mean run) = {rc10:+.3f}")
    print(f"Spearman(RCMPLX, periodic run_max)        = {rc:+.3f}  (noisy -- one repeat")
    print("                                              anywhere in the flank counts)")
    print("\nRCMPLX decile breakdown (n, median top10-mean run, median run_max):")
    srt = sorted(rows, key=lambda r: r["rcmplx"])
    qsize = len(srt) // a.bins
    for q in range(a.bins):
        chunk = srt[q * qsize: (q + 1) * qsize] if q < a.bins - 1 else srt[q * qsize:]
        lo, hi = chunk[0]["rcmplx"], chunk[-1]["rcmplx"]
        rtop10 = statistics.median(c["run_top10_mean"] for c in chunk)
        rmax = statistics.median(c["run_max"] for c in chunk)
        print(f"  Q{q + 1:<2} [{lo:.3f}-{hi:.3f}] n={len(chunk):4d}  "
              f"run_top10_mean={rtop10:6.2f}  run_max={rmax:6.1f}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
