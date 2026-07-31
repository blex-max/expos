"""
Generate synthetic variant-call pileups and query-position visualisations
for the expos paper walkthrough.

Run with:  python generate.py --output-dir <dir>

To add a new case:
  1. add a gen_*(sam_path: Path) func.
  2. Append Case("<name>", gen_*) to CASES.
"""

import argparse
import random
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pysam
import seaborn as sb
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import htsgen

# pyright: reportAny = false
# pyright: reportUnusedCallResult = false
# pyright: reportUnknownMemberType = false

# ---------------------------------------------------------------------------
# Shared constants
# ---------------------------------------------------------------------------

READ_LEN = 50
PILEUP_POS = READ_LEN - 1      # genomic position of interest (0-based)
REF_SEQ = "T" * (READ_LEN * 2)
TID = 0
NREADS_ALT = 20
NREADS_REF = 60

coords = htsgen.PileupCoordinates(
    gstart = 0,
    gend = len(REF_SEQ),
    gpos = PILEUP_POS,
    tid = TID,
)
ppars = htsgen.PileupParams(
    coordinates = coords,
    refseq = REF_SEQ,
    readlen = READ_LEN,
)


# ---------------------------------------------------------------------------
# Visualisation
# ---------------------------------------------------------------------------

@dataclass
class _ReadSet:
    base: str
    qpos: list[int] = field(default_factory=list)


def visualise(
    aln_path: Path,
    ref_path: Path,
    output_path: Path,
    read_length: int = READ_LEN,
) -> None:
    aln = pysam.AlignmentFile(str(aln_path), "r")
    fasta = pysam.FastaFile(str(ref_path))

    ref_base = fasta.fetch(aln.references[0], PILEUP_POS, PILEUP_POS + 1).upper()
    qpos_range = (0, read_length - 1)

    set_ref: _ReadSet = _ReadSet(ref_base)
    sets_alt: dict[str, _ReadSet] = {}

    for col in aln.pileup():
        if col.reference_pos != PILEUP_POS:
            continue
        for read in col.pileups:
            qseq = read.alignment.query_sequence
            if qseq is None:
                continue
            if read.is_del:
                key = "DEL"
                if key not in sets_alt:
                    sets_alt[key] = _ReadSet(key)
                sets_alt[key].qpos.append(read.query_position_or_next)
                continue
            qpos = read.query_position
            if qpos is None:
                continue
            qbase = qseq[qpos].upper()
            if qbase == ref_base:
                set_ref.qpos.append(qpos)
            else:
                if qbase not in sets_alt:
                    sets_alt[qbase] = _ReadSet(qbase)
                sets_alt[qbase].qpos.append(qpos)

    n_plots = 1 + len(sets_alt)
    fig, axes = plt.subplots(n_plots, 1, figsize=(10, 2 * n_plots), sharex=True)
    if n_plots == 1:
        axes = [axes]

    axes[0].set_title(f"Ref ({set_ref.base})")
    sb.kdeplot(x=set_ref.qpos, fill=True, bw_adjust=0.5, ax=axes[0], clip=qpos_range)
    sb.rugplot(x=set_ref.qpos, height=0.15, color="black", ax=axes[0])
    axes[0].set_yticks([])
    axes[0].set_xlim(qpos_range)
    axes[0].text(0.97, 0.95, f"n = {len(set_ref.qpos)}", transform=axes[0].transAxes,
                 ha="right", va="top", fontsize=12)

    for i, alt_set in enumerate(sets_alt.values()):
        ax = axes[i + 1]
        ax.set_title(f"Alt ({alt_set.base})")
        sb.kdeplot(x=alt_set.qpos, fill=True, bw_adjust=0.5, ax=ax, clip=qpos_range)
        sb.rugplot(x=alt_set.qpos, height=0.15, color="black", ax=ax)
        ax.set_yticks([])
        ax.set_xlim(qpos_range)
        ax.text(0.97, 0.95, f"n = {len(alt_set.qpos)}", transform=ax.transAxes,
                ha="right", va="top", fontsize=12)

    axes[-1].set_xlabel("Query position")
    plt.tight_layout()
    plt.savefig(str(output_path))
    plt.close(fig)


# ---------------------------------------------------------------------------
# Generation
# ---------------------------------------------------------------------------

def gen_clustered_variant(sam_path: Path) -> None:
    """Alt qpos clustered around read midpoint; ref qpos broadly distributed."""

    midpoint = (READ_LEN // 2) - 1
    wobble = int(READ_LEN * 0.05)
    set_a = htsgen.PileupReadSet(
        event = htsgen.EventSpec(htsgen.BaseEvents.A),
        qpos_cb = lambda: random.randint(midpoint - wobble, midpoint + wobble),
    )
    set_ref = htsgen.PileupReadSet(
        event = htsgen.EventSpec(htsgen.BaseEvents.ref),
        qpos_cb = lambda: random.randint(0, READ_LEN - 1),
    )

    pileup = htsgen.generate_pileup(ppars, [(NREADS_ALT, set_a), (NREADS_REF, set_ref)])
    htsgen.write_pileup(pileup, str(sam_path))


def gen_broad_variant(sam_path: Path) -> None:
    """Both alt and ref qpos broadly distributed across query positions."""

    broad = lambda: random.randint(0, READ_LEN - 1)
    set_a = htsgen.PileupReadSet(
        event = htsgen.EventSpec(htsgen.BaseEvents.A),
        qpos_cb = broad,
    )
    set_ref = htsgen.PileupReadSet(
        event = htsgen.EventSpec(htsgen.BaseEvents.ref),
        qpos_cb = broad,
    )

    pileup = htsgen.generate_pileup(ppars, [(NREADS_ALT, set_a), (NREADS_REF, set_ref)])
    htsgen.write_pileup(pileup, str(sam_path))


def gen_mixed_clustering(sam_path: Path) -> None:
    """Both ref and alt clustered, but at different reference loci, and to differing degrees"""

    variant_center  = READ_LEN // 4
    set_a = htsgen.PileupReadSet(
        event = htsgen.EventSpec(htsgen.BaseEvents.A),
        qpos_cb = lambda: random.randint(variant_center - 2, variant_center + 2),
    )

    ref_center = (READ_LEN * 3) // 4
    nref_clustered = (NREADS_REF * 3) // 4
    nref_broad = NREADS_REF - nref_clustered
    set_ref_clustered = htsgen.PileupReadSet(
        event = htsgen.EventSpec(htsgen.BaseEvents.ref),
        qpos_cb = lambda: random.randint(ref_center - 10, ref_center + 10),
    )
    set_ref_broad = htsgen.PileupReadSet(
        event = htsgen.EventSpec(htsgen.BaseEvents.ref),
        qpos_cb = lambda: random.randint(0, READ_LEN - 1),
    )

    pileup = htsgen.generate_pileup(
        ppars,
        [(NREADS_ALT, set_a), (nref_clustered, set_ref_clustered), (nref_broad, set_ref_broad)],
    )
    htsgen.write_pileup(pileup, str(sam_path))


@dataclass
class Case:
    name: str  # used as the output filename stem
    generate: Callable[[Path], None]  # called with the SAM output path
CASES: list[Case] = [
    Case("clustered-variant", gen_clustered_variant),
    Case("broad-variant",     gen_broad_variant),
    Case("mixed-clustering",  gen_mixed_clustering),
]




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate synthetic variant-call pileups and visualisations.")
    parser.add_argument("outdir", type=Path, metavar="OUTDIR", help="Directory to write all output files into.")
    args = parser.parse_args()

    out: Path = args.outdir
    try:
        out.mkdir(exist_ok=True)
    except Exception as ex:
        print ("Could not use OUTDIR, reporting: {}".format(ex), file=sys.stderr)
        sys.exit(1)

    ref_path = out / "ref.fa"
    ref_path.write_text(f">chr1\n{REF_SEQ}\n")

    for case in CASES:
        sam_path = out / f"{case.name}.pileup.sam"
        tmp_path = out / f"{case.name}.pileup.tmp.sam"
        case.generate(sam_path)
        pysam.sort("-o", str(tmp_path), str(sam_path))
        sam_path.unlink()
        tmp_path.rename(sam_path)
        visualise(sam_path, ref_path, out / f"{case.name}.svg")

    sys.exit(0)
