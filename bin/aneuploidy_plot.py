#!/usr/bin/env python3
"""
Mosaic aneuploidy profiling from read depth and B-allele frequency.

In Mosaic Variegated Aneuploidy the defining cellular phenotype is that a
*fraction* of cells carry whole-chromosome gains or losses, and the chromosome
involved differs between cells. Bulk WGS therefore does not show clean copy
number steps; it shows small, chromosome-wide shifts in mean depth whose size
is proportional to the fraction of cells affected.

This script measures those shifts directly from the patient's own reads:

  1. GC-correct per-bin depth (GC bias is the dominant technical confounder at
     the scale of a whole chromosome).
  2. Normalise to the autosomal median -> per-bin copy-number ratio.
  3. Per chromosome, take the median ratio with a bootstrap CI over bins.
  4. Convert the shift to an estimated mosaic fraction:
         gain  (one extra copy in fraction f):  ratio = 1 + f/2  ->  f = 2(ratio - 1)
         loss  (one lost copy in fraction f):   ratio = 1 - f/2  ->  f = 2(1 - ratio)
  5. Cross-check with B-allele frequency at heterozygous sites, which shifts
     away from 0.5 in the same direction and is independent of depth bias.

Honesty about the detection limit is part of the output: the reported CI comes
from a bootstrap over bins, which assumes bins are independent. They are not
(GC waves and mappability are spatially correlated), so the true limit is
higher than the nominal CI. The table reports both the nominal CI and an
inflated limit using a block bootstrap.

Usage:
  aneuploidy_plot.py --bins SAMPLE.regions.bed.gz \\
                     --fasta /opt/refs/GRCh38.fa \\
                     --vcf  SAMPLE.dv.vcf.gz \\
                     --bin-size 100000 \\
                     --out-figure aneuploidy_karyotype.pdf \\
                     --out-table  mosaic_fractions.tsv
"""

import argparse
import gzip
import sys

import numpy as np

# Diverging palette: two poles + neutral gray midpoint.
# Validated colourblind-safe on a light surface (CVD dE 21.1 protan / 29.8 tritan).
C_LOSS = "#2166AC"
C_NEUTRAL = "#8C8C8C"
C_GAIN = "#B2182B"
# Flagged but not confirmed: a muted bronze, deliberately outside the red/blue
# result axis so a confounded shift never reads as a call. The QC gate is
# sample-wide, so bronze and red never share one figure.
C_UNCONF = "#A07A45"
C_INK = "#1a1a1a"
C_MUTED = "#6b6b6b"
C_GRID = "#e0e0e0"

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
CHROMS = AUTOSOMES + ["chrX"]


# ---------------------------------------------------------------- input

def read_bins(path):
    """mosdepth --by N regions.bed.gz -> chrom, start, end, mean_depth."""
    chrom, start, end, depth = [], [], [], []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] not in CHROMS:
                continue
            chrom.append(f[0])
            start.append(int(f[1]))
            end.append(int(f[2]))
            depth.append(float(f[3]))
    return (np.array(chrom), np.array(start), np.array(end), np.array(depth))


def gc_fraction(fasta_path, chrom, start, end):
    """GC fraction and N fraction per bin, read from the reference."""
    import pysam
    fa = pysam.FastaFile(fasta_path)
    gc = np.full(len(chrom), np.nan)
    nf = np.full(len(chrom), np.nan)
    for i in range(len(chrom)):
        try:
            seq = fa.fetch(chrom[i], int(start[i]), int(end[i])).upper()
        except (KeyError, ValueError):
            continue
        if not seq:
            continue
        n = seq.count("N")
        usable = len(seq) - n
        nf[i] = n / len(seq)
        if usable > 0:
            gc[i] = (seq.count("G") + seq.count("C")) / usable
    fa.close()
    return gc, nf


def read_baf(vcf_path, min_dp=15, max_dp=200):
    """B-allele frequency at heterozygous SNVs, from AD."""
    from cyvcf2 import VCF
    pos_by_chrom = {c: [] for c in CHROMS}
    baf_by_chrom = {c: [] for c in CHROMS}
    vcf = VCF(vcf_path)
    for v in vcf:
        if v.CHROM not in pos_by_chrom or not v.is_snp or len(v.ALT) != 1:
            continue
        gt = v.gt_types[0]
        if gt != 1:                      # cyvcf2: 1 == HET
            continue
        ad = v.format("AD")
        if ad is None:
            continue
        ref, alt = int(ad[0][0]), int(ad[0][1])
        dp = ref + alt
        if dp < min_dp or dp > max_dp:
            continue
        pos_by_chrom[v.CHROM].append(v.POS)
        baf_by_chrom[v.CHROM].append(alt / dp)
    vcf.close()
    return ({c: np.array(p) for c, p in pos_by_chrom.items()},
            {c: np.array(b) for c, b in baf_by_chrom.items()})


# ------------------------------------------------------------ correction

def gc_correct(depth, gc, keep, n_buckets=40):
    """Divide each bin by the median depth of bins with similar GC.

    GC buckets are quantiles, not fixed percentage points, so every bucket
    holds a usable number of bins regardless of how many bins exist.
    """
    corrected = depth.astype(float).copy()
    global_med = np.median(depth[keep])
    gk = gc[keep]
    n_buckets = max(4, min(n_buckets, keep.sum() // 25))
    edges = np.quantile(gk, np.linspace(0, 1, n_buckets + 1))
    edges[0] -= 1e-9
    edges = np.unique(edges)
    idx = np.digitize(gc, edges[1:-1], right=True)
    for b in np.unique(idx[keep]):
        sel = keep & (idx == b)
        if sel.sum() < 10:
            continue
        med = np.median(depth[sel])
        if med > 0:
            corrected[sel] = depth[sel] * (global_med / med)
    return corrected


def regional_gc(chrom, gc, half_window=10):
    """GC averaged over the surrounding +/- half_window bins (~2 Mb at 100 kb).

    Bin-level GC misses isochore-scale structure: gene-dense, early-replicating,
    GC-rich regions carry extra depth beyond what their own 100 kb GC predicts.
    Those regions are concentrated on chr16, 17, 19 and 22, so uncorrected they
    masquerade as whole-chromosome gains.
    """
    out = np.full(len(gc), np.nan)
    for c in np.unique(chrom):
        idx = np.where(chrom == c)[0]
        g = gc[idx]
        ok = np.isfinite(g)
        num = np.convolve(np.where(ok, g, 0.0), np.ones(2 * half_window + 1), "same")
        den = np.convolve(ok.astype(float), np.ones(2 * half_window + 1), "same")
        with np.errstate(invalid="ignore", divide="ignore"):
            out[idx] = np.where(den > 0, num / den, np.nan)
    return out


def _cell_index(x, train, n):
    edges = np.unique(np.quantile(x[train], np.linspace(0, 1, n + 1)))
    return np.digitize(x, edges[1:-1], right=True)


def gc_correct_2d_loo(depth, gc, rgc, chrom, keep, n_local=20, n_regional=8,
                      min_cell=25):
    """Two-scale GC correction, fitted leave-one-chromosome-out.

    Cells are (local-GC quantile x regional-GC quantile). For each chromosome
    the cell medians come only from OTHER chromosomes, so a genuine
    whole-chromosome gain cannot be absorbed into its own correction factor --
    it is judged against GC-matched regions elsewhere in the genome. An
    artefact shared with GC-rich regions on other chromosomes is removed; a
    real chromosome-specific shift survives.

    Cells with too few training bins fall back to the 1D (local GC) median.
    """
    li = _cell_index(gc, keep, n_local)
    ri = _cell_index(rgc, keep, n_regional)
    cell = li * (n_regional + 1) + ri
    corrected = depth.astype(float).copy()

    for c in np.unique(chrom[keep]):
        target = keep & (chrom == c)
        train = keep & (chrom != c)
        global_med = np.median(depth[train])
        for k in np.unique(cell[target]):
            sel_t = target & (cell == k)
            sel_tr = train & (cell == k)
            if sel_tr.sum() >= min_cell:
                med = np.median(depth[sel_tr])
            else:
                fb = train & (li == li[sel_t][0])
                med = np.median(depth[fb]) if fb.sum() >= min_cell else global_med
            if med > 0:
                corrected[sel_t] = depth[sel_t] * (global_med / med)
    return corrected


def gc_correct_poly_loo(depth, gc, rgc, chrom, keep, trim=0.02):
    """Smooth two-scale GC model, fitted leave-one-chromosome-out.

    log(depth) ~ quadratic(local GC) + quadratic(regional GC) + interaction,
    least squares on the other chromosomes' bins (tails trimmed for
    robustness). A parametric surface extrapolates sensibly to the GC extremes
    that only chr19 and chr22 reach, where a cell-based model runs out of
    training bins.
    """
    def design(g, r):
        g0, r0 = g - 0.41, r - 0.41
        return np.column_stack([np.ones_like(g0), g0, g0**2, r0, r0**2, g0 * r0])

    ld = np.log(np.where(depth > 0, depth, np.nan))
    corrected = depth.astype(float).copy()
    for c in np.unique(chrom[keep]):
        target = keep & (chrom == c)
        train = keep & (chrom != c) & np.isfinite(ld)
        lo, hi = np.quantile(ld[train], [trim, 1 - trim])
        train &= (ld >= lo) & (ld <= hi)
        X = design(gc[train], rgc[train])
        beta, *_ = np.linalg.lstsq(X, ld[train], rcond=None)
        pred = design(gc[target], rgc[target]) @ beta
        corrected[target] = depth[target] / np.exp(pred - beta[0])
    return corrected


def chrom_gc_rho(chrom, gc, ratio, keep):
    """Spearman correlation of per-chromosome mean GC vs median ratio.

    The direct test for the artefact: a real mosaic gain has no reason to
    favour GC-rich chromosomes. Values near 0 are what a clean correction
    should produce; the uncorrected profile of this sample gave ~0.6.
    """
    xs, ys = [], []
    for c in AUTOSOMES:
        s = keep & (chrom == c)
        if s.sum() < 50:
            continue
        xs.append(np.nanmean(gc[s]))
        ys.append(np.median(ratio[s]))
    rx = np.argsort(np.argsort(xs))
    ry = np.argsort(np.argsort(ys))
    return float(np.corrcoef(rx, ry)[0, 1])


def block_bootstrap_ci(values, n_boot=2000, block=20, rng=None):
    """Median CI with a moving-block bootstrap (bins are spatially correlated)."""
    rng = rng or np.random.default_rng(0)
    n = len(values)
    if n < block * 2:
        return (np.nan, np.nan)
    n_blocks = int(np.ceil(n / block))
    starts_pool = np.arange(0, n - block + 1)
    meds = np.empty(n_boot)
    for i in range(n_boot):
        st = rng.choice(starts_pool, size=n_blocks, replace=True)
        idx = (st[:, None] + np.arange(block)[None, :]).ravel()[:n]
        meds[i] = np.median(values[idx])
    return tuple(np.percentile(meds, [2.5, 97.5]))


# ------------------------------------------------------------- analysis

def analyse(chrom, ratio, baf_by_chrom, rng):
    """Per-chromosome median ratio, CI, and a call against an empirical floor.

    The baseline is the autosomal median, so roughly half of all chromosomes
    sit above 1.0 by construction and a CI excluding 1.0 is not on its own
    evidence of anything. A chromosome is called only when its shift exceeds
    BOTH its own block-bootstrap CI AND a noise floor derived from how much
    the other chromosomes' medians scatter (3 x MAD). That floor absorbs
    residual GC and mappability structure that survives correction.
    """
    # Leave-one-out baseline: a chromosome is compared against the other
    # autosomes, not against a genome median it helps define. Without this an
    # affected chromosome deflates its own apparent shift.
    stats = {}
    for c in CHROMS:
        sel = chrom == c
        if sel.sum() < 50:
            continue
        others = np.isin(chrom, [a for a in AUTOSOMES if a != c])
        base = float(np.median(ratio[others])) if others.sum() > 100 else 1.0
        vals = ratio[sel] / base
        stats[c] = {
            "n_bins": int(sel.sum()),
            "median": float(np.median(vals)),
            "ci": block_bootstrap_ci(vals, rng=rng),
        }

    auto_meds = np.array([s["median"] for c, s in stats.items() if c in AUTOSOMES])
    mad = float(np.median(np.abs(auto_meds - np.median(auto_meds))))
    floor = max(3 * 1.4826 * mad, 0.01)      # never claim better than 1% in ratio

    rows = []
    for c, s in stats.items():
        med = s["median"]
        lo, hi = s["ci"]
        dev = abs(med - 1.0)
        ci_excludes_1 = np.isfinite(lo) and (lo > 1.0 or hi < 1.0)
        sig = bool(ci_excludes_1 and dev > floor)

        call_dir = "gain" if med >= 1 else "loss"
        frac = 2 * dev                        # one copy in fraction f -> ratio 1 +/- f/2

        b = baf_by_chrom.get(c, np.array([]))
        baf_dev = float(np.median(np.abs(b - 0.5))) if len(b) > 100 else np.nan

        rows.append({
            "chrom": c,
            "n_bins": s["n_bins"],
            "median_ratio": med,
            "ci_lo": lo,
            "ci_hi": hi,
            "direction": call_dir if sig else "none",
            "mosaic_fraction": frac if sig else 0.0,
            "baf_median_abs_dev": baf_dev,
            "n_het_snv": int(len(b)),
            "significant": sig,
        })
    return rows, floor


# ------------------------------------------------------------- plotting

def plot(chrom, start, ratio, baf_pos, baf_val, rows, out_pdf, bin_size):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    call = {r["chrom"]: (r["direction"] if r.get("qc") == "pass"
                         else "unconfirmed" if r.get("qc") == "gc_confounded"
                         else "none")
            for r in rows}
    palette = {"gain": C_GAIN, "loss": C_LOSS, "none": C_NEUTRAL,
               "unconfirmed": C_UNCONF}

    # cumulative genome coordinates
    offsets, cum, ticks, labels = {}, 0, [], []
    for c in CHROMS:
        sel = chrom == c
        if not sel.any():
            continue
        offsets[c] = cum
        span = start[sel].max() + bin_size
        ticks.append(cum + span / 2)
        labels.append(c.replace("chr", ""))
        cum += span
    genome_len = cum

    x = np.array([offsets.get(c, 0) for c in chrom]) + start
    colours = np.array([
        palette[call.get(c, "none")]
        for c in chrom
    ])

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(11, 6.2), sharex=True,
        gridspec_kw={"height_ratios": [1.25, 1], "hspace": 0.12},
    )
    fig.patch.set_facecolor("white")

    # alternating chromosome bands for readability
    for i, c in enumerate(CHROMS):
        if c not in offsets:
            continue
        nxt = offsets[CHROMS[i + 1]] if i + 1 < len(CHROMS) and CHROMS[i + 1] in offsets else genome_len
        if i % 2 == 0:
            for ax in (ax1, ax2):
                ax.axvspan(offsets[c], nxt, color="#f5f5f4", lw=0, zorder=0)

    ax1.scatter(x, ratio, s=1.2, c=colours, alpha=0.45, linewidths=0, zorder=2)
    ax1.axhline(1.0, color=C_MUTED, lw=1, zorder=3)
    ax1.set_ylim(0.70, 1.30)
    ax1.set_ylabel("Copy-number ratio\n(GC-corrected depth)", fontsize=9, color=C_INK)
    ax1.set_title(
        "Mosaic aneuploidy profile from whole-genome sequencing depth and allele balance",
        fontsize=11, color=C_INK, pad=20, loc="left",
    )
    if any(r.get("qc") == "gc_confounded" for r in rows):
        ax1.text(0, 1.015, "Depth shifts track chromosome GC content; elevated "
                 "chromosomes are not called as mosaic gains (see methods).",
                 transform=ax1.transAxes, fontsize=7.5, color=C_MUTED, va="bottom")

    # direct labels on called chromosomes: identity is never colour-alone.
    # Small chromosomes sit close together, so labels are centred on their
    # chromosome and staggered vertically to avoid overprinting.
    called = [r for r in rows if r["significant"]]
    for k, r in enumerate(called):
        c = r["chrom"]
        sel = chrom == c
        centre = offsets[c] + (start[sel].max() + bin_size) / 2
        ax1.annotate(
            (f"{c.replace('chr','')} {r['direction']} ~{r['mosaic_fraction']*100:.0f}%"
             if r.get("qc") == "pass" else
             f"{c.replace('chr','')}: +{(r['median_ratio']-1)*100:.1f}% depth, unconfirmed"),
            xy=(centre, 1.27 - 0.035 * (k % 3)),
            ha="right" if centre > 0.8 * genome_len else "center",
            fontsize=7.5, color=C_INK, zorder=4,
        )

    for c in CHROMS:
        p, b = baf_pos.get(c, np.array([])), baf_val.get(c, np.array([]))
        if len(p) == 0:
            continue
        if len(p) > 40000:                       # thin for legibility
            idx = np.random.default_rng(0).choice(len(p), 40000, replace=False)
            p, b = p[idx], b[idx]
        ax2.scatter(p + offsets.get(c, 0), b, s=0.6,
                    c=palette[call.get(c, "none")],
                    alpha=0.18, linewidths=0, zorder=2)
    ax2.axhline(0.5, color=C_MUTED, lw=1, zorder=3)
    ax2.set_ylim(0, 1)
    ax2.set_ylabel("B-allele frequency\n(heterozygous SNVs)", fontsize=9, color=C_INK)
    ax2.set_xlabel("Chromosome", fontsize=9, color=C_INK)

    for ax in (ax1, ax2):
        ax.set_xlim(0, genome_len)
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels, fontsize=7.5, color=C_MUTED)
        ax.tick_params(axis="y", labelsize=7.5, colors=C_MUTED, length=3)
        ax.grid(axis="y", color=C_GRID, lw=0.6, zorder=1)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(C_GRID)

    from matplotlib.lines import Line2D
    ax1.legend(
        handles=[Line2D([], [], marker="o", ls="", ms=5, color=c, label=l)
                 for c, l in (((C_UNCONF, "flagged, unconfirmed"), (C_NEUTRAL, "no call"))
                              if "unconfirmed" in call.values() else
                              ((C_GAIN, "gain"), (C_NEUTRAL, "no call"), (C_LOSS, "loss")))],
        loc="lower left", frameon=False, fontsize=8, ncol=3,
        handletextpad=0.4, columnspacing=1.2, labelcolor=C_INK,
    )

    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(str(out_pdf).rsplit(".", 1)[0] + ".png", dpi=200, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------- main

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--bins", required=True, help="mosdepth regions.bed.gz")
    p.add_argument("--fasta", required=True)
    p.add_argument("--vcf", help="VCF with AD, for B-allele frequency")
    p.add_argument("--bin-size", type=int, default=100000)
    p.add_argument("--out-figure", required=True)
    p.add_argument("--out-table", required=True)
    p.add_argument("--gc-model", choices=["local", "2d", "poly"], default="local",
                   help="local: bin GC only (default, original model). "
                        "2d: local x regional GC cells, leave-one-chromosome-out. "
                        "poly: smooth quadratic surface in local and regional GC, "
                        "leave-one-chromosome-out. Run all three as a sensitivity "
                        "analysis: a call that depends on the model is not a finding.")
    a = p.parse_args()

    rng = np.random.default_rng(0)

    chrom, start, end, depth = read_bins(a.bins)
    if len(chrom) == 0:
        sys.exit("No bins on primary chromosomes - check contig naming in the mosdepth output.")
    print(f"{len(chrom)} bins, median depth {np.median(depth):.1f}x", file=sys.stderr)

    gc, nf = gc_fraction(a.fasta, chrom, start, end)

    auto = np.isin(chrom, AUTOSOMES)
    keep = auto & np.isfinite(gc) & (nf < 0.10) & (depth > 0)
    # drop extreme bins (segmental duplications, collapsed repeats)
    if keep.sum() > 100:
        lo, hi = np.percentile(depth[keep], [0.5, 99.5])
        keep &= (depth >= lo) & (depth <= hi)
    print(f"{keep.sum()} bins pass GC/N/outlier filters", file=sys.stderr)

    rgc = regional_gc(chrom, gc)
    keep &= np.isfinite(rgc)

    local_only = gc_correct(depth, gc, keep)
    if a.gc_model == "2d":
        corrected = gc_correct_2d_loo(depth, gc, rgc, chrom, keep)
    elif a.gc_model == "poly":
        corrected = gc_correct_poly_loo(depth, gc, rgc, chrom, keep)
    else:
        corrected = local_only
    baseline = np.median(corrected[keep])
    if baseline <= 0:
        sys.exit("Median corrected depth is zero - nothing to normalise against.")
    ratio = corrected / baseline

    # report GC correction effectiveness, at bin and at chromosome scale
    r_raw = np.corrcoef(gc[keep], depth[keep])[0, 1]
    r_cor = np.corrcoef(gc[keep], corrected[keep])[0, 1]
    print(f"bin depth~GC correlation: {r_raw:+.3f} raw -> {r_cor:+.3f} corrected",
          file=sys.stderr)
    rho_local = chrom_gc_rho(chrom, gc, local_only / np.median(local_only[keep]), keep)
    rho_used = chrom_gc_rho(chrom, gc, ratio, keep)
    print(f"chromosome GC~ratio Spearman: {rho_local:+.2f} (local-GC model) -> "
          f"{rho_used:+.2f} ({a.gc_model} model)", file=sys.stderr)

    baf_pos, baf_val = ({c: np.array([]) for c in CHROMS},
                        {c: np.array([]) for c in CHROMS})
    if a.vcf:
        baf_pos, baf_val = read_baf(a.vcf)
        print(f"{sum(len(v) for v in baf_val.values())} heterozygous SNVs for BAF",
              file=sys.stderr)

    rows, floor = analyse(chrom[keep], ratio[keep], baf_val, rng)

    # QC gate. With 22 autosomes, |Spearman rho| >= 0.43 is p < 0.05: chromosome
    # GC content then predicts the depth ratio, and calls on GC-rich chromosomes
    # cannot be distinguished from residual GC / replication-timing bias. Calls
    # are kept in the table but marked unconfirmed and drawn neutral.
    RHO_GATE = 0.43
    gc_confounded = abs(rho_used) >= RHO_GATE
    for r in rows:
        r["qc"] = ("gc_confounded" if (r["significant"] and gc_confounded)
                   else ("pass" if r["significant"] else "-"))
    if gc_confounded:
        print(f"QC: chromosome GC~ratio |rho| {abs(rho_used):.2f} >= {RHO_GATE}; "
              f"all calls marked UNCONFIRMED", file=sys.stderr)

    with open(a.out_table, "w") as fh:
        fh.write("chrom\tn_bins\tmedian_ratio\tci_lo\tci_hi\tdirection\t"
                 "mosaic_fraction\tbaf_median_abs_dev\tn_het_snv\tsignificant\tqc\n")
        for r in rows:
            fh.write(f"{r['chrom']}\t{r['n_bins']}\t{r['median_ratio']:.5f}\t"
                     f"{r['ci_lo']:.5f}\t{r['ci_hi']:.5f}\t{r['direction']}\t"
                     f"{r['mosaic_fraction']:.4f}\t{r['baf_median_abs_dev']:.5f}\t"
                     f"{r['n_het_snv']}\t{r['significant']}\t{r['qc']}\n")

    # Detection limit: the larger of the inter-chromosome noise floor and the
    # widest per-chromosome CI, expressed as a mosaic fraction.
    widths = [r["ci_hi"] - r["ci_lo"] for r in rows
              if r["chrom"] in AUTOSOMES and np.isfinite(r["ci_lo"])]
    limit = 2 * max(floor, (max(widths) / 2) if widths else 0.0)
    print(f"detection limit for a mosaic fraction: ~{limit*100:.1f}% of cells "
          f"(noise floor {floor:.4f} in ratio units)", file=sys.stderr)
    with open(str(a.out_table) + ".limit", "w") as fh:
        fh.write(f"noise_floor_ratio\t{floor:.6f}\n")
        fh.write(f"detection_limit_mosaic_fraction\t{limit:.6f}\n")
        fh.write(f"gc_model\t{a.gc_model}\n")
        fh.write(f"chrom_gc_rho_local_model\t{rho_local:.4f}\n")
        fh.write(f"chrom_gc_rho_used_model\t{rho_used:.4f}\n")

    n_sig = sum(1 for r in rows if r["significant"])
    print(f"{n_sig} chromosome(s) called above the floor", file=sys.stderr)

    plot(chrom[keep], start[keep], ratio[keep], baf_pos, baf_val, rows,
         a.out_figure, a.bin_size)
    print(f"wrote {a.out_figure} and {a.out_table}", file=sys.stderr)


if __name__ == "__main__":
    main()
