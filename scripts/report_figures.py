#!/usr/bin/env python3
"""
Figures for the Track 2 report (docs/smarin00_track2_report.md), drawn from
the committed result tables in findings/ only. No patient data is read.

  fig1_overview.svg     genotype -> mechanism -> goal -> outputs
  fig2_mechanism.svg    calibration alleles under ESM-2, AlphaMissense, burial
  fig3_screen.svg       direction-aware Open Targets screen, by panel gene
  fig4_oncology.svg     candidate x chemotherapy flags from FDA labels
  fig6_clinvar.svg      atlas bins: all reachable substitutions vs ClinVar
  fig7_related.svg      constraint and NMD maps of the related MVA genes

Usage: python3 scripts/report_figures.py [--root .] [--out docs/figures]
"""
import argparse
import csv
import os
import re
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                     # noqa: E402
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle  # noqa: E402
import numpy as np                                  # noqa: E402

# Reference palette (first three categorical slots validate all-pairs) and ink
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
WASH, PANEL = "#f3f2ee", "#fbfaf8"
W = 6.85                                            # A4 text width, inches

plt.rcParams.update({
    "font.family": "Inter", "font.size": 7.5, "svg.fonttype": "none",
    "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "axes.labelcolor": INK2,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2,
    "ytick.labelcolor": INK2, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "white", "axes.facecolor": "white",
})


def tsv(path):
    with open(path) as fh:
        rows = [l for l in fh if not l.startswith("#")]
    return list(csv.DictReader(rows, delimiter="\t"))


def save(fig, out, name):
    fig.savefig(os.path.join(out, name), bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


# ------------------------------------------------------------------ figure 1
def fig1(root, out):
    fig = plt.figure(figsize=(W, 3.05))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 44)
    ax.axis("off")

    def box(x, y, w, h, title, body, edge=AXIS, fill=PANEL):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.0",
                                    fc=fill, ec=edge, lw=0.7))
        ax.text(x + 1.3, y + h - 1.4, title, ha="left", va="top", fontsize=7.4,
                fontweight="semibold", color=INK, linespacing=1.2)
        nt = title.count("\n") + 1
        ax.text(x + 1.3, y + h - 1.6 - 2.9 * nt, body, ha="left", va="top", fontsize=6.6,
                color=INK2, linespacing=1.35)

    for x, t in [(1, "GENOTYPE (TRACK 1)"), (26, "MECHANISM"), (51, "THERAPEUTIC GOAL"),
                 (75.5, "OUTPUTS")]:
        ax.text(x, 42.6, t, fontsize=6.4, color=MUTED, fontweight="semibold", va="center")

    box(1, 21, 22, 19, "p.Leu737*", "stop-gain in exon 17 of 23\nchr15:40209701 T>G")
    box(1, 1, 22, 19, "p.Asn1002Lys", "missense in the pseudo-\nkinase domain, 10 residues\nfrom MVA allele L1012P\nchr15:40220612 T>G")
    box(26, 21, 22, 19, "Null allele", "stop lies 747 nt upstream\nof the last exon junction:\nmRNA decay (NMD),\nno protein")
    box(26, 1, 22, 19, "Unstable but functional",
        "ESM-2 calibrated on six\nmeasured alleles; AlphaFold\ncontacts; AlphaMissense;\npublished biology")
    box(51, 1, 22, 39, "Too little BubR1,\nnot broken BubR1",
        "All of the child's BubR1\ncomes from the unstable\nallele.\n\n"
        "Segregation fails at about\n6% of normal BubR1 and is\nnear normal at about 13%.\n\n"
        "Goal: raise the level of the\nprotein the child already\nmakes, about two-fold.",
        edge=INK2, fill=WASH)
    outs = [("Candidates", "niacin, arimoclomol;\ndasatinib + quercetin"),
            ("Do-no-harm list", "HSP90, NAMPT, Aurora B\nand MPS1 inhibitors"),
            ("Cancer-treatment check", "candidates vs chemotherapy,\nfrom FDA label sentences"),
            ("Resource for families", "atlas of every BUB1B variant;\nsix related MVA genes")]
    for i, (t, b) in enumerate(outs):
        box(75.5, 31 - 10 * i, 23.5, 9, t, b)

    def arrow(x0, y0, x1, y1):
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=7,
                                     lw=0.8, color=INK2, shrinkA=0, shrinkB=0))
    arrow(23.3, 30.5, 25.7, 30.5)
    arrow(23.3, 10.5, 25.7, 10.5)
    arrow(48.3, 30.5, 50.7, 25.0)
    arrow(48.3, 10.5, 50.7, 16.0)
    arrow(73.3, 20.5, 75.2, 20.5)
    save(fig, out, "fig1_overview.svg")


# ------------------------------------------------------------------ figure 2
def fig2(root, out):
    t2 = os.path.join(root, "findings/track2")
    esm = {r["variant"]: (r["class"], float(r["esm_llr"])) for r in tsv(f"{t2}/esm_scores.tsv")}
    st = {r["variant"]: r for r in tsv(f"{t2}/structure_context.tsv")}
    am = {}
    for r in tsv(f"{t2}/alphamissense/am_vs_esm.tsv"):
        am[f"{r['wt']}{r['pos']}{r['mt']}"] = float(r["am_pathogenicity"])
    order = ["N1002K", "L1012P", "I909T", "L844F", "R727C", "Q921H", "R550Q"]
    ys = {v: len(order) - i for i, v in enumerate(order)}
    cls = {v: esm[v][0] for v in order}

    fig, axes = plt.subplots(1, 3, figsize=(W, 2.5), sharey=True,
                             gridspec_kw={"wspace": 0.14})

    def head(ax, letter, title, sub):
        ax.text(0, 1.17, f"{letter}   {title}", transform=ax.transAxes, fontsize=7.4,
                fontweight="semibold", color=INK, va="bottom")
        ax.text(0, 1.06, f"      {sub}", transform=ax.transAxes, fontsize=6.6, color=INK2,
                va="bottom")
    d_max = max(esm[v][1] for v in order if cls[v] == "destabilising")
    s_min = min(esm[v][1] for v in order if cls[v] == "stable")

    def marks(ax, val, note=None):
        for v in order:
            x = val(v)
            if x is None:
                continue
            y = ys[v]
            c = cls[v]
            if c == "PATIENT":
                ax.scatter(x, y, marker="*", s=150, c=INK, edgecolors="white", linewidths=1.0, zorder=4)
            elif c == "destabilising":
                ax.scatter(x, y, s=34, c=ORANGE, edgecolors="white", linewidths=1.2, zorder=3)
            else:
                ax.scatter(x, y, s=34, c=BLUE, edgecolors="white", linewidths=1.2, zorder=3)
            if note and v in note:
                ax.text(x - 0.04, y, note[v], ha="right", va="center", fontsize=5.8, color=MUTED,
                        bbox=dict(fc="white", ec="none", pad=0.6))
        for y in ys.values():
            ax.axhline(y, color=GRID, lw=0.5, zorder=0)
        ax.set_ylim(0.3, len(order) + 0.8)

    # a: ESM-2
    a = axes[0]
    a.axvspan(-9, d_max, color=ORANGE, alpha=0.07, lw=0)
    a.axvspan(s_min, 1.5, color=BLUE, alpha=0.07, lw=0)
    marks(a, lambda v: esm[v][1])
    a.set_xlim(-8, 1.5)
    a.set_xlabel("ESM-2 log-likelihood ratio")
    head(a, "a", "ESM-2, calibrated", "measured classes separate")
    a.text(-7.8, 7.62, "like destabilisers", fontsize=5.8, color=INK2, va="center")
    a.text(1.4, 7.62, "like stable", fontsize=5.8, color=INK2, va="center", ha="right")

    # b: AlphaMissense (class cut-offs from Cheng et al. 2023)
    b = axes[1]
    b.axvspan(0, 0.34, color=BLUE, alpha=0.07, lw=0)
    b.axvspan(0.564, 1, color=ORANGE, alpha=0.07, lw=0)
    marks(b, lambda v: am.get(v))
    b.set_xlim(-0.03, 1.03)
    b.set_xticks([0, 0.25, 0.5, 0.75, 1])
    b.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"])
    b.set_xlabel("AlphaMissense pathogenicity")
    head(b, "b", "AlphaMissense", "measured classes overlap")
    b.text(0.02, 7.62, "likely benign", fontsize=5.8, color=INK2, va="center")
    b.text(0.98, 7.62, "likely pathogenic", fontsize=5.8, color=INK2, va="center", ha="right")

    # c: burial
    c = axes[2]
    c.axvspan(0, 0.20, color=ORANGE, alpha=0.07, lw=0)
    rsa = {v: float(st[v]["rsa"]) for v in order}
    conf = {v: float(st[v]["plddt"]) >= 70 for v in order}
    marks(c, lambda v: rsa[v], note={v: f"low confidence (pLDDT {float(st[v]['plddt']):.0f})"
                                     for v in order if not conf[v]})
    c.set_xlim(-0.03, 1.03)
    c.set_xticks([0, 0.25, 0.5, 0.75, 1])
    c.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"])
    c.set_xlabel("Relative solvent accessibility")
    head(c, "c", "Burial in the AlphaFold model", "measured classes overlap")
    c.text(0.02, 7.62, "buried (< 0.20)", fontsize=5.8, color=INK2, va="center")

    a.set_yticks(list(ys.values()))
    a.set_yticklabels([("N1002K (patient)" if v == "N1002K" else v) for v in order])
    for lab in a.get_yticklabels():
        if lab.get_text().startswith("N1002K"):
            lab.set_fontweight("semibold")
            lab.set_color(INK)
    for ax in axes:
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)

    h = [plt.Line2D([], [], marker="o", ls="", ms=5.5, mfc=ORANGE, mec="white", label="destabilising (measured)"),
         plt.Line2D([], [], marker="o", ls="", ms=5.5, mfc=BLUE, mec="white", label="stable (measured)"),
         plt.Line2D([], [], marker="*", ls="", ms=10, mfc=INK, mec="white", label="patient allele")]
    for ax in axes:
        ax.tick_params(axis="x", labelsize=6.4)
        ax.xaxis.label.set_size(6.8)
    fig.legend(handles=h, loc="upper center", ncol=3, frameon=False, fontsize=6.6,
               bbox_to_anchor=(0.5, -0.02), handletextpad=0.3, columnspacing=1.6)
    save(fig, out, "fig2_mechanism.svg")


# ------------------------------------------------------------------ figure 3
SALT = re.compile(r"\s+(HYDROCHLORIDE|LACTATE|MESYLATE|CHLORIDE|SODIUM|SULFATE)$")


def parent(name):
    return SALT.sub("", name.upper()).strip()


def fig3(root, out):
    import yaml
    cfg = yaml.safe_load(open(os.path.join(root, "config/genes/bub1b.yaml")))
    rows = tsv(os.path.join(root, "findings/track2/opentargets/ot_target_drugs.tsv"))
    good, bad = defaultdict(set), defaultdict(set)
    for r in rows:
        if r["verdict"] == "candidate":
            good[r["target"]].add(parent(r["drug"]))
        elif r["verdict"] == "do_no_harm":
            bad[r["target"]].add(parent(r["drug"]))
    want = {"up": "raise", "down": "lower", "no_inhibit": "do not inhibit"}
    outcome = {
        "BUB1B": "no drug in Open Targets",
        "SIRT2": "no drug; reached through NAD\u207a (niacin)",
        "NAMPT": "do-no-harm: NAMPT inhibitors deplete NAD\u207a",
        "HSF1": "no drug; reached through arimoclomol",
        "HSPA1A": "do-no-harm: CD24Fc (annotated HSP70/90 inhibitor)",
        "CREBBP": "rejected: PRI-724 blocks CBP\u2013\u03b2-catenin, not K668",
        "HDAC2": "not ranked: acts on mitotic timing, not level",
        "HDAC3": "not ranked: as HDAC2",
        "HSP90AA1": "do-no-harm: mutant BubR1 is an HSP90 client",
        "HSP90AB1": "do-no-harm: as HSP90AA1",
        "TTK": "do-no-harm: weakens the checkpoint",
        "AURKB": "do-no-harm: impairs error correction",
        "BCL2L1": "rejected: navitoclax, obatoclax (thrombocytopenia)",
    }
    panel = cfg["pathway"]
    n = len(panel)
    fig = plt.figure(figsize=(W * 0.92, 0.215 * n + 0.75))
    ax = fig.add_axes([0.205, 0.16, 0.33, 0.74])
    from matplotlib.transforms import blended_transform_factory as blend
    tr = blend(ax.transAxes, ax.transData)
    for i, p in enumerate(panel):
        y = n - i
        g, b_ = len(good.get(p["symbol"], ())), len(bad.get(p["symbol"], ()))
        if g:
            ax.barh(y, g, height=0.58, color=BLUE, lw=0)
            ax.text(g + 0.4, y, str(g), va="center", fontsize=6.4, color=INK2)
        if b_:
            ax.barh(y, -b_, height=0.58, color=ORANGE, lw=0)
            ax.text(-b_ - 0.4, y, str(b_), va="center", ha="right", fontsize=6.4, color=INK2)
        if not g and not b_:
            ax.text(0.4, y, "0", va="center", fontsize=6.4, color=MUTED)
        ax.text(-0.62, y, p["symbol"], transform=tr, va="center", ha="left", fontsize=7.2,
                color=INK, fontweight="semibold")
        ax.text(-0.33, y, want[p["desired"]], transform=tr, va="center", ha="left",
                fontsize=6.8, color=INK2)
        ax.text(1.04, y, outcome.get(p["symbol"], ""), transform=tr, va="center", ha="left",
                fontsize=6.7, color=INK2)
        if i:
            ax.plot([-0.62, 2.95], [y + 0.5, y + 0.5], transform=tr, color=GRID, lw=0.4,
                    clip_on=False, zorder=0)
    ax.axvline(0, color=AXIS, lw=0.8)
    ax.set_xlim(-20, 14)
    ax.set_ylim(0.45, n + 0.55)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xticks([-18, -12, -6, 0, 6, 12])
    ax.set_xticklabels(["18", "12", "6", "0", "6", "12"])
    ax.tick_params(axis="x", labelsize=6.2)
    ax.set_xlabel("compounds (parent molecules; salts merged)", fontsize=6.4)
    yh = n + 1.05
    for x, t, ha in [(-0.62, "PANEL GENE", "left"), (-0.33, "DESIRED", "left"),
                     (0.585, "\u2190 ACTS AGAINST", "right"), (0.605, "IN DIRECTION \u2192", "left"),
                     (1.04, "CURATION OUTCOME", "left")]:
        ax.text(x, yh, t, transform=tr, fontsize=5.9, color=MUTED, fontweight="semibold",
                ha=ha, va="center", clip_on=False)
    save(fig, out, "fig3_screen.svg")
    return good, bad


# ------------------------------------------------------------------ figure 4
def fig4(root, out):
    import yaml
    cfg = yaml.safe_load(open(os.path.join(root, "config/genes/bub1b.yaml")))["oncology"]
    comp = tsv(os.path.join(root, "findings/track2/oncology/compatibility.tsv"))
    labels = {r["drug"]: r["status"] for r in
              tsv(os.path.join(root, "findings/track2/oncology/labels.tsv"))}
    chemo = []
    regs = defaultdict(list)
    for reg, drugs in cfg["regimens"].items():
        for d in drugs:
            if d not in chemo:
                chemo.append(d)
            regs[d].append(reg.replace("RMS ", "").replace("Wilms ", ""))
    cands = ["niacin", "arimoclomol", "dasatinib", "valproic acid", "quercetin",
             "nicotinamide riboside"]
    flags = defaultdict(list)
    for r in comp:
        f = "CYP3A4" if r["check"] == "cyp3a4" else (r["organ"] if r["check"] == "toxicity"
                                                       else "named")
        flags[(r["candidate"], r["chemotherapy"])].append(f)
    tint = ["#ffffff", "#fbe6dc", "#f7cdb9", "#f3b496"]
    order = ["CYP3A4", "marrow", "liver", "heart/QT", "nerves", "kidney", "pancreas", "named"]
    nC, nR = len(chemo), len(cands)
    nolab = {c: labels.get(c, "").startswith("no US") for c in cands}
    nl = {c: max([len(set(flags.get((c, d), []))) for d in chemo] + [1]) for c in cands}
    rh = {c: (0.62 if nolab[c] else 0.55 + 0.33 * nl[c]) for c in cands}    # row heights
    total = sum(rh.values())
    fig, ax = plt.subplots(figsize=(W, 0.36 * total + 0.55))
    x0 = 1.45
    for j, d in enumerate(chemo):
        ax.text(x0 + j + 0.5, total + 0.42, d, ha="center", va="bottom", fontsize=6.9,
                fontweight="semibold", color=INK)
        ax.text(x0 + j + 0.5, total + 0.12, " \u00b7 ".join(regs[d]) if len(regs[d]) < 4
                else "all regimens", ha="center", va="bottom", fontsize=5.6, color=MUTED)
    y = total
    for c in cands:
        h = rh[c]
        y -= h
        name = "NR" if c == "nicotinamide riboside" else c
        ax.text(x0 - 0.08, y + h / 2, name, ha="right", va="center", fontsize=6.9,
                fontweight="semibold", color=INK)
        for j, d in enumerate(chemo):
            f = sorted(set(flags.get((c, d), [])), key=order.index)
            fill = WASH if nolab[c] else tint[min(len(f), 3)]
            edge = GRID if (not f and not nolab[c]) else "none"
            ax.add_patch(Rectangle((x0 + j + 0.03, y + 0.05), 0.94, h - 0.1, fc=fill, ec=edge,
                                   lw=0.5))
            if nolab[c]:
                continue
            ax.text(x0 + j + 0.5, y + h / 2, "\n".join(f) if f else "\u2013", ha="center",
                    va="center", fontsize=6.0 if f else 7, color=INK if f else MUTED,
                    linespacing=1.25)
        if nolab[c]:
            ax.text(x0 + nC / 2, y + h / 2, "no US prescription label: not assessable from labels",
                    ha="center", va="center", fontsize=6.2, color=INK2)
    ax.set_xlim(0, x0 + nC + 0.02)
    ax.set_ylim(-0.02, total + 0.95)
    ax.axis("off")
    save(fig, out, "fig4_oncology.svg")


# ------------------------------------------------------------------ figure 6
def fig6(root, out):
    t2 = os.path.join(root, "findings/track2")
    am = tsv(f"{t2}/alphamissense/am_vs_esm.tsv")
    cv = tsv(f"{t2}/clinvar_validation.tsv")
    test = [r for r in cv if r["calibration_allele"] == "False" and r["esm_llr"]]

    def uniq(cls):
        u = {}
        for r in test:
            if r["class"] == cls:
                u[r["variant"]] = r["atlas_bin"]
        return list(u.values())
    groups = [("All single-nucleotide-reachable\nmissense changes", [r["esm_class"] for r in am]),
              ("ClinVar benign / likely benign", uniq("B/LB")),
              ("ClinVar uncertain / conflicting", uniq("VUS/conflicting"))]
    bins = [("like_destabilisers", "like destabilisers", ORANGE),
            ("indeterminate", "indeterminate", "#d6d5cf"),
            ("like_tolerated", "like tolerated", BLUE)]
    fig, ax = plt.subplots(figsize=(W * 0.9, 1.55))
    for i, (name, vals) in enumerate(groups):
        y = len(groups) - 1 - i
        n = len(vals)
        left = 0
        for key, lab, col in bins:
            k = sum(v == key for v in vals)
            w = 100 * k / n
            ax.barh(y, w - 0.35, left=left, height=0.55, color=col, lw=0)
            txt = f"{w:.0f}%" if w >= 6 else ""
            if 0 < w < 6:
                ax.text(left + w + 0.8, y, f"{w:.0f}%", ha="left", va="center", fontsize=6.6,
                        color=INK, fontweight="semibold")
            if txt:
                ax.text(left + w / 2, y, txt, ha="center", va="center", fontsize=6.6,
                        color="white" if col != "#d6d5cf" else INK, fontweight="semibold")
            left += w
        ax.text(-1.5, y, f"{name}  (n = {n:,})", ha="right", va="center", fontsize=6.7,
                color=INK, linespacing=1.15)
    ax.set_xlim(0, 100)
    ax.set_ylim(-0.6, len(groups) - 0.35)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xticklabels(["0", "25", "50", "75", "100%"])
    h = [Rectangle((0, 0), 1, 1, color=c, label=l) for _, l, c in bins]
    ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, -0.3), ncol=3,
              frameon=False, fontsize=6.4, handlelength=1.2)
    save(fig, out, "fig6_clinvar.svg")


# ------------------------------------------------------------------ figure 7
def fig7(root, out):
    rel = tsv(os.path.join(root, "findings/track2/related_genes.tsv"))
    dis = defaultdict(list)
    for r in rel:
        if r["selected"] == "yes":
            m = re.search(r"syndrome (\d+)", r["disease_name"])
            if r["route"] == "series" and m:
                dis[r["gene"]].append(f"MVA{m.group(1)} ({r['disease_id']})")
            elif r["route"] == "orphanet" and not any("MVA" in x for x in dis[r["gene"]]):
                pass
    genes = ["CEP57", "TRIP13", "CENATAC", "MAD1L1", "BUB1", "BUB3"]
    fig, axes = plt.subplots(3, 2, figsize=(W, 3.9), gridspec_kw={"hspace": 0.75, "wspace": 0.12})
    for ax, g in zip(axes.flat, genes):
        d = os.path.join(root, "findings/genes", g)
        mis = tsv(f"{d}/atlas_missense.tsv")
        tr = tsv(f"{d}/atlas_truncating.tsv")
        L = max(int(r["pos"]) for r in mis)
        s = np.zeros(L)
        k = np.zeros(L)
        for r in mis:
            p = int(r["pos"]) - 1
            s[p] += float(r["esm_llr"])
            k[p] += 1
        mean = s / np.maximum(k, 1)
        win = 9
        sm = (np.convolve(mean, np.ones(win), mode="same")
              / np.convolve(np.ones(L), np.ones(win), mode="same"))     # edge-corrected
        esc = [int(r["codon"]) for r in tr if not r["outcome"].startswith("null")]
        if esc:
            ax.axvspan(min(esc), L, color=AQUA, alpha=0.18, lw=0)
        ax.axhline(0, color=GRID, lw=0.5)
        ax.plot(np.arange(1, L + 1), sm, color=INK, lw=0.9, solid_joinstyle="round")
        ax.set_xlim(1, L)
        ax.set_ylim(-13, 1.2)
        ax.set_yticks([-12, -8, -4, 0])
        tag = ", ".join(dis.get(g, [])) or "Orphanet MVA group (ORPHA:1052)"
        ax.set_title(f"{g}", loc="left", fontsize=7.4, fontweight="semibold", color=INK, pad=3)
        ax.text(1.0, 1.04, f"{tag} · {L} aa", transform=ax.transAxes, ha="right",
                va="bottom", fontsize=6.0, color=INK2)
        ax.tick_params(labelsize=6.2)
    for ax in axes[:, 0]:
        ax.set_ylabel("ESM-2 LLR", fontsize=6.4)
    for ax in axes[-1, :]:
        ax.set_xlabel("residue / codon", fontsize=6.4)
    for ax in axes[:, 1]:
        ax.set_yticklabels([])
    h = [plt.Line2D([], [], color=INK, lw=0.9, label="ESM-2 constraint (mean of 19 substitutions, 9-residue window)"),
         Rectangle((0, 0), 1, 1, color=AQUA, alpha=0.35, lw=0, label="a premature stop here escapes NMD")]
    fig.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=2, frameon=False,
               fontsize=6.4, handlelength=1.6)
    save(fig, out, "fig7_related.svg")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", default=".")
    p.add_argument("--out", default="docs/figures")
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    fig1(a.root, a.out)
    fig2(a.root, a.out)
    good, bad = fig3(a.root, a.out)
    fig4(a.root, a.out)
    fig6(a.root, a.out)
    fig7(a.root, a.out)
    g = set().union(*good.values()) if good else set()
    b = set().union(*bad.values()) if bad else set()
    print(f"screen: {len(g)} candidate compounds, {len(b)} do-no-harm compounds")
    print(f"wrote figures to {a.out}")


if __name__ == "__main__":
    main()
