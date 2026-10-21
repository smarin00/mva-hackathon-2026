#!/usr/bin/env python3
"""
Variant atlas: every possible missense change and every possible premature
stop in one protein, each mapped to a predicted mechanism and, where the gene
is calibrated, the therapeutic direction that mechanism implies.

Gene-agnostic: everything case-specific comes from a config file
(config/genes/<gene>.yaml). Built first for BUB1B / MVA1, where every reported
patient carries one truncating and one missense allele (Hanks 2004), so a new
family always needs the same two answers:

  * Truncating allele: is the transcript degraded (null), or does it escape
    NMD -- and if it escapes, does the truncated protein keep the features
    the protein needs?
  * Missense allele: does it look like experimentally verified stability
    mutants (where raising protein level restores function), or does it hit
    an essential functional site (where more protein may not help)?

Missense: ESM-2 masked-marginal log-likelihood ratio for all 19 substitutions
at every position. With `calibration` alleles in the config, classes are set
from them, not from a fixed cut-off:
  like_destabilisers  LLR <= the weakest verified destabiliser
  like_tolerated      LLR >= the most damaging verified stable variant
  indeterminate       in between
Without calibration, every change is 'uncalibrated' and no therapy direction
is given: constraint is reported, mechanism is not inferred.

Truncating: NMD geometry (50/55-nt rule, last-exon escape, start-proximal
caveat) at every codon, from the Ensembl transcript structure, plus which
curated essential features the truncated product would keep. Without a
curated domain map, domains come from UniProt and essentiality is reported as
not assessed.

Not modelled: frameshifts (their PTC lies downstream of the variant -- run the
PTC position), splice variants, in-frame indels, complexes.

The ESM scan (one forward pass per residue) is cached, so reruns of the
tables and figure take seconds.

Usage (or `make atlas`, `make gene-atlas G=CEP57`)
  variant_atlas.py --config config/genes/bub1b.yaml --fasta resources/O60566.fasta \\
      --gtf /opt/refs/Homo_sapiens.GRCh38.116.gtf.gz \\
      --pdb resources/AF-O60566-F1-model.pdb --outdir findings/track2/atlas
"""

import argparse
import hashlib
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import caseconfig                    # noqa: E402
import esm_variant_score as evs      # noqa: E402
import nmd_check as nmd              # noqa: E402

AA = evs.AA

# Categorical slots 1-3 of the reference palette, validated all-pairs
# (light mode). Aqua is below 3:1 on the surface, so every segment is also
# direct-labelled and the TSV is the table view.
C_NULL, C_NONFUNC, C_FUNC = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb"


# --------------------------------------------------------------------- ESM
def esm_logprobs(seq, model_name, threads, cache):
    key = hashlib.sha1((seq + model_name).encode()).hexdigest()
    if cache and os.path.exists(cache):
        z = np.load(cache, allow_pickle=False)
        if str(z["key"]) == key:
            print(f"ESM: loaded cached scan from {cache}", file=sys.stderr)
            return z["lp"]
        print("ESM: cache is for a different sequence/model; recomputing",
              file=sys.stderr)

    import torch
    import esm
    torch.set_num_threads(threads)
    model, alphabet = getattr(esm.pretrained, model_name)()
    model.eval()
    conv = alphabet.get_batch_converter()
    idx = [alphabet.get_idx(a) for a in AA]

    L = len(seq)
    lp = np.zeros((L, len(AA)), dtype=np.float32)
    for pos in range(1, L + 1):
        s0 = evs.window(seq, pos)
        _, _, toks = conv([("x", seq[s0:s0 + evs.WINDOW])])
        i = pos - 1 - s0 + 1
        toks[0, i] = alphabet.mask_idx
        with torch.no_grad():
            out = torch.log_softmax(model(toks)["logits"][0, i], dim=-1)
        lp[pos - 1] = out[idx].numpy()
        if pos % 100 == 0:
            print(f"  ESM {pos}/{L}", file=sys.stderr)
    if cache:
        np.savez_compressed(cache, lp=lp, key=key)
    return lp


# --------------------------------------------------------------- structure
def per_residue_structure(pdb, L):
    plddt = np.full(L, np.nan)
    rsa = np.full(L, np.nan)
    if not pdb:
        return plddt, rsa
    from Bio.PDB import PDBParser
    from Bio.PDB.SASA import ShrakeRupley
    import structure_context as sc
    s = PDBParser(QUIET=True).get_structure("m", pdb)
    ShrakeRupley().compute(s[0], level="R")
    for res in next(s[0].get_chains()):
        i = res.id[1]
        if 1 <= i <= L and res.get_resname() in sc.MAX_ASA:
            plddt[i - 1] = np.mean([a.get_bfactor() for a in res])
            rsa[i - 1] = res.sasa / sc.MAX_ASA[res.get_resname()]
    return plddt, rsa


# ----------------------------------------------------------------- features
UNIPROT_TYPES = {"Domain", "Region", "Motif", "Repeat", "Zinc finger",
                 "Coiled coil", "DNA binding"}


def uniprot_features(acc, out_tsv):
    """Domain annotations from UniProt, in the curated-map format, with
    essentiality 'untested'. Written to out_tsv for provenance."""
    import json
    import urllib.request
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
    with urllib.request.urlopen(url, timeout=60) as r:
        entry = json.load(r)
    rows = []
    for f in entry.get("features", []):
        if f.get("type") not in UNIPROT_TYPES:
            continue
        loc = f.get("location", {})
        s, e = loc.get("start", {}).get("value"), loc.get("end", {}).get("value")
        if s is None or e is None:
            continue
        name = (f.get("description") or f["type"]).replace("\t", " ")
        rows.append({"start": int(s), "end": int(e), "feature": name,
                     "role": f["type"], "dispensable": "untested",
                     "source": f"UniProt {acc}"})
    with open(out_tsv, "w") as fh:
        fh.write(f"# Domain annotations from UniProt {acc} ({url}); essentiality not assessed\n")
        fh.write("start\tend\tfeature\trole\tdispensable\tsource\n")
        for d in rows:
            fh.write(f"{d['start']}\t{d['end']}\t{d['feature']}\t{d['role']}\t"
                     f"{d['dispensable']}\t{d['source']}\n")
    return rows


def feature_at(domains, pos):
    hits = [d for d in domains if d["start"] <= pos <= d["end"]]
    return min(hits, key=lambda d: d["end"] - d["start"]) if hits else None


def essential_features(domains):
    return [d for d in domains if d["dispensable"] == "no"]


# ---------------------------------------------------------------- missense
def missense_rows(seq, lp, domains, plddt, rsa, thresholds, stability_features):
    d_max, s_min = thresholds if thresholds else (None, None)
    flat = []
    rows = []
    for p in range(1, len(seq) + 1):
        wt = seq[p - 1]
        wi = AA.index(wt)
        feat = feature_at(domains, p)
        for mt in AA:
            if mt == wt:
                continue
            llr = float(lp[p - 1, AA.index(mt)] - lp[p - 1, wi])
            flat.append(llr)
            if thresholds is None:
                cls = "uncalibrated"
            elif llr <= d_max:
                cls = "like_destabilisers"
            elif llr >= s_min:
                cls = "like_tolerated"
            else:
                cls = "indeterminate"

            fname = feat["feature"] if feat else "linker/other"
            ess = feat["dispensable"] if feat else "untested"
            if cls == "uncalibrated":
                direction = "n/a (gene not calibrated)"
            elif cls == "like_tolerated":
                direction = "likely tolerated"
            elif ess == "no":
                direction = "functional-site defect: raising dose may not compensate"
            elif fname in stability_features and cls == "like_destabilisers":
                direction = "stability defect: dose-raising strategy applies"
            else:
                direction = "uncertain"
            rows.append([p, wt, mt, llr, fname, ess, plddt[p - 1], rsa[p - 1],
                         cls, direction])
    flat = np.sort(np.array(flat))
    for r in rows:
        r.insert(4, 100.0 * np.searchsorted(flat, r[3], side="right") / len(flat))
    return rows


# --------------------------------------------------------------- truncating
def transcript_geometry(gtf, gene):
    tx = nmd.scan(gtf, gene=gene)
    tid, rec, how = nmd.pick(tx)
    strand = rec["strand"]
    exons = nmd.ordered(rec["exons"], strand)
    cds = nmd.ordered(rec["cds"], strand)
    first = cds[0][0] if strand == "+" else cds[0][1]
    cds_tx0 = nmd.to_transcript_coord(first, rec["exons"], strand)   # tx coord of CDS base 1
    bounds, acc = [], 0
    for s, e in exons:
        acc += e - s + 1
        bounds.append(acc)                                           # exon ends, tx coords
    cds_len = sum(e - s + 1 for s, e in cds)
    return {"tid": tid, "how": how, "bounds": bounds, "cds_tx0": cds_tx0,
            "n_aa": nmd.aa_count(rec, cds_len), "last_jct": bounds[-2]}


def truncating_rows(seq, geo, domains, curated=True):
    ess = essential_features(domains)
    rows = []
    for c in range(1, geo["n_aa"] + 1):
        cds_pos = 3 * (c - 1) + 1
        tpos = geo["cds_tx0"] + cds_pos - 1
        exon = next(i for i, b in enumerate(geo["bounds"], 1) if tpos <= b)
        dist = geo["last_jct"] - tpos
        n_ex = len(geo["bounds"])
        if exon == n_ex:
            nmd_call = "escapes (last exon)"
        elif dist > nmd.NMD_WINDOW_STRICT:
            nmd_call = "NMD"
        elif dist > nmd.NMD_WINDOW:
            nmd_call = "NMD (borderline)"
        else:
            nmd_call = "escapes (<=50 nt from last junction)"

        lost = [d["feature"] for d in ess if d["end"] >= c]
        if nmd_call.startswith("NMD"):
            outcome = "null (NMD)"
        elif not curated:
            outcome = "escapes NMD; features not assessed"
        elif lost:
            outcome = "escapes NMD; loses essential feature"
        else:
            outcome = "escapes NMD; keeps all essential features"
        rows.append([c, seq[c - 1] if c <= len(seq) else "?", cds_pos, tpos, exon,
                     dist, nmd_call, ",".join(lost) or "-", outcome,
                     "start-proximal: reinitiation possible" if cds_pos <= nmd.REINIT_WINDOW else ""])
    return rows


# ------------------------------------------------------------------- figure
LABEL_OFFSETS = {"I909T": (-8, 0), "L1012P": (8, 0), "L844F": (8, 6),
                 "R727C": (-8, 0), "Q921H": (8, 0), "R550Q": (8, 0)}


def plot(seq, lp, domains, trunc, calib, patient_mis, patient_ptc, geo, out_pdf,
         gene="", curated=True):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    L = len(seq)
    x = np.arange(1, L + 1)
    wi = np.array([AA.index(a) for a in seq])
    mean_llr = (lp.sum(axis=1) - lp[np.arange(L), wi] - 19 * lp[np.arange(L), wi]) / 19
    # = mean over the 19 substitutions of (lp[mt] - lp[wt])
    k = 9
    smooth = np.convolve(np.pad(mean_llr, k // 2, mode="edge"), np.ones(k) / k, "valid")

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(11, 5.6), sharex=True,
                                 gridspec_kw={"height_ratios": [2.4, 1], "hspace": 0.22})
    fig.patch.set_facecolor("white")

    # --- panel A: positional constraint + features + calibrated variants
    shade = [d for d in domains if d["end"] - d["start"] >= 5 or d["feature"].startswith(("KEN", "ABBA"))]
    for d in shade:
        a1.axvspan(d["start"] - 0.5, d["end"] + 0.5, color="#f1f0ec", lw=0, zorder=0)
    if curated:
        for d in domains:
            if d["dispensable"] == "no" or d["feature"] == "PSEUDOKINASE":
                req = d["dispensable"] == "no"
                a1.text((d["start"] + d["end"]) / 2, 1.9,
                        d["feature"] if req else "pseudokinase (dispensable)",
                        ha="center", va="bottom", fontsize=7.5,
                        color=INK if req else INK2, fontweight="bold" if req else "normal")
    else:
        # UniProt annotations: label the longer ones, staggered, never bold
        # (essentiality is not assessed).
        big = sorted([d for d in domains if d["end"] - d["start"] >= 0.04 * len(seq)],
                     key=lambda d: d["start"])
        for k, d in enumerate(big):
            a1.text((d["start"] + d["end"]) / 2, 1.9 - 0.7 * (k % 2), d["feature"][:28],
                    ha="center", va="bottom", fontsize=7, color=INK2)
    a1.plot(x, smooth, color=INK, lw=1.2, zorder=2)
    a1.axhline(0, color=GRID, lw=0.8, zorder=1)

    for v, cls in calib:
        wt, pos, mt = evs.parse(v)
        llr = lp[pos - 1, AA.index(mt)] - lp[pos - 1, AA.index(wt)]
        filled = cls == "destabilising"
        a1.scatter([pos], [llr], s=46, zorder=4, marker="o",
                   facecolor=INK2 if filled else SURFACE, edgecolor=INK2, linewidth=1.2)
        dx, dy = LABEL_OFFSETS.get(v, (6, 4))
        # A light backing keeps labels legible where the constraint line runs
        # through them (the real line is far from flat).
        a1.annotate(v, (pos, llr), xytext=(dx, dy), textcoords="offset points",
                    fontsize=7, color=INK2, ha="left" if dx >= 0 else "right",
                    va="center", zorder=6,
                    bbox=dict(boxstyle="round,pad=0.15", facecolor=SURFACE, edgecolor="none", alpha=0.9))
    if patient_mis:
        wt, pos, mt = evs.parse(patient_mis)
        llr = lp[pos - 1, AA.index(mt)] - lp[pos - 1, AA.index(wt)]
        # Patient in ink, not a palette hue: blue means "null" in panel B.
        a1.scatter([pos], [llr], s=170, marker="*", color=INK, edgecolor=SURFACE,
                   linewidth=1.2, zorder=5)
        a1.annotate(f"{patient_mis} (patient)", (pos, llr), xytext=(-10, -16),
                    textcoords="offset points", ha="right", va="center", fontsize=8,
                    color=INK, fontweight="bold", zorder=6,
                    bbox=dict(boxstyle="round,pad=0.15", facecolor=SURFACE, edgecolor="none", alpha=0.9))
    a1.set_ylabel("ESM-2 log-likelihood ratio\n(lower = less tolerated)", fontsize=8.5, color=INK)
    a1.set_ylim(min(-9, np.nanmin(smooth) - 1), 2.7)
    a1.set_title(f"{gene} variant atlas: missense constraint (top) and fate of a "
                 "premature stop at each codon (bottom)", loc="left", fontsize=10.5,
                 color=INK, pad=8)
    if curated:
        caption = ("Line: mean over 19 substitutions, 9-residue window. Filled "
                   "circles: verified destabilising alleles; open: verified stable.\n"
                   "Bold: features required for function; light bands: other annotated motifs.")
    else:
        caption = ("Line: mean over 19 substitutions, 9-residue window. Bands: UniProt "
                   "domain annotations (essentiality not assessed). No calibration alleles.")
    a1.text(0.0, -0.02, caption, transform=a1.transAxes, fontsize=7, color=MUTED, va="top")

    # --- panel B: truncation outcome strip
    colour = {"null (NMD)": C_NULL,
              "escapes NMD; loses essential feature": C_NONFUNC,
              "escapes NMD; keeps all essential features": C_FUNC,
              "escapes NMD; features not assessed": C_FUNC}
    runs, start = [], 0
    for i in range(1, len(trunc) + 1):
        if i == len(trunc) or trunc[i][8] != trunc[start][8]:
            runs.append((trunc[start][0], trunc[i - 1][0], trunc[start][8]))
            start = i
    for s, e, o in runs:
        a2.axvspan(s - 0.5, e + 0.5, ymin=0.35, ymax=0.85, color=colour[o], lw=0)
        if e - s > 25:
            short = {"null (NMD)": "null (NMD)",
                     "escapes NMD; loses essential feature": "escapes, non-functional",
                     "escapes NMD; keeps all essential features": "escapes, may function",
                     "escapes NMD; features not assessed": "escapes NMD"}[o]
            near_end = e >= L - 1 and (e - s) < 0.2 * L     # short run at the C-terminus
            a2.text(e if near_end else (s + e) / 2, 0.92, short,
                    ha="right" if near_end else "center", va="bottom",
                    fontsize=7, color=INK, transform=a2.get_xaxis_transform())
    # exon boundaries in codon units
    for b in geo["bounds"][:-1]:
        cod = (b - geo["cds_tx0"] + 1) / 3 + 1
        if 1 < cod < L:
            a2.plot([cod, cod], [0.25, 0.35], color=MUTED, lw=0.6,
                     transform=a2.get_xaxis_transform())
    if patient_ptc:
        import structure_context as sc
        three = sc.THREE.get(seq[patient_ptc - 1], "Xaa").title()
        a2.annotate(f"p.{three}{patient_ptc}* (patient)", (patient_ptc, 0.35),
                    xycoords=("data", "axes fraction"), xytext=(patient_ptc, 0.02),
                    textcoords=("data", "axes fraction"), ha="center", fontsize=8,
                    color=INK, fontweight="bold",
                    arrowprops=dict(arrowstyle="-", color=INK, lw=0.8))
    a2.set_yticks([])
    a2.set_ylim(0, 1)
    a2.set_xlabel(f"{gene} residue / codon", fontsize=8.5, color=INK)
    present = {o for _, _, o in runs}
    entries = [("null (NMD)", C_NULL, "null: NMD"),
               ("escapes NMD; loses essential feature", C_NONFUNC,
                "escapes NMD, loses an essential feature"),
               ("escapes NMD; keeps all essential features", C_FUNC,
                "escapes NMD, keeps all essential features"),
               ("escapes NMD; features not assessed", C_FUNC, "escapes NMD")]
    a2.legend(handles=[Patch(color=c, label=l) for o, c, l in entries if o in present],
              loc="upper left", bbox_to_anchor=(0, -0.42), ncol=3, frameon=False,
              fontsize=7.5, labelcolor=INK)

    for ax in (a1, a2):
        ax.set_xlim(0, L + 1)
        ax.tick_params(labelsize=7.5, colors=INK2, length=3)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(GRID)
    a2.spines["left"].set_visible(False)

    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_pdf.rsplit(".", 1)[0] + ".png", dpi=200, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------- main
def write_tsv(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(f"{v:.3f}" if isinstance(v, float) and not np.isnan(v)
                               else ("" if isinstance(v, float) else str(v))
                               for v in r) + "\n")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", required=True, help="resolved case (findings/track2/case.yaml) or gene file")
    p.add_argument("--fasta", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--pdb", help="AlphaFold model for per-residue pLDDT / RSA")
    p.add_argument("--outdir", required=True)
    p.add_argument("--model", default="esm2_t33_650M_UR50D")
    p.add_argument("--threads", type=int, default=16)
    p.add_argument("--esm-cache", help="default: <outdir>/esm_scan.npz")
    a = p.parse_args()

    cfg = caseconfig.load(a.config)
    gene, acc = cfg["gene"], cfg["uniprot"]
    calib = caseconfig.calibration(cfg)
    patient_mis = cfg["patient"].get("missense")
    patient_ptc = cfg["patient"].get("ptc_codon")

    os.makedirs(a.outdir, exist_ok=True)
    cache = a.esm_cache or os.path.join(a.outdir, "esm_scan.npz")
    seq = evs.read_fasta(a.fasta)
    header = open(a.fasta).readline()
    if f"GN={gene} " not in header + " " and acc not in header:
        sys.exit(f"{a.fasta} does not look like {gene} / {acc}: {header.strip()}")

    curated = bool(cfg["domains"])
    if curated:
        domains = nmd.load_domains(cfg["domains"])
    else:
        ft = os.path.join(a.outdir, "uniprot_features.tsv")
        if os.path.exists(ft):
            domains = nmd.load_domains(ft)
        else:
            try:
                domains = uniprot_features(acc, ft)
            except OSError as e:
                print(f"WARNING: UniProt not reachable ({e}); atlas without domain "
                      f"annotation", file=sys.stderr)
                domains = []
        print(f"{len(domains)} UniProt domain annotations (essentiality not assessed)",
              file=sys.stderr)

    for v, _ in ([(patient_mis, "")] if patient_mis else []) + calib:
        wt, pos, _ = evs.parse(v)
        if seq[pos - 1] != wt:
            sys.exit(f"{v}: sequence has {seq[pos-1]} at {pos}. Stopping.")

    geo = transcript_geometry(a.gtf, gene)
    if geo["n_aa"] != len(seq):
        sys.exit(f"transcript {geo['tid']} encodes {geo['n_aa']} aa but the "
                 f"protein sequence has {len(seq)}: wrong isoform. Stopping.")
    print(f"{gene}: transcript {geo['tid']} ({geo['how']}), {len(geo['bounds'])} exons, "
          f"{geo['n_aa']} aa; matches protein", file=sys.stderr)

    lp = esm_logprobs(seq, a.model, a.threads, cache)

    def llr(v):
        wt, pos, mt = evs.parse(v)
        return float(lp[pos - 1, AA.index(mt)] - lp[pos - 1, AA.index(wt)])
    d = [llr(v) for v, c in calib if c == "destabilising"]
    s = [llr(v) for v, c in calib if c == "stable"]
    thresholds = (max(d), min(s)) if d and s else None
    if thresholds:
        print(f"calibrated bins: like_destabilisers <= {thresholds[0]:.2f}; "
              f"like_tolerated >= {thresholds[1]:.2f}", file=sys.stderr)
    else:
        print("no calibration alleles: missense reported as uncalibrated constraint",
              file=sys.stderr)

    plddt, rsa = per_residue_structure(a.pdb, len(seq))

    mis = missense_rows(seq, lp, domains, plddt, rsa, thresholds,
                        set(cfg["stability_features"]))
    write_tsv(os.path.join(a.outdir, "atlas_missense.tsv"),
              ["pos", "wt", "mt", "esm_llr", "pctile_in_protein", "feature",
               "feature_essential", "plddt", "rsa", "esm_class", "therapy_direction"], mis)

    tr = truncating_rows(seq, geo, domains, curated)
    write_tsv(os.path.join(a.outdir, "atlas_truncating.tsv"),
              ["codon", "wt", "cds_pos", "tx_pos", "exon", "nt_to_last_junction",
               "nmd", "essential_features_lost", "outcome", "caveat"], tr)

    plot(seq, lp, domains, tr, calib, patient_mis, patient_ptc, geo,
         os.path.join(a.outdir, "atlas_figure.pdf"), gene=gene, curated=curated)

    # ---- summary ----------------------------------------------------------
    from collections import Counter
    print(f"\nMissense (all 19 x {len(seq)}):", file=sys.stderr)
    for k, n in Counter(r[10] for r in mis).most_common():
        print(f"  {n:6d}  {k}", file=sys.stderr)
    print("Premature stops by codon:", file=sys.stderr)
    for k, n in Counter(r[8] for r in tr).most_common():
        span = [r[0] for r in tr if r[8] == k]
        print(f"  {n:6d}  {k}   (codons {min(span)}-{max(span)})", file=sys.stderr)
    if patient_mis:
        pm = next(r for r in mis if f"{r[1]}{r[0]}{r[2]}" == patient_mis)
        print(f"\nPatient {patient_mis}: LLR {pm[3]:.2f}, {pm[4]:.1f} pctile, "
              f"{pm[5]}, {pm[9]} -> {pm[10]}", file=sys.stderr)
    if patient_ptc:
        pt = tr[patient_ptc - 1]
        print(f"Patient stop at codon {patient_ptc}: exon {pt[4]}, {pt[5]} nt to "
              f"last junction -> {pt[8]}", file=sys.stderr)
    print(f"wrote {a.outdir}/atlas_missense.tsv, atlas_truncating.tsv, "
          f"atlas_figure.pdf/.png", file=sys.stderr)


if __name__ == "__main__":
    main()
