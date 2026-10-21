#!/usr/bin/env python3
"""
AlphaMissense as a second, independent missense predictor next to ESM-2.

AlphaMissense (Cheng et al., Science 2023) was trained on primate population
frequency and structural context, not on the alleles used here to calibrate
ESM-2, so it is a test of the ESM-2 calls rather than a repeat of them.

  1. Calibration: AlphaMissense scores of the experimentally classed alleles
     (gene file) and the patient's missense; do the classes separate?
  2. Agreement with ESM-2 over every missense change AlphaMissense scores
     (Spearman), overall and per curated feature.
  3. ClinVar (if clinvar_validation.tsv exists): AUROC of AlphaMissense for
     P/LP vs B/LB, next to ESM-2's on the same variants; agreement of the two
     predictors on variants of uncertain significance.

Scores are taken for the gene's UniProt accession from the genome-wide file
(region from the gene's CDS in the GTF, tabix). Higher AlphaMissense score =
more likely pathogenic; lower ESM-2 LLR = less tolerated.

Usage
  alphamissense_compare.py --config findings/track2/case.yaml \\
      --am /opt/refs/alphamissense/AlphaMissense_hg38.tsv.gz --gtf GTF \\
      --atlas findings/track2/atlas/atlas_missense.tsv \\
      [--clinvar findings/track2/clinvar_validation.tsv] --outdir findings/track2/alphamissense
"""
import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import caseconfig              # noqa: E402
import nmd_check as nmd        # noqa: E402

COLS = ["#CHROM", "POS", "REF", "ALT", "genome", "uniprot_id", "transcript_id",
        "protein_variant", "am_pathogenicity", "am_class"]


def gene_span(gtf, gene):
    tx = nmd.scan(gtf, gene=gene)
    _, rec, _ = nmd.pick(tx)
    starts = [s for s, _ in rec["cds"]]
    ends = [e for _, e in rec["cds"]]
    chrom = rec["chrom"]
    return (chrom if chrom.startswith("chr") else f"chr{chrom}"), min(starts), max(ends)


def load_am(path, uniprot, region):
    """{(pos, wt, mt): (score, class)} for one protein."""
    import pysam
    tbx = pysam.TabixFile(path)
    header = [h for h in tbx.header if h.startswith("#CHROM")]
    if not header or header[-1].split("\t")[:len(COLS)] != COLS:
        sys.exit(f"{path}: unexpected AlphaMissense columns {header[-1:] or '(none)'}")
    chrom, lo, hi = region
    if chrom not in tbx.contigs:
        chrom = chrom.removeprefix("chr")
    out = {}
    for line in tbx.fetch(chrom, lo - 1, hi):
        f = line.split("\t")
        if f[5].split("-")[0] != uniprot:
            continue
        pv = f[7]
        try:
            key = (int(pv[1:-1]), pv[0], pv[-1])
        except ValueError:
            continue
        out[key] = (float(f[8]), f[9])
    return out


def ranks(x):
    from scipy.stats import rankdata
    return rankdata(x)


def auroc(pos, neg):
    """P(score of a positive > score of a negative); ties count half."""
    if not pos or not neg:
        return None
    from scipy.stats import mannwhitneyu
    u = mannwhitneyu(pos, neg, alternative="two-sided").statistic
    return u / (len(pos) * len(neg))


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--config", required=True)
    a.add_argument("--am", required=True)
    a.add_argument("--gtf", required=True)
    a.add_argument("--atlas", required=True)
    a.add_argument("--clinvar")
    a.add_argument("--outdir", required=True)
    args = a.parse_args()
    from scipy.stats import spearmanr

    cfg = caseconfig.load(args.config)
    gene, acc = cfg["gene"], cfg["uniprot"]
    am = load_am(args.am, acc, gene_span(args.gtf, gene))
    if not am:
        sys.exit(f"no AlphaMissense scores for {gene} / {acc}")
    os.makedirs(args.outdir, exist_ok=True)
    out = []

    def say(s=""):
        print(s, file=sys.stderr)
        out.append(s)

    # ---- atlas join ------------------------------------------------------
    atlas = {}
    with open(args.atlas) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            atlas[(int(r["pos"]), r["wt"], r["mt"])] = r
    joined = [(k, am[k], atlas[k]) for k in sorted(am) if k in atlas]
    # AlphaMissense numbers residues on its canonical transcript; a different
    # isoform would silently mis-join, so check the reference residues.
    wt_at = {p: w for (p, w, _) in atlas}
    mism = sum(1 for (p, w, _) in am if p in wt_at and wt_at[p] != w)
    if mism > 0.02 * len(am):
        sys.exit(f"{mism} of {len(am)} AlphaMissense variants have a different reference "
                 f"residue from the protein sequence: isoform mismatch. Stopping.")
    with open(os.path.join(args.outdir, "am_vs_esm.tsv"), "w") as fh:
        fh.write("pos\twt\tmt\tam_pathogenicity\tam_class\tesm_llr\tesm_class\tfeature\n")
        for (p, w, m), (s, c), r in joined:
            fh.write(f"{p}\t{w}\t{m}\t{s:.4f}\t{c}\t{r['esm_llr']}\t{r['esm_class']}\t{r['feature']}\n")
    say(f"{gene} ({acc}): {len(am)} AlphaMissense missense scores, {len(joined)} matched to the atlas")

    # ---- 1. calibration --------------------------------------------------
    calib = caseconfig.calibration(cfg)
    patient = cfg.get("patient", {}).get("missense")
    rows = [(v, cls) for v, cls in calib] + ([(patient, "PATIENT")] if patient else [])
    if rows:
        say("\n1. Calibration alleles")
        say(f"   {'variant':<9}{'class':<15}{'AlphaMissense':>14}  {'AM class':<18}{'ESM LLR':>8}")
        for v, cls in rows:
            k = (int(v[1:-1]), v[0], v[-1])
            s, c = am.get(k, (None, "not scored"))
            e = atlas.get(k, {}).get("esm_llr", "")
            say(f"   {v:<9}{cls:<15}{'' if s is None else f'{s:.3f}':>14}  {c:<18}{e:>8}")
        d = [am[(int(v[1:-1]), v[0], v[-1])][0] for v, c in calib
             if c == "destabilising" and (int(v[1:-1]), v[0], v[-1]) in am]
        st = [am[(int(v[1:-1]), v[0], v[-1])][0] for v, c in calib
              if c == "stable" and (int(v[1:-1]), v[0], v[-1]) in am]
        if d and st:
            sep = min(d) > max(st)
            say(f"   destabilising min {min(d):.3f} vs stable max {max(st):.3f}: classes "
                f"{'SEPARATE' if sep else 'OVERLAP'}")

    # ---- 2. agreement with ESM-2 -----------------------------------------
    say("\n2. Agreement with ESM-2 (Spearman; negative = agree, since the scales run opposite)")
    s_am = [j[1][0] for j in joined]
    s_esm = [float(j[2]["esm_llr"]) for j in joined]
    rho = spearmanr(s_am, s_esm).statistic
    say(f"   all matched missense (n = {len(joined)}): rho = {rho:+.2f}")
    feats = {}
    for _, (s, _), r in joined:
        feats.setdefault(r["feature"], []).append((s, float(r["esm_llr"])))
    for f, xs in sorted(feats.items(), key=lambda kv: -len(kv[1])):
        if len(xs) >= 100:
            say(f"   {f:<24} n = {len(xs):5d}  rho = "
                f"{spearmanr([x for x, _ in xs], [y for _, y in xs]).statistic:+.2f}")

    # ---- 3. ClinVar --------------------------------------------------------
    if args.clinvar and os.path.exists(args.clinvar):
        say("\n3. ClinVar (calibration alleles excluded)")
        cv = []
        with open(args.clinvar) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["calibration_allele"] == "True" or not r["esm_llr"]:
                    continue
                v = r["variant"]
                k = (int(v[1:-1]), v[0], v[-1])
                if k in am:
                    cv.append((r["class"], am[k][0], am[k][1], float(r["esm_llr"]), r["atlas_bin"]))
        p = [x for x in cv if x[0] == "P/LP"]
        b = [x for x in cv if x[0] == "B/LB"]
        a_am = auroc([x[1] for x in p], [x[1] for x in b])
        a_esm = auroc([-x[3] for x in p], [-x[3] for x in b])
        say(f"   P/LP n = {len(p)}, B/LB n = {len(b)}")
        if a_am is not None:
            say(f"   AUROC AlphaMissense {a_am:.2f}, ESM-2 {a_esm:.2f} (same variants)")
            if min(len(p), len(b)) < 10:
                say("   small groups: indicative only")
        for name, grp in (("B/LB", b), ("P/LP", p)):
            if grp:
                lb = sum(x[2] == "likely_benign" for x in grp)
                lp = sum(x[2] == "likely_pathogenic" for x in grp)
                say(f"   {name}: AlphaMissense likely_benign {lb}, likely_pathogenic {lp}, "
                    f"ambiguous {len(grp) - lb - lp}")
        vus = [x for x in cv if x[0] == "VUS/conflicting"]
        if vus:
            both_bad = sum(x[2] == "likely_pathogenic" and x[4] == "like_destabilisers" for x in vus)
            both_ok = sum(x[2] == "likely_benign" and x[4] == "like_tolerated" for x in vus)
            say(f"   VUS/conflicting n = {len(vus)}: both predictors damaging {both_bad}, "
                f"both tolerated {both_ok}, rest discordant or intermediate "
                f"{len(vus) - both_bad - both_ok}")

    with open(os.path.join(args.outdir, "summary.txt"), "w") as fh:
        fh.write("\n".join(out) + "\n")
    print(f"wrote {args.outdir}/am_vs_esm.tsv, summary.txt", file=sys.stderr)


if __name__ == "__main__":
    main()
