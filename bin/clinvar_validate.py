#!/usr/bin/env python3
"""
Independent validation of the BUB1B variant atlas against ClinVar.

The atlas's missense bins are set by six experimentally characterised MVA
alleles, so those six cannot also test it. ClinVar provides variants the
atlas never saw. For every BUB1B missense SNV in ClinVar this script:

  1. maps the genomic change to its codon and amino-acid change from the
     reference genome and the Ensembl transcript (and checks the reference
     residue against the UniProt sequence -- a failed check stops the run);
  2. looks the change up in the atlas (ESM-2 score, calibrated bin, therapy
     direction);
  3. compares pathogenic/likely pathogenic with benign/likely benign variants:
     cross-tabulation by atlas bin and a rank-based AUROC of the ESM score.

It also triages ClinVar's variants of uncertain significance (VUS) and
conflicting classifications by atlas bin -- the use case the atlas is for.

Caveats reported with the result:
  * ClinVar assertions sometimes use in-silico predictors as supporting
    evidence (ACMG PP3/BP4), so agreement is not fully independent.
  * BUB1B benign variants are mostly population polymorphisms and pathogenic
    ones mostly MVA alleles: the test asks whether the atlas separates those
    two groups, not whether it predicts disease in general.

Usage
  clinvar_validate.py --clinvar /opt/refs/clinvar/clinvar.chr.vcf.gz \\
      --fasta /opt/refs/GRCh38.fa --gtf /opt/refs/Homo_sapiens.GRCh38.116.gtf.gz \\
      --protein resources/O60566.fasta \\
      --atlas findings/track2/atlas/atlas_missense.tsv \\
      --out findings/track2/clinvar_validation.tsv
"""

import argparse
import csv
import os
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import nmd_check as nmd            # noqa: E402
import esm_variant_score as evs    # noqa: E402

CODON = {}
_b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a in enumerate(_b):
    for j, b in enumerate(_b):
        for k, c in enumerate(_b):
            CODON[a + b + c] = _aa[16 * i + 4 * j + k]
COMP = str.maketrans("ACGT", "TGCA")

PLP = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}
BLB = {"Benign", "Likely_benign", "Benign/Likely_benign"}
UNC = {"Uncertain_significance", "Conflicting_classifications_of_pathogenicity",
       "Conflicting_interpretations_of_pathogenicity"}


def cds_map(gtf, gene):
    """Ordered list of genomic CDS positions (transcription order) and strand."""
    tx = nmd.scan(gtf, gene=gene)
    tid, rec, how = nmd.pick(tx)
    strand = rec["strand"]
    pos = []
    for s, e in nmd.ordered(rec["cds"], strand):
        rng = range(s, e + 1) if strand == "+" else range(e, s - 1, -1)
        pos.extend(rng)
    return tid, rec["chrom"], strand, pos


def clnsig_class(sig):
    parts = set(sig.replace(",", "|").split("|"))
    if parts & PLP and not parts & BLB:
        return "P/LP"
    if parts & BLB and not parts & PLP:
        return "B/LB"
    if parts & UNC or (parts & PLP and parts & BLB):
        return "VUS/conflicting"
    return "other"


def auroc(pos_scores, neg_scores):
    """P(score of a pathogenic variant < score of a benign one): lower LLR =
    more damaging. Ties count half. Returns None if either group is empty."""
    if not pos_scores or not neg_scores:
        return None
    wins = 0.0
    for p in pos_scores:
        for n in neg_scores:
            wins += 1.0 if p < n else (0.5 if p == n else 0.0)
    return wins / (len(pos_scores) * len(neg_scores))


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--clinvar", required=True)
    a.add_argument("--fasta", required=True)
    a.add_argument("--gtf", required=True)
    a.add_argument("--gene", default="BUB1B")
    a.add_argument("--protein", required=True)
    a.add_argument("--atlas", required=True)
    a.add_argument("--out", required=True)
    args = a.parse_args()

    import pysam

    protein = evs.read_fasta(args.protein)
    tid, chrom, strand, cds = cds_map(args.gtf, args.gene)
    chrom_vcf = chrom if chrom.startswith("chr") else f"chr{chrom}"
    idx = {g: i for i, g in enumerate(cds)}          # genomic pos -> CDS index (0-based)
    fa = pysam.FastaFile(args.fasta)
    ref_cds = "".join(fa.fetch(chrom_vcf, p - 1, p).upper() for p in cds)
    if strand == "-":
        ref_cds = ref_cds.translate(COMP)
    translated = "".join(CODON.get(ref_cds[i:i + 3], "X") for i in range(0, len(ref_cds) - 2, 3))
    if translated != protein[:len(translated)]:
        mism = sum(x != y for x, y in zip(translated, protein))
        sys.exit(f"reference CDS of {tid} does not translate to the UniProt sequence "
                 f"({mism} mismatches): wrong transcript or assembly. Stopping.")
    print(f"{tid}: CDS translates to the UniProt sequence ({len(translated)} aa)",
          file=sys.stderr)

    atlas = {}
    with open(args.atlas) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            atlas[(int(r["pos"]), r["wt"], r["mt"])] = r

    calib = {v for v, _ in evs.CALIBRATION}
    lo, hi = min(cds), max(cds)
    vcf = pysam.VariantFile(args.clinvar)
    rows = []
    for rec in vcf.fetch(chrom_vcf, lo - 1, hi):
        if rec.pos not in idx or len(rec.ref) != 1 or not rec.alts:
            continue
        gene_info = str(rec.info.get("GENEINFO", ""))
        if f"{args.gene}:" not in gene_info:
            continue
        sig = rec.info.get("CLNSIG")
        sig = "|".join(sig) if isinstance(sig, tuple) else str(sig or "")
        rev = rec.info.get("CLNREVSTAT")
        rev = ",".join(rev) if isinstance(rev, tuple) else str(rev or "")
        for alt in rec.alts:
            if len(alt) != 1 or alt == rec.ref:
                continue
            i = idx[rec.pos]
            c0 = i - i % 3
            codon = list(ref_cds[c0:c0 + 3])
            b = alt.upper()
            codon[i % 3] = b.translate(COMP) if strand == "-" else b
            wt = CODON[ref_cds[c0:c0 + 3]]
            mt = CODON.get("".join(codon), "X")
            aa_pos = c0 // 3 + 1
            if mt == wt or mt in "*X":
                continue                      # synonymous or nonsense: not in the missense atlas
            v = f"{wt}{aa_pos}{mt}"
            at = atlas.get((aa_pos, wt, mt))
            rows.append({
                "variant": v, "chrom": rec.chrom, "pos": rec.pos, "ref": rec.ref,
                "alt": alt, "clinvar_id": rec.id or "", "clnsig": sig,
                "review_status": rev, "class": clnsig_class(sig),
                "calibration_allele": v in calib,
                "esm_llr": float(at["esm_llr"]) if at else None,
                "atlas_bin": at["esm_class"] if at else "",
                "therapy_direction": at["therapy_direction"] if at else "",
                "feature": at["feature"] if at else "",
            })

    cols = ["variant", "chrom", "pos", "ref", "alt", "clinvar_id", "clnsig",
            "review_status", "class", "calibration_allele", "esm_llr", "atlas_bin",
            "therapy_direction", "feature"]
    with open(args.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda r: (r["class"], r["esm_llr"] if r["esm_llr"] is not None else 99)):
            fh.write("\t".join("" if r[c] is None else
                               (f"{r[c]:.3f}" if isinstance(r[c], float) else str(r[c]))
                               for c in cols) + "\n")

    # ---- summary -----------------------------------------------------------
    test = [r for r in rows if not r["calibration_allele"] and r["esm_llr"] is not None]
    print(f"\nClinVar {args.gene} missense SNVs: {len(rows)} "
          f"({sum(r['calibration_allele'] for r in rows)} calibration alleles excluded from the test)",
          file=sys.stderr)
    by = Counter(r["class"] for r in test)
    print("  " + ", ".join(f"{k}: {n}" for k, n in sorted(by.items())), file=sys.stderr)

    bins = ["like_destabilisers", "indeterminate", "like_tolerated"]
    print(f"\n{'atlas bin':<20}{'P/LP':>7}{'B/LB':>7}{'VUS/confl':>11}", file=sys.stderr)
    for b in bins:
        c = Counter(r["class"] for r in test if r["atlas_bin"] == b)
        print(f"{b:<20}{c['P/LP']:>7}{c['B/LB']:>7}{c['VUS/conflicting']:>11}", file=sys.stderr)

    p = [r["esm_llr"] for r in test if r["class"] == "P/LP"]
    n = [r["esm_llr"] for r in test if r["class"] == "B/LB"]
    auc = auroc(p, n)
    if auc is None:
        print("\nAUROC: not computable (need both P/LP and B/LB variants)", file=sys.stderr)
    else:
        print(f"\nAUROC of ESM score, P/LP vs B/LB: {auc:.2f} "
              f"(n = {len(p)} vs {len(n)}; 0.5 = no separation)", file=sys.stderr)
        if min(len(p), len(n)) < 10:
            print("  small groups: report as indicative, with the counts", file=sys.stderr)

    vus = Counter(r["therapy_direction"] for r in test if r["class"] == "VUS/conflicting")
    if vus:
        print("\nVUS / conflicting variants by atlas therapy direction:", file=sys.stderr)
        for k, c in vus.most_common():
            print(f"  {c:4d}  {k}", file=sys.stderr)
    print(f"\nwrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
