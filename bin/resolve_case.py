#!/usr/bin/env python3
"""
Turn the case file into the resolved case that Track 2 reads.

  track2.gene / track2.variants in params.yaml:
    auto      take them from the top-ranked row of the Track 1 submission
    explicit  gene symbol; variants as a list of {chrom, pos, ref, alt}

Each SNV is placed on the gene's MANE/canonical transcript (Ensembl GTF) and
its consequence derived from the reference genome: premature stop (ptc),
missense, synonymous, or outside the coding sequence. The resolved case names
the first ptc and the first missense as the patient's alleles for the
NMD, ESM, structure and splice steps.

If the gene has no config/genes/<gene>.yaml yet, one is generated from
UniProt (make_gene_config.py), so a new gene needs no hand-written config.

Usage
  resolve_case.py --case params.yaml --gtf GTF --fasta GRCh38.fa \\
      [--submission results_vcf/submission/submission.csv] \\
      --out findings/track2/case.yaml
"""
import argparse
import csv
import os
import sys

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from clinvar_validate import CODON, COMP, cds_map   # noqa: E402
import make_gene_config                             # noqa: E402

GENES_DIR = "config/genes"


def from_submission(path, gene=None):
    """Gene and variants of the first submission row (for `gene`, if given)."""
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            g = (r.get("notes") or "").split(":", 1)[0].strip()
            if gene and g.upper() != gene.upper():
                continue
            vs = []
            for i in ("1", "2"):
                if r.get(f"pos_{i}"):
                    vs.append({"chrom": r[f"chrom_{i}"], "pos": int(r[f"pos_{i}"]),
                               "ref": r[f"ref_{i}"], "alt": r[f"alt_{i}"]})
            return g, vs
    sys.exit(f"{path}: no row{' for ' + gene if gene else ''}")


def consequence(v, chrom, strand, cds, ref_cds, idx):
    """Annotate one variant in place."""
    if len(v["ref"]) != 1 or len(v["alt"]) != 1:
        v["type"] = "indel (not classified)"
        return
    if v["pos"] not in idx:
        v["type"] = "non-coding"
        return
    i = idx[v["pos"]]
    c0 = i - i % 3
    codon = list(ref_cds[c0:c0 + 3])
    ref_base = v["ref"].upper().translate(COMP) if strand == "-" else v["ref"].upper()
    if codon[i % 3] != ref_base:
        sys.exit(f"{v['chrom']}:{v['pos']} REF {v['ref']} does not match the reference "
                 f"genome: wrong assembly or contig naming. Stopping.")
    alt = v["alt"].upper()
    codon[i % 3] = alt.translate(COMP) if strand == "-" else alt
    wt, mt = CODON[ref_cds[c0:c0 + 3]], CODON.get("".join(codon), "X")
    n = c0 // 3 + 1
    v["codon"] = n
    v["protein"] = f"{wt}{n}{mt}"
    v["type"] = "ptc" if mt == "*" else ("synonymous" if mt == wt else "missense")


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--case", required=True)
    a.add_argument("--gtf", required=True)
    a.add_argument("--fasta", required=True)
    a.add_argument("--submission")
    a.add_argument("--out", required=True)
    args = a.parse_args()

    case = yaml.safe_load(open(args.case)) or {}
    t2 = case.get("track2") or {}
    gene = t2.get("gene", "auto")
    variants = t2.get("variants", "auto")

    if gene == "auto" or variants == "auto":
        if not args.submission or not os.path.exists(args.submission):
            sys.exit("track2 gene/variants are 'auto' but there is no Track 1 "
                     f"submission ({args.submission}): run Track 1 first, or set "
                     "them explicitly in the case file.")
        g, vs = from_submission(args.submission, None if gene == "auto" else gene)
        source = f"{args.submission}, rank 1" + ("" if gene == "auto" else f" row for {g}")
        gene = g if gene == "auto" else gene
        variants = vs if variants == "auto" else variants
    else:
        source = args.case
    variants = [dict(v) for v in variants]

    gene_config = os.path.join(GENES_DIR, f"{gene.lower()}.yaml")
    if not os.path.exists(gene_config):
        print(f"no {gene_config}: generating it from UniProt", file=sys.stderr)
        make_gene_config.write(gene, gene_config)

    import pysam
    tid, chrom, strand, cds = cds_map(args.gtf, gene)
    fa = pysam.FastaFile(args.fasta)
    contig = chrom if chrom in fa.references else f"chr{chrom}"
    ref_cds = "".join(fa.fetch(contig, p - 1, p).upper() for p in cds)
    if strand == "-":
        ref_cds = ref_cds.translate(COMP)
    idx = {g: i for i, g in enumerate(cds)}
    for v in variants:
        consequence(v, chrom, strand, cds, ref_cds, idx)

    ptc = next((v for v in variants if v["type"] == "ptc"), None)
    mis = next((v for v in variants if v["type"] == "missense"), None)
    patient = {"chrom": contig}
    if ptc:
        patient.update(ptc_pos=ptc["pos"], ptc_codon=ptc["codon"],
                       ptc_allele=f"{ptc['pos']} {ptc['ref']} {ptc['alt']}")
    if mis:
        patient.update(missense=mis["protein"],
                       missense_allele=f"{mis['pos']} {mis['ref']} {mis['alt']}")

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write(f"# Resolved by bin/resolve_case.py from {source}.\n"
                 f"# Generated file: edit {args.case} instead.\n")
        yaml.safe_dump({"gene": gene, "gene_config": gene_config, "transcript": tid,
                        "variants": variants, "patient": patient},
                       fh, sort_keys=False, default_flow_style=None)

    print(f"{gene} ({tid}, {strand} strand), from {source}:", file=sys.stderr)
    for v in variants:
        print(f"  {v['chrom']}:{v['pos']} {v['ref']}>{v['alt']}  "
              f"{v.get('protein', '')}  {v['type']}", file=sys.stderr)
    if not (ptc or mis):
        print("  no premature stop or missense: patient-specific Track 2 steps "
              "will be skipped", file=sys.stderr)
    print(f"wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
