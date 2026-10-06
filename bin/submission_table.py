#!/usr/bin/env python3
"""
Turn a Track 1 submission CSV into a Markdown results table for the report,
joining each submitted variant with its Exomiser annotation.

Default output is gene-level (no genomic coordinates), suitable for sharing
and for a public repo. Add --coords to include positions and alleles.

Usage:
  submission_table.py --submission results/submission/submission.csv \
      --variants results/ranking/proband_vendor.variants.tsv [--coords]
"""

import argparse
import csv


def chr_name(c):
    c = str(c)
    if c.startswith("chr"):
        return c
    return "chrM" if c in ("MT", "M") else f"chr{c}"


def load_exomiser(path):
    """Index Exomiser variant rows by (chrom, pos, ref, alt); keep all MOIs."""
    idx = {}
    with open(path, newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        r.fieldnames = [f.lstrip("#") for f in r.fieldnames]
        for row in r:
            key = (chr_name(row["CONTIG"]), int(row["START"]), row["REF"], row["ALT"])
            e = idx.setdefault(key, {"row": row, "mois": set()})
            e["mois"].add(row.get("MOI", "?"))
    return idx


def describe(key, idx, coords):
    e = idx.get(key)
    if e is None:
        return {"gene": "?", "cons": "not in Exomiser output", "gt": "?",
                "freq": "?", "path": "?", "clinvar": "?", "mois": "?",
                "loc": f"{key[0]}:{key[1]} {key[2]}>{key[3]}" if coords else ""}
    row = e["row"]
    freq = row.get("MAX_FREQ", "")
    try:
        freq = f"{float(freq):.4f}%" if freq not in ("", ".") else "absent"
    except ValueError:
        pass
    path = row.get("MAX_PATH", "")
    src = row.get("MAX_PATH_SOURCE", "")
    return {
        "gene": row.get("GENE_SYMBOL", "?"),
        "cons": row.get("FUNCTIONAL_CLASS", "?").lower().replace("_", " "),
        "gt": row.get("GENOTYPE", "?"),
        "freq": freq,
        "path": f"{path} ({src})" if path else "—",
        "clinvar": row.get("CLINVAR_PRIMARY_INTERPRETATION", "") or "—",
        "mois": "/".join(sorted(e["mois"])),
        "loc": f"{key[0]}:{key[1]} {key[2]}>{key[3]}" if coords else "",
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--submission", required=True)
    p.add_argument("--variants", required=True, help="Exomiser *.variants.tsv")
    p.add_argument("--coords", action="store_true", help="include genomic positions")
    a = p.parse_args()

    idx = load_exomiser(a.variants)

    cols = ["Rank", "Row", "Gene", "EPCR", "Exomiser MOI", "Consequence", "Genotype",
            "Max pop. freq", "Top pathogenicity", "ClinVar"]
    if a.coords:
        cols.insert(3, "Variant")
    print("| " + " | ".join(cols) + " |")
    print("|" + "---|" * len(cols))

    with open(a.submission, newline="") as fh:
        for i, r in enumerate(csv.DictReader(fh), 1):
            keys = [(r["chrom_1"], int(r["pos_1"]), r["ref_1"], r["alt_1"])]
            if r["chrom_2"]:
                keys.append((r["chrom_2"], int(r["pos_2"]), r["ref_2"], r["alt_2"]))
            ds = [describe(k, idx, a.coords) for k in keys]
            join = lambda f: " + ".join(d[f] for d in ds)   # noqa: E731
            cells = [str(i), "pair" if len(ds) == 2 else "single", ds[0]["gene"],
                     r["epcr"], join("mois"), join("cons"), join("gt"),
                     join("freq"), join("path"), join("clinvar")]
            if a.coords:
                cells.insert(3, join("loc"))
            print("| " + " | ".join(cells) + " |")


if __name__ == "__main__":
    main()
