#!/usr/bin/env python3
"""
Genes that cause other forms of the same disease as GENE, from the HPO
project's gene-disease annotations (OMIM and Orphanet).

Two routes, both reported per gene:
  series    OMIM diseases numbered as one series by name, e.g.
            "Mosaic variegated aneuploidy syndrome 1 ... 7"
  orphanet  Orphanet disorders that group several genes, e.g. ORPHA:1052
            "Mosaic variegated aneuploidy syndrome"

Orphanet groups larger than --max-group genes are skipped, and genes are
ranked by the number of disease links they share with GENE; the top
--max-genes are selected (broad cancer genes would otherwise pull in dozens).
The table lists every candidate with a `selected` column.

Usage
  related_genes.py --gene BUB1B --g2d genes_to_disease.txt --hpoa phenotype.hpoa \\
      --out findings/track2/related_genes.tsv
  -> prints the related gene symbols, space-separated
"""
import argparse
import re
import sys
from collections import defaultdict

SERIES = re.compile(r"^(.*?)\s+(\d+)\b")      # "<name> <number>[ ...]"


def read_g2d(path):
    gene2dis, dis2gene = defaultdict(set), defaultdict(set)
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4:
                gene2dis[f[1]].add(f[3])
                dis2gene[f[3]].add(f[1])
    return gene2dis, dis2gene


def read_names(path):
    names = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith(("#", "database_id")):
                continue
            f = line.split("\t", 2)
            if len(f) >= 2:
                names.setdefault(f[0], f[1])
    return names


def related(gene, gene2dis, dis2gene, names, max_group):
    own = gene2dis.get(gene, set())
    hits = defaultdict(list)                   # gene -> [(route, disease id, name)]
    stems = {m.group(1).lower() for d in own if d.startswith("OMIM:")
             for m in [SERIES.match(names.get(d, ""))] if m}
    for d, name in names.items():
        m = SERIES.match(name)
        if d.startswith("OMIM:") and m and m.group(1).lower() in stems and d not in own:
            for g in dis2gene.get(d, ()):
                hits[g].append(("series", d, name))
    for d in own:
        if d.startswith("ORPHA:") and 1 < len(dis2gene[d]) <= max_group:
            for g in dis2gene[d]:
                hits[g].append(("orphanet", d, names.get(d, "")))
    hits.pop(gene, None)
    return own, stems, hits


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--gene", required=True)
    a.add_argument("--g2d", required=True, help="HPO genes_to_disease.txt")
    a.add_argument("--hpoa", required=True, help="HPO phenotype.hpoa (disease names)")
    a.add_argument("--max-group", type=int, default=15)
    a.add_argument("--max-genes", type=int, default=8)
    a.add_argument("--out")
    args = a.parse_args()

    gene2dis, dis2gene = read_g2d(args.g2d)
    names = read_names(args.hpoa)
    if args.gene not in gene2dis:
        print(f"{args.gene}: no disease association in {args.g2d}", file=sys.stderr)
    own, stems, hits = related(args.gene, gene2dis, dis2gene, names, args.max_group)

    print(f"{args.gene} diseases: " + "; ".join(f"{d} {names.get(d, '')}" for d in sorted(own)),
          file=sys.stderr)
    if stems:
        print(f"series: {', '.join(sorted(stems))}", file=sys.stderr)
    ranked = sorted(hits, key=lambda g: (-len(hits[g]), g))
    chosen = set(ranked[:args.max_genes])
    if args.out:
        with open(args.out, "w") as fh:
            fh.write(f"# genes related to {args.gene} by disease (HPO annotations)\n")
            fh.write("gene\tshared_links\tselected\troute\tdisease_id\tdisease_name\n")
            for g in ranked:
                for route, d, name in hits[g]:
                    fh.write(f"{g}\t{len(hits[g])}\t{'yes' if g in chosen else 'no'}\t"
                             f"{route}\t{d}\t{name}\n")
    for g in ranked:
        print(f"  {'*' if g in chosen else ' '} {g:<10} " +
              "; ".join(f"{r} {d} {n}" for r, d, n in hits[g]), file=sys.stderr)
    if len(ranked) > args.max_genes:
        print(f"{len(ranked)} related genes; selected the {args.max_genes} with most shared "
              f"links (--max-genes)", file=sys.stderr)
    print(" ".join(sorted(chosen)))


if __name__ == "__main__":
    main()
