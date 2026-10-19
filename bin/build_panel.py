#!/usr/bin/env python3
"""
Build the candidate-gene panel BED from an Ensembl GTF.

Reads a gene-symbol list (one per line, '#' comments and '##' section headers
ignored), extracts each gene's span from the GTF, pads it, merges overlaps and
writes a sorted, chr-prefixed BED.

Genes absent from the GTF are reported on stderr and written to a .missing file
rather than silently dropped -- a renamed symbol must be fixed, not ignored.

Usage:
  build_panel.py --genes docs/panel_genes.txt \\
                 --gtf  /opt/refs/Homo_sapiens.GRCh38.116.gtf.gz \\
                 --pad  10000 \\
                 --out  resources/mva_panel.bed
"""

import argparse
import gzip
import re
import sys

PRIMARY = {str(c) for c in range(1, 23)} | {"X", "Y", "MT"}
NAME_RE = re.compile(r'gene_name "([^"]+)"')
BIOTYPE_RE = re.compile(r'gene_biotype "([^"]+)"')


def load_genes(path):
    genes = []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            genes.append(s)
    if len(genes) != len(set(genes)):
        dupes = {g for g in genes if genes.count(g) > 1}
        print(f"note: duplicate symbols collapsed: {sorted(dupes)}", file=sys.stderr)
    return sorted(set(genes))


def chr_name(c):
    return "chrM" if c in ("MT", "M") else f"chr{c}"


def scan_gtf(path, wanted):
    """Return {symbol: (chrom, start, end)} using the widest span per symbol."""
    found = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if f[2] != "gene" or f[0] not in PRIMARY:
                continue
            m = NAME_RE.search(f[8])
            if not m or m.group(1) not in wanted:
                continue
            bt = BIOTYPE_RE.search(f[8])
            if bt and bt.group(1) not in ("protein_coding", "protein_coding_LoF"):
                continue
            sym = m.group(1)
            start, end = int(f[3]) - 1, int(f[4])   # GTF is 1-based inclusive
            if sym in found:
                c, s0, e0 = found[sym]
                if c != f[0]:
                    print(f"warning: {sym} on both {c} and {f[0]}; keeping {c}",
                          file=sys.stderr)
                    continue
                found[sym] = (c, min(s0, start), max(e0, end))
            else:
                found[sym] = (f[0], start, end)
    return found


def merge(regions):
    """regions: list of (chrom, start, end, name) -> merged, sorted BED rows."""
    def key(r):
        c = r[0][3:]
        order = {"X": 23, "Y": 24, "M": 25}.get(c, None)
        return (order if order else int(c), r[1])

    regions.sort(key=key)
    out = []
    for chrom, start, end, name in regions:
        if out and out[-1][0] == chrom and start <= out[-1][2]:
            prev = out[-1]
            names = prev[3].split(",")
            if name not in names:
                names.append(name)
            out[-1] = (chrom, prev[1], max(prev[2], end), ",".join(sorted(names)))
        else:
            out.append((chrom, start, end, name))
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--genes", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--pad", type=int, default=10000,
                   help="bp added either side of each gene span (default 10000)")
    p.add_argument("--out", required=True)
    a = p.parse_args()

    wanted = load_genes(a.genes)
    found = scan_gtf(a.gtf, set(wanted))

    missing = [g for g in wanted if g not in found]
    if missing:
        with open(a.out + ".missing", "w") as fh:
            fh.write("\n".join(missing) + "\n")
        print(f"MISSING {len(missing)} symbol(s), written to {a.out}.missing:",
              file=sys.stderr)
        print("  " + " ".join(missing), file=sys.stderr)
        print("  -> check for HGNC renaming before using this panel.", file=sys.stderr)

    regions = [(chr_name(c), max(0, s - a.pad), e + a.pad, sym)
               for sym, (c, s, e) in found.items()]
    merged = merge(regions)

    with open(a.out, "w") as fh:
        for chrom, start, end, name in merged:
            fh.write(f"{chrom}\t{start}\t{end}\t{name}\n")

    bp = sum(e - s for _, s, e, _ in merged)
    print(f"{len(found)}/{len(wanted)} genes -> {len(merged)} regions, "
          f"{bp/1e6:.2f} Mb (pad {a.pad} bp) -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
