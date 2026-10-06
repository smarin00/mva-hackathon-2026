#!/usr/bin/env python3
"""
Write the Track 1 submission TSV.

The scoring rules live here, in one place, so they are never re-derived by
hand at 23:00 on 24 October. Three things the scorer does that are easy to
get wrong:

  1. A "full match" requires the row's variant SET to equal the answer key's
     variant set exactly. If the answer is a compound heterozygote and you
     submit the two alleles on separate rows, you can never score a full
     match -- only partial credit. Paired candidates therefore go on ONE row
     using the *_2 columns, with the singles kept as lower-ranked hedges.

  2. Rank points are tiered and steep (1 -> 100, 2-3 -> 50, 4-5 -> 25,
     6-10 -> 10, then nothing). Ties are broken by file order, so every EPCR
     must be strictly distinct and descending. Never tie your best candidate
     with anything.

  3. F-max is computed over variants above an EPCR threshold, and the
     finding_type column does not exclude anything from it. A secondary
     finding sitting above your primary candidate counts as a false positive
     and drags precision down. Secondary findings get the lowest EPCRs.

Output format matches the Space's scorer (evaluation.py): a comma-separated
.csv, read with csv.DictReader, requiring proband_id. Compare the header with
the Space's track1_submission_template.csv before the first real submission.
"""

import argparse
import csv
import sys
from dataclasses import dataclass, field

# Exactly the columns of the Space's Track 1 template / evaluation.py.
COLUMNS = [
    "proband_id",
    "chrom_1", "pos_1", "ref_1", "alt_1",
    "chrom_2", "pos_2", "ref_2", "alt_2",
    "epcr", "finding_type", "notes",
]
DEFAULT_PROBAND = "PROBAND01"   # the only proband_id the Space accepts

import re

MAX_ROWS = 10          # rows past 10 earn zero rank points

# Hypervariable / reference-problem loci: short-read calls here are dominated
# by mapping artefacts, so they are excluded from candidates.
ARTEFACT_GENES = re.compile(r"^(HLA-|KIR\d|MUC\d|OR\d+[A-Z]|PRAMEF|NBPF|FAM231|TRB[VJC]|TRA[VJC]|IG[HKL][VDJC])")

RECESSIVE_MOIS = {"AR", "XR"}

# Two variants closer than this are treated as one event / same haplotype
# (complex indels split into records, repeat-length alleles) and never paired.
MIN_PAIR_DISTANCE = 50

# Primary rows need at least this Exomiser gene combined score; unused slots
# stay empty rather than being filled with noise. Curated secondaries exempt.
DEFAULT_MIN_GENE_SCORE = 0.1

# Variants ClinVar classifies as benign are never submitted.
BENIGN = re.compile(r"BENIGN", re.I)


def load_curation(path):
    """Gene-level curation decisions: gene <TAB> action <TAB> note.

    action is 'exclude' (drop all rows for the gene) or 'secondary'
    (keep as a secondary/incidental finding with the given note).
    Lines starting with '#' are comments.
    """
    cur = {}
    if not path:
        return cur
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            gene, action = parts[0].strip(), parts[1].strip().lower()
            note = parts[2].strip() if len(parts) > 2 else ""
            if action not in ("exclude", "secondary"):
                sys.exit(f"curation: unknown action {action!r} for {gene}")
            cur[gene] = (action, note)
    return cur


def apply_curation(candidates, curation, min_gene_score):
    kept = []
    for c in candidates:
        action, note = curation.get(c.gene, (None, ""))
        if action == "exclude":
            continue
        if action == "secondary":
            c.finding_type = "secondary"
            c.rationale = note or c.rationale
        elif c.gene_score < min_gene_score:
            continue
        kept.append(c)
    return kept
EPCR_TOP = 0.95        # highest confidence assigned to the leading candidate
EPCR_FLOOR = 0.01      # secondary findings sit below anything primary


@dataclass
class Candidate:
    """One submission row. Two variants means a compound-het pair."""
    variants: list = field(default_factory=list)   # list of (chrom, pos, ref, alt)
    gene: str = ""
    finding_type: str = "primary"
    rationale: str = ""
    score: float = 0.0                              # internal ranking score
    mois: set = field(default_factory=set)          # Exomiser modes this variant passed
    gene_score: float = 0.0                         # Exomiser gene combined score
    homozygous: bool = False

    @property
    def is_pair(self):
        return len(self.variants) == 2


def chr_name(contig):
    """Submission uses chr-prefixed names; Exomiser may report '1' or 'MT'."""
    c = str(contig)
    if c.startswith("chr"):
        return c
    return "chrM" if c in ("MT", "M") else f"chr{c}"


def load_ranked(path):
    """Read the Exomiser 15 TSV_VARIANT table into Candidate objects.

    Exomiser writes one row per variant per mode of inheritance, so the same
    variant can appear several times. Keep the best-scoring row per variant.
    Ranking score = gene combined score (phenotype + variant), with the
    variant score as a small tiebreak; contributing variants get a bonus.

    Check the header of your first real output against these column names;
    the TSV schema has changed between Exomiser releases.
    """
    best = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        # header may start with '#RANK'
        reader.fieldnames = [f.lstrip("#") for f in reader.fieldnames]
        for row in reader:
            if ARTEFACT_GENES.match(row.get("GENE_SYMBOL", "")):
                continue
            if BENIGN.search(row.get("CLINVAR_PRIMARY_INTERPRETATION", "") or ""):
                continue
            key = (
                chr_name(row["CONTIG"]),
                int(row["START"]),
                row["REF"],
                row["ALT"],
            )
            gene_score = float(row.get("EXOMISER_GENE_COMBINED_SCORE") or 0.0)
            var_score = float(row.get("EXOMISER_VARIANT_SCORE") or 0.0)
            contributing = str(row.get("CONTRIBUTING_VARIANT", "0")) in ("1", "true", "True")
            score = gene_score + 0.01 * var_score + (0.05 if contributing else 0.0)

            moi = row.get("MOI", "?")
            seen_mois = best[key].mois if key in best else set()
            if key not in best or score > best[key].score:
                best[key] = Candidate(
                    variants=[key],
                    gene=row.get("GENE_SYMBOL", ""),
                    score=score,
                    rationale=f"Exomiser MOI {moi}; "
                              f"gene score {gene_score:.3f}, variant score {var_score:.3f}.",
                    mois=seen_mois,
                    gene_score=gene_score,
                    homozygous=row.get("GENOTYPE", "").replace("|", "/") in ("1/1", "2/2"),
                )
            best[key].mois.add(moi)
    return list(best.values())


def pair_compound_hets(candidates, phase_blocks=None):
    """Collapse same-gene candidates into compound-het pairs.

    In a singleton with no parents, phase is only resolvable by read-backed
    phasing when the two variants fall within one phase block. Where phase
    IS resolved and the variants are in cis, do not pair them -- a cis pair
    is not a compound heterozygote. Where phase is unresolved, pair anyway
    and say so in the rationale; an unproven pair is still worth the row.
    """
    in_cis = load_cis_pairs(phase_blocks) if phase_blocks else set()

    # Only variants Exomiser scored under a recessive mode can form a pair.
    by_gene = {}
    for c in candidates:
        if c.mois & RECESSIVE_MOIS:
            by_gene.setdefault(c.gene, []).append(c)

    paired, used = [], set()
    for gene, group in by_gene.items():
        group.sort(key=lambda c: -c.score)
        for i in range(len(group)):
            for j in range(i + 1, len(group)):
                a, b = group[i], group[j]
                key = (a.variants[0], b.variants[0])
                va, vb = a.variants[0], b.variants[0]
                if va[0] == vb[0] and abs(va[1] - vb[1]) < MIN_PAIR_DISTANCE:
                    continue   # same event / same haplotype, not compound het
                if key in in_cis or tuple(reversed(key)) in in_cis:
                    continue
                if id(a) in used or id(b) in used:
                    continue
                phase_note = (
                    "phase unresolved (variants outside a shared phase block)"
                    if phase_blocks else "phase not assessed"
                )
                paired.append(
                    Candidate(
                        variants=[a.variants[0], b.variants[0]],
                        gene=gene,
                        gene_score=max(a.gene_score, b.gene_score),
                        mois=a.mois | b.mois,
                        # Pairs score above either single: if the answer key
                        # is a pair, only this row can score a full match.
                        score=(a.score + b.score) / 2 + 0.10,
                        rationale=f"Biallelic candidate in {gene}; {phase_note}.",
                    )
                )
                used.update({id(a), id(b)})

    # A heterozygous variant that only fits a recessive model cannot explain the
    # disease alone: keep it only as a hedge for a gene that has a valid pair.
    paired_genes = {p.gene for p in paired}
    singles = [
        c for c in candidates
        if not (c.mois <= RECESSIVE_MOIS and not c.homozygous and c.gene not in paired_genes)
    ]
    return paired + singles


def load_cis_pairs(path):
    """Parse whatshap output for variant pairs proven to be in cis."""
    # Left as a stub: fill in once you have seen your phase block format.
    return set()


def assign_epcr(rows):
    """Strictly distinct, descending EPCR values.

    Primaries occupy the top of the range. Secondary findings are pushed
    below EPCR_FLOOR so they fall under any sensible F-max threshold and
    cannot count as false positives against your primary candidates.
    """
    primaries = [r for r in rows if r.finding_type == "primary"]
    secondaries = [r for r in rows if r.finding_type != "primary"]

    primaries.sort(key=lambda r: -r.score)
    secondaries.sort(key=lambda r: -r.score)

    step = (EPCR_TOP - 0.20) / max(len(primaries), 1)
    for i, r in enumerate(primaries):
        r.epcr = round(EPCR_TOP - i * step, 4)

    for i, r in enumerate(secondaries):
        r.epcr = round(EPCR_FLOOR - i * 0.001, 4)

    ordered = primaries + secondaries

    # Hard guarantee: no ties anywhere in the file.
    seen = set()
    for r in ordered:
        while r.epcr in seen:
            r.epcr = round(r.epcr - 0.0001, 4)
        seen.add(r.epcr)

    return ordered


def write_csv(rows, path, proband_id):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)   # comma-separated, quoted as needed
        w.writeheader()
        for r in rows:
            v1 = r.variants[0]
            v2 = r.variants[1] if r.is_pair else ("", "", "", "")
            w.writerow({
                "proband_id": proband_id,
                "chrom_1": v1[0], "pos_1": v1[1], "ref_1": v1[2], "alt_1": v1[3],
                "chrom_2": v2[0], "pos_2": v2[1], "ref_2": v2[2], "alt_2": v2[3],
                "epcr": r.epcr,
                "finding_type": r.finding_type,
                "notes": f"{r.gene}: {r.rationale}".strip(": "),
            })


def validate(path, proband_id):
    """Re-read the file with the same checks the Space's evaluation.py applies."""
    errors = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        missing = [c for c in COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            errors.append(f"missing columns: {missing}")
        rows = list(reader)
    if not 1 <= len(rows) <= MAX_ROWS:
        errors.append(f"{len(rows)} rows; must be 1-{MAX_ROWS}")
    epcrs = []
    for n, row in enumerate(rows, start=2):
        if row["proband_id"].strip() != proband_id:
            errors.append(f"line {n}: proband_id {row['proband_id']!r}")
        if not (row["chrom_1"] and row["pos_1"] and row["ref_1"] and row["alt_1"]):
            errors.append(f"line {n}: incomplete first variant")
        if not row["chrom_1"].startswith("chr"):
            errors.append(f"line {n}: chrom_1 {row['chrom_1']!r} lacks 'chr'")
        second = [row["chrom_2"], row["pos_2"], row["ref_2"], row["alt_2"]]
        if any(second) and not all(second):
            errors.append(f"line {n}: second variant partially filled")
        try:
            int(row["pos_1"])
            if row["pos_2"]:
                int(row["pos_2"])
        except ValueError:
            errors.append(f"line {n}: non-integer position")
        e = float(row["epcr"])
        if not 0 < e <= 1:
            errors.append(f"line {n}: epcr {e} outside (0, 1]")
        epcrs.append(e)
        if row["finding_type"] not in ("primary", "secondary"):
            errors.append(f"line {n}: finding_type {row['finding_type']!r}")
    if len(set(epcrs)) != len(epcrs):
        errors.append("duplicate EPCR values")
    if epcrs != sorted(epcrs, reverse=True):
        errors.append("rows not sorted by descending EPCR")
    return errors


def check_in_vcf(rows, vcf_path):
    """Every submitted variant must exist, exactly, in the callset it came from."""
    import pysam
    missing = []
    vcf = pysam.VariantFile(vcf_path)
    for r in rows:
        for chrom, pos, ref, alt in r.variants:
            hit = any(
                rec.pos == pos and rec.ref == ref and alt in (rec.alts or ())
                for rec in vcf.fetch(chrom, pos - 1, pos)
            )
            if not hit:
                missing.append(f"{chrom}:{pos} {ref}>{alt}")
    return missing


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ranked", required=True, help="Exomiser *.variants.tsv")
    p.add_argument("--phase", help="whatshap phase block TSV")
    p.add_argument("--vcf", help="normalized VCF the candidates came from (existence check)")
    p.add_argument("--curation", help="gene-level curation TSV (gene, action, note)")
    p.add_argument("--min-gene-score", type=float, default=DEFAULT_MIN_GENE_SCORE,
                   help="minimum Exomiser gene combined score for primary rows")
    p.add_argument("--proband-id", default=DEFAULT_PROBAND)
    p.add_argument("--out", required=True, help="output .csv")
    args = p.parse_args()

    candidates = load_ranked(args.ranked)
    if not candidates:
        sys.exit("No candidates in ranked input")

    candidates = pair_compound_hets(candidates, args.phase)
    candidates = apply_curation(candidates, load_curation(args.curation), args.min_gene_score)
    candidates.sort(key=lambda c: -c.score)
    rows = assign_epcr(candidates[:MAX_ROWS])
    write_csv(rows, args.out, args.proband_id)

    errors = validate(args.out, args.proband_id)
    if args.vcf:
        errors += [f"not in VCF: {m}" for m in check_in_vcf(rows, args.vcf)]
    if errors:
        print("VALIDATION FAILED:\n  " + "\n  ".join(errors), file=sys.stderr)
        sys.exit(1)

    pairs = sum(1 for r in rows if r.is_pair)
    print(f"OK: wrote {len(rows)} rows ({pairs} paired) to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()