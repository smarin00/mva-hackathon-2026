#!/usr/bin/env python3
"""
Predict nonsense-mediated decay (NMD) sensitivity of a premature termination
codon, from transcript geometry in an Ensembl GTF.

The question this answers is not "is the variant damaging" -- it is "does the
mutant transcript survive long enough to make protein". Those have opposite
therapeutic implications, and for BUB1B specifically they invert the sign of
the intervention (see --domains and the BUB1B note at the bottom).

Rules applied, in the order they are reported:

  1. 50-nt rule. A PTC more than 50 nt upstream of the last exon-exon junction
     triggers NMD. This is the dominant determinant. Boundaries in the
     literature range 50-55 nt; both are reported.
  2. Last-exon escape. A PTC in the final exon has no downstream junction, so
     no exon-junction complex remains to license decay. Escapes.
  3. Last-junction proximity. A PTC within 50 nt of the last junction escapes,
     for the same reason in weaker form.
  4. Start-proximal escape. A PTC within ~150 nt of the start codon can escape
     via translation reinitiation at a downstream AUG. Reported as a caveat,
     not a call, because reinitiation depends on downstream AUG context that
     this script does not model.
  5. Single-exon transcripts. No junctions at all, so never NMD-sensitive.

Everything reported as "ESCAPES" means the transcript persists; it does NOT
mean the protein is functional. Conversely "NMD-SENSITIVE" means the allele
contributes little or no protein, which for a recessive disorder usually makes
it the null allele and places the entire phenotype on the other allele.

Usage
  nmd_check.py --gtf /opt/refs/Homo_sapiens.GRCh38.116.gtf.gz \\
               --gene BUB1B --pos 40463318 \\
               --domains docs/bub1b_domains.tsv

  # explicit transcript instead of auto-selected canonical
  nmd_check.py --gtf ... --transcript ENST00000287598 --pos 40463318

Coordinates are 1-based genomic, matching VCF POS. Chromosome is taken from
the transcript, so a 'chr' prefix on input is neither needed nor rejected.
"""

import argparse
import gzip
import re
import sys

ATTR = {
    k: re.compile(rf'{k} "([^"]+)"')
    for k in ("gene_name", "gene_id", "transcript_id", "transcript_biotype", "tag")
}
NMD_WINDOW = 50          # canonical 50-nt rule
NMD_WINDOW_STRICT = 55   # conservative upper bound seen in the literature
REINIT_WINDOW = 150      # start-proximal reinitiation caveat


def opener(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def attr(field, key):
    m = ATTR[key].search(field)
    return m.group(1) if m else None


def tags(field):
    return set(re.findall(r'tag "([^"]+)"', field))


def scan(gtf, gene=None, transcript=None):
    """Collect exon and CDS intervals per transcript for the target gene."""
    tx = {}
    with opener(gtf) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] not in ("exon", "CDS", "transcript", "stop_codon"):
                continue
            a = f[8]
            if transcript:
                if attr(a, "transcript_id") != transcript:
                    continue
            else:
                if attr(a, "gene_name") != gene:
                    continue
            tid = attr(a, "transcript_id")
            if tid is None:
                continue
            rec = tx.setdefault(tid, {
                "chrom": f[0], "strand": f[6], "exons": [], "cds": [], "stop": False,
                "biotype": None, "tags": set(), "gene": attr(a, "gene_name"),
            })
            if f[2] == "transcript":
                rec["biotype"] = attr(a, "transcript_biotype")
                rec["tags"] = tags(a)
            elif f[2] == "stop_codon":
                rec["stop"] = True
            elif f[2] == "exon":
                rec["exons"].append((int(f[3]), int(f[4])))
            else:
                rec["cds"].append((int(f[3]), int(f[4])))
    return tx


def pick(tx):
    """Prefer MANE Select, then Ensembl canonical, then longest CDS."""
    coding = {t: r for t, r in tx.items()
              if r["biotype"] == "protein_coding" and r["cds"]}
    if not coding:
        raise SystemExit("error: no protein-coding transcript with CDS found")
    for want in ("MANE_Select", "Ensembl_canonical"):
        hits = [t for t, r in coding.items() if want in r["tags"]]
        if len(hits) == 1:
            return hits[0], coding[hits[0]], want
        if hits:  # tie: fall through to CDS length among the tagged set
            coding = {t: coding[t] for t in hits}
            break
    best = max(coding, key=lambda t: sum(e - s + 1 for s, e in coding[t]["cds"]))
    return best, coding[best], "longest-CDS"


def ordered(intervals, strand):
    """Exons in transcription order, 5' -> 3'."""
    return sorted(intervals, key=lambda x: x[0], reverse=(strand == "-"))


def to_transcript_coord(pos, exons, strand):
    """1-based position within the spliced transcript, or None if intronic."""
    offset = 0
    for s, e in ordered(exons, strand):
        if s <= pos <= e:
            within = (pos - s + 1) if strand == "+" else (e - pos + 1)
            return offset + within
        offset += e - s + 1
    return None


def aa_count(rec, cds_len):
    """Protein length. Ensembl GTF CDS features exclude the stop codon (it is
    a separate stop_codon feature); other GTFs may include it."""
    return cds_len // 3 if rec.get("stop") else cds_len // 3 - 1


def load_domains(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] == "start":
                continue
            rows.append({
                "start": int(f[0]), "end": int(f[1]), "feature": f[2],
                "role": f[3], "dispensable": f[4],
                "source": f[5] if len(f) > 5 else "",
            })
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gtf", required=True)
    p.add_argument("--gene", help="gene symbol (ignored if --transcript given)")
    p.add_argument("--transcript", help="explicit Ensembl transcript ID")
    p.add_argument("--pos", type=int, required=True,
                   help="1-based genomic position of the PTC (VCF POS)")
    p.add_argument("--domains", help="TSV of protein domains to intersect")
    a = p.parse_args()

    if not (a.gene or a.transcript):
        raise SystemExit("error: need --gene or --transcript")

    tx = scan(a.gtf, a.gene, a.transcript)
    if not tx:
        raise SystemExit(f"error: {a.transcript or a.gene} not found in {a.gtf}")

    if a.transcript:
        tid = a.transcript
        rec = tx[tid]
        how = "user-specified"
        if not rec["cds"]:
            raise SystemExit(f"error: {tid} has no CDS")
    else:
        tid, rec, how = pick(tx)

    strand = rec["strand"]
    exons = ordered(rec["exons"], strand)
    cds = ordered(rec["cds"], strand)
    n_exons = len(exons)
    tx_len = sum(e - s + 1 for s, e in exons)
    cds_len = sum(e - s + 1 for s, e in cds)

    print(f"gene        {rec['gene']}")
    print(f"transcript  {tid}  ({how})")
    print(f"locus       {rec['chrom']}:{min(s for s,_ in exons)}-"
          f"{max(e for _,e in exons)}  strand {strand}")
    print(f"structure   {n_exons} exons, {tx_len} nt transcript, "
          f"{cds_len} nt CDS ({aa_count(rec, cds_len)} aa)")
    print()

    # --- locate the variant ------------------------------------------------
    tpos = to_transcript_coord(a.pos, rec["exons"], strand)
    if tpos is None:
        print(f"position {a.pos} is INTRONIC in {tid} (or outside it).")
        print("A splice-region variant needs its effect on the transcript")
        print("resolved first -- rerun with the position of the resulting PTC.")
        return 2

    exon_idx = next(i for i, (s, e) in enumerate(exons, 1) if s <= a.pos <= e)

    cpos = to_transcript_coord(a.pos, rec["cds"], strand)
    if cpos is None:
        print(f"position {a.pos} is in exon {exon_idx} but OUTSIDE the CDS "
              f"(UTR). Not a PTC.")
        return 2
    codon = (cpos + 2) // 3

    if n_exons == 1:
        print("VERDICT  ESCAPES NMD -- single-exon transcript, no junctions.")
        return 0

    # last exon-exon junction, in transcript coordinates
    last_junction = sum(e - s + 1 for s, e in exons[:-1])
    dist = last_junction - tpos

    print(f"PTC in exon {exon_idx} of {n_exons}")
    print(f"  transcript coord     {tpos} / {tx_len}")
    print(f"  CDS coord            {cpos} / {cds_len}  -> codon {codon}")
    print(f"  last exon-exon jct   {last_junction}")
    print(f"  distance to last jct {dist} nt "
          f"({'upstream' if dist > 0 else 'downstream -- in final exon'})")
    print()

    # --- rules -------------------------------------------------------------
    if exon_idx == n_exons:
        verdict = "ESCAPES NMD"
        why = "PTC lies in the final exon; no downstream exon-junction complex."
    elif dist > NMD_WINDOW_STRICT:
        verdict = "NMD-SENSITIVE"
        why = (f"PTC is {dist} nt upstream of the last junction, beyond both the "
               f"{NMD_WINDOW}-nt and the conservative {NMD_WINDOW_STRICT}-nt "
               f"boundary.")
    elif dist > NMD_WINDOW:
        verdict = "NMD-SENSITIVE (borderline)"
        why = (f"PTC is {dist} nt upstream of the last junction: past the "
               f"{NMD_WINDOW}-nt rule but inside the {NMD_WINDOW_STRICT}-nt "
               f"conservative bound. Treat as uncertain and confirm on RNA.")
    else:
        verdict = "ESCAPES NMD"
        why = (f"PTC is only {dist} nt upstream of the last junction, inside "
               f"the {NMD_WINDOW}-nt window.")

    print(f"VERDICT  {verdict}")
    print(f"         {why}")

    if cpos <= REINIT_WINDOW:
        print()
        print(f"CAVEAT   PTC is {cpos} nt into the CDS, within the "
              f"{REINIT_WINDOW}-nt start-proximal window where translation")
        print("         reinitiation at a downstream AUG can rescue the "
              "transcript from decay. Not modelled here.")

    # --- consequence for the protein --------------------------------------
    aa_total = aa_count(rec, cds_len)
    lost = aa_total - (codon - 1)
    print()
    print(f"If translated: {codon - 1} of {aa_total} aa retained "
          f"({100*(codon-1)/aa_total:.1f}%), {lost} aa lost.")

    if a.domains:
        doms = load_domains(a.domains)
        ko = [d for d in doms if d["end"] >= codon]
        kept = [d for d in doms if d["end"] < codon]
        print()
        print("Features RETAINED (end before the PTC):")
        for d in kept or [None]:
            print(f"  {d['feature']:<14} {d['start']}-{d['end']}" if d
                  else "  (none)")
        print("Features LOST or TRUNCATED (span or follow the PTC):")
        for d in ko or [None]:
            if d is None:
                print("  (none)")
                continue
            flag = {"yes": "dispensable for SAC",
                    "no": "REQUIRED for SAC",
                    "partial": "partially required",
                    "untested": "untested"}.get(d["dispensable"], d["dispensable"])
            print(f"  {d['feature']:<14} {d['start']}-{d['end']}   [{flag}]")

        essential_lost = [d["feature"] for d in ko if d["dispensable"] == "no"]
        print()
        if essential_lost:
            print(f"=> Truncated product would be NON-FUNCTIONAL: loses "
                  f"{', '.join(essential_lost)}.")
        else:
            print("=> Truncated product retains every SAC-essential feature. "
                  "If it escaped decay")
            print("   and were stable, it could be partially functional "
                  "(cf. BUB1B 731X, which")
            print("   fully rescues checkpoint and alignment despite lacking "
                  "the pseudokinase domain).")

    return 0


if __name__ == "__main__":
    sys.exit(main())