#!/usr/bin/env python3
"""
ESM-2 masked-marginal scoring of BUB1B missense variants, calibrated against
experimentally characterised MVA alleles.

Why calibration: a protein language model score on its own is a number with no
scale. Suijkerbuijk et al. (Cancer Res 2010) measured which MVA missense alleles
reduce BubR1 protein stability and which do not. Scoring those alongside the
patient's allele tells us whether the model separates the two classes on this
protein -- and if it does, where the patient's allele falls.

Score: log P(mut) - log P(wt) at the masked position (Meier et al. 2021).
More negative = less compatible with the protein family. This measures
evolutionary/fitness constraint, NOT folding stability directly; the
calibration set is what licenses reading it as a stability proxy here.

Context: every possible substitution in the pseudokinase domain (766-1050) is
also scored, so the patient's allele gets a percentile within its domain.

ESM-2 accepts at most 1022 residues; BubR1 is 1050. Each position is scored in
a 1022-residue window placed as close to centred on it as the sequence ends
allow -- in practice the window drops at most 28 residues from one terminus.

Runs on CPU. Model weights (~2.5 GB) download once from dl.fbaipublicfiles.com.

Usually run by `make esm`. Directly:
  esm_variant_score.py --config findings/track2/case.yaml \\
      --fasta resources/O60566.fasta --out findings/track2/esm_scores.tsv
  (a gene file instead of the resolved case needs --patient, e.g. N1002K)
"""

import argparse
import os
import sys

AA = "ACDEFGHIKLMNPQRSTVWY"

# Calibration alleles and the domain map come from the case configuration
# (config/genes/<gene>.yaml); for BUB1B these are the six MVA alleles characterised
# by Suijkerbuijk et al. 2010 (Cancer Res 70:4891).
WINDOW = 1022


def read_fasta(path):
    seq = "".join(l.strip() for l in open(path) if not l.startswith(">"))
    if not seq:
        sys.exit(f"empty sequence in {path}")
    return seq


def parse(v):
    return v[0], int(v[1:-1]), v[-1]


def window(seq, pos):
    """0-based slice start for a WINDOW-residue window centred on 1-based pos."""
    if len(seq) <= WINDOW:
        return 0
    start = pos - 1 - WINDOW // 2
    return max(0, min(start, len(seq) - WINDOW))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", required=True, help="resolved case or gene file (YAML)")
    p.add_argument("--fasta", required=True)
    p.add_argument("--patient", help="e.g. N1002K (default: patient.missense in the case)")
    p.add_argument("--out", required=True)
    p.add_argument("--model", default="esm2_t33_650M_UR50D")
    p.add_argument("--threads", type=int, default=16)
    p.add_argument("--no-scan", action="store_true",
                   help="skip the domain-wide scan (calibration + patient only)")
    a = p.parse_args()

    import torch
    import esm
    torch.set_num_threads(a.threads)

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import caseconfig
    import nmd_check
    cfg = caseconfig.load(a.config)
    CALIBRATION = caseconfig.calibration(cfg)
    a.patient = a.patient or cfg["patient"].get("missense")
    if not a.patient:
        sys.exit("needs a patient missense (--patient or patient.missense in the case)")
    if not CALIBRATION:
        print("no calibration alleles: reporting the patient score and percentile "
              "only, without a class", file=sys.stderr)

    seq = read_fasta(a.fasta)
    print(f"sequence length {len(seq)}", file=sys.stderr)

    # Percentile context: the curated domain holding the patient residue,
    # else the whole protein.
    ppos = parse(a.patient)[1]
    hits = [d for d in (nmd_check.load_domains(cfg["domains"]) if cfg["domains"] else [])
            if d["start"] <= ppos <= d["end"] and d["end"] - d["start"] >= 20]
    DOMAIN = (hits[0]["start"], hits[0]["end"]) if hits else (1, len(seq))
    print(f"percentile context: residues {DOMAIN[0]}-{DOMAIN[1]}"
          f"{' (' + hits[0]['feature'] + ')' if hits else ' (whole protein)'}", file=sys.stderr)

    variants = [(a.patient, "PATIENT")] + CALIBRATION
    for v, _ in variants:
        wt, pos, mt = parse(v)
        if seq[pos - 1] != wt:
            sys.exit(f"{v}: reference has {seq[pos-1]} at {pos}, not {wt} -- "
                     f"wrong isoform or numbering. Stopping rather than scoring "
                     f"the wrong residue.")
    print("all variant reference residues match the sequence", file=sys.stderr)

    model, alphabet = getattr(esm.pretrained, a.model)()
    model.eval()
    conv = alphabet.get_batch_converter()

    cache = {}

    def logprobs_at(pos):
        """Log-probabilities over AA at 1-based pos, with pos masked."""
        if pos in cache:
            return cache[pos]
        s0 = window(seq, pos)
        sub = seq[s0:s0 + WINDOW]
        _, _, toks = conv([("x", sub)])
        i = pos - 1 - s0 + 1                    # +1 for BOS token
        toks[0, i] = alphabet.mask_idx
        with torch.no_grad():
            lp = torch.log_softmax(model(toks)["logits"][0, i], dim=-1)
        cache[pos] = {aa: float(lp[alphabet.get_idx(aa)]) for aa in AA}
        return cache[pos]

    def score(v):
        wt, pos, mt = parse(v)
        lp = logprobs_at(pos)
        return lp[mt] - lp[wt]

    rows = [(v, cls, score(v)) for v, cls in variants]

    domain_scores = []
    if not a.no_scan:
        lo, hi = DOMAIN
        hi = min(hi, len(seq))
        print(f"scanning {hi - lo + 1} positions x 19 substitutions ...",
              file=sys.stderr)
        for pos in range(lo, hi + 1):
            lp = logprobs_at(pos)
            wt = seq[pos - 1]
            domain_scores += [lp[m] - lp[wt] for m in AA if m != wt]
            if (pos - lo) % 50 == 0:
                print(f"  {pos}/{hi}", file=sys.stderr)

    with open(a.out, "w") as fh:
        fh.write("variant\tclass\tesm_llr\tdomain_percentile\n")
        for v, cls, s in rows:
            pct = (100 * sum(x <= s for x in domain_scores) / len(domain_scores)
                   if domain_scores else float("nan"))
            fh.write(f"{v}\t{cls}\t{s:.3f}\t{pct:.1f}\n")

    print("\nvariant   class          ESM LLR   pctile in domain", file=sys.stderr)
    for v, cls, s in sorted(rows, key=lambda r: r[2]):
        pct = (100 * sum(x <= s for x in domain_scores) / len(domain_scores)
               if domain_scores else float("nan"))
        print(f"{v:<9} {cls:<14} {s:+8.3f}   {pct:5.1f}", file=sys.stderr)

    d = [s for _, c, s in rows if c == "destabilising"]
    st = [s for _, c, s in rows if c == "stable"]
    if not (d and st):
        print(f"wrote {a.out}", file=sys.stderr)
        return
    sep = max(d) < min(st)
    print(f"\ncalibration: destabilising max {max(d):+.2f}, stable min "
          f"{min(st):+.2f} -> classes {'SEPARATE' if sep else 'OVERLAP'}",
          file=sys.stderr)
    if not sep:
        print("  The model does not cleanly separate the experimental classes on "
              "this protein;\n  the patient score should be reported as "
              "supporting context only.", file=sys.stderr)
    print(f"wrote {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
