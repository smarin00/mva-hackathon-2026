#!/usr/bin/env python3
"""
Load configuration. Three kinds of YAML file:

  gene file      config/genes/<gene>.yaml: gene knowledge (uniprot, domains,
                 calibration, pathway panel). Has `uniprot`.
  resolved case  findings/track2/case.yaml, written by resolve_case.py: the
                 patient's variants with their consequences, plus a pointer
                 to the gene file. Has `gene_config`.
  case file      params.yaml: everything about one patient, for both tracks.

load() returns the gene knowledge merged with the patient block (empty for a
gene file), which is what the Track 2 scripts read.

As a module:   cfg = caseconfig.load("config/genes/bub1b.yaml")
From make:     caseconfig.py FILE KEY      e.g. caseconfig.py params.yaml mode
               (lists are space-joined; a missing key prints nothing)
               caseconfig.py --make findings/track2/case.yaml > case.mk
"""
import sys

import yaml


def _read(path):
    with open(path) as fh:
        return yaml.safe_load(fh) or {}


def load(path):
    y = _read(path)
    if "gene_config" in y:                      # resolved case
        cfg = load(y["gene_config"])
        cfg["patient"] = y.get("patient") or {}
        cfg["variants"] = y.get("variants") or []
        return cfg
    if "uniprot" not in y:
        sys.exit(f"{path}: not a gene file (no 'uniprot') or resolved case "
                 f"(no 'gene_config')")
    if "gene" not in y:
        sys.exit(f"{path}: required key 'gene' missing")
    for k, v in (("calibration", []), ("pathway", []), ("named_drugs", []),
                 ("stability_features", []), ("patient", {}), ("domains", ""),
                 ("variants", [])):
        if y.get(k) is None:
            y[k] = v
    return y


def calibration(cfg):
    """[(variant, class), ...] in the order given."""
    return [(c["variant"], c["class"]) for c in cfg["calibration"]]


def panel(cfg):
    """[(symbol, desired, tier, rationale), ...]"""
    return [(p["symbol"], p["desired"], str(p.get("tier", "1")), p.get("rationale", ""))
            for p in cfg["pathway"]]


def get(cfg, dotted):
    cur = cfg
    for part in dotted.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return None
        cur = cur[part]
    return cur


MAKE_VARS = [("GENE", "gene"), ("UNIPROT", "uniprot"), ("DOMAINS", "domains"),
             ("PATIENT_CHROM", "patient.chrom"), ("PATIENT_PTC_POS", "patient.ptc_pos"),
             ("PATIENT_PTC_ALLELE", "patient.ptc_allele"),
             ("PATIENT_MISSENSE", "patient.missense"),
             ("PATIENT_MIS_ALLELE", "patient.missense_allele")]


def make_vars(path):
    """Resolved case as make assignments (included by the Makefile)."""
    cfg = load(path)
    out = [f"# generated from {path}"]
    for var, key in MAKE_VARS:
        v = get(cfg, key)
        out.append(f"{var} := {'' if v is None else v}")
    out.append(f"PANEL := {'yes' if cfg['pathway'] else ''}")
    return "\n".join(out) + "\n"


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "--make":
        sys.stdout.write(make_vars(sys.argv[2]))
        sys.exit(0)
    if len(sys.argv) != 3:
        sys.exit("usage: caseconfig.py FILE KEY | caseconfig.py --make RESOLVED_CASE")
    raw = _read(sys.argv[1])
    cfg = load(sys.argv[1]) if ("uniprot" in raw or "gene_config" in raw) else raw
    v = get(cfg, sys.argv[2])
    if v is None or v == "":
        sys.exit(0)
    print(" ".join(map(str, v)) if isinstance(v, list) else v)
