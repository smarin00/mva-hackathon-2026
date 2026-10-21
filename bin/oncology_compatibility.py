#!/usr/bin/env python3
"""
Are the drug candidates compatible with the chemotherapy this child would
receive if one of the cancers the disease predisposes to develops?

Every check is taken from the current US prescribing information (FDA labels
via openFDA), not from memory:

  named      one label names the other drug
  cyp3a4     the chemotherapy drug is a CYP3A4 substrate and the candidate's
             label describes CYP3A4 inhibition or induction (vincristine
             neurotoxicity rises with CYP3A4 inhibitors)
  toxicity   both labels warn of the same organ toxicity (marrow, liver,
             heart/QT, nerves, pancreas, kidney) in their boxed warning or
             warnings sections

Each flag carries the label sentences it rests on. Flags are prompts for
review, not verdicts: an overlapping warning may be clinically trivial, and a
missing flag is not proof of safety. Candidates without a US drug label
(dietary supplements) are reported as such.

Tumours, regimens and candidates come from the `oncology` block of the gene
file. Labels are cached in resources/labels/; their set id, version and date
are written to labels.tsv so the exact versions can be retrieved (DailyMed).

Usage
  oncology_compatibility.py --config findings/track2/case.yaml \\
      --outdir findings/track2/oncology [--labels resources/labels] [--offline]
"""
import argparse
import json
import os
import re
import sys
import urllib.error
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import caseconfig  # noqa: E402

API = "https://api.fda.gov/drug/label.json"
WARN_FIELDS = ["boxed_warning", "warnings_and_cautions", "warnings", "precautions"]
ALL_FIELDS = WARN_FIELDS + ["drug_interactions", "clinical_pharmacology",
                            "pharmacokinetics", "adverse_reactions", "contraindications"]

TOXICITY = {
    "marrow": r"myelosuppression|bone marrow suppression|neutropenia|thrombocytopenia|"
              r"pancytopenia|agranulocytosis",
    "liver": r"hepatotoxicity|hepatic failure|liver injury|liver failure|"
             r"veno-occlusive|hepatic injury|sinusoidal obstruction",
    "heart/QT": r"qt prolongation|qtc prolongation|prolong the qt|cardiomyopathy|"
                r"heart failure|cardiotoxicity|cardiac toxicity",
    "nerves": r"peripheral neuropathy|neurotoxicity|neuropathy",
    "pancreas": r"pancreatitis",
    "kidney": r"nephrotoxicity|renal failure|renal toxicity|acute kidney injury",
}
CYP_SUB = r"cyp3a4?[^.]{0,120}(substrate|metaboli[sz])|(substrate|metaboli[sz])[^.]{0,120}cyp3a4?"
CYP_MOD = r"(inhibitor|inducer|inhibits|induces)[^.]{0,60}cyp3a4?|cyp3a4?[^.]{0,30}(inhibitor|inducer)"


def about(drug, pattern):
    """Sentences stating `pattern` of the drug itself: they name the drug and
    are not advice about other drugs ("Avoid strong CYP3A4 inhibitors")."""
    return re.compile(rf"^(?!\s*(avoid|do not|consider|if|when|concomitant))"
                      rf"(?=.*\b{re.escape(drug)}\b).*(?:{pattern})", re.I)


def sentences(text):
    return [s.strip() for s in re.split(r"(?<=[.;])\s+", text) if s.strip()]


def section_text(label, fields):
    return " ".join(" ".join(label.get(f, [])) for f in fields)


def find(label, fields, pattern):
    """First label sentence (in `fields`) matching pattern, shortened."""
    rx = pattern if hasattr(pattern, "search") else re.compile(pattern, re.I)
    for s in sentences(section_text(label, fields)):
        if rx.search(s):
            return s if len(s) <= 300 else s[:297] + "..."
    return ""


SALTS = {"hydrochloride", "hcl", "sulfate", "sodium", "citrate", "acetate", "phosphate",
         "mesylate", "maleate", "tartrate", "monohydrate", "anhydrous", "extended", "release",
         "for", "injection", "usp", "and", "calcium", "potassium"}


def plain(label, drug):
    """Generic name is the drug itself (plus salt words): not a liposomal
    formulation and not a combination product."""
    want = set(drug.lower().split())
    for g in label.get("openfda", {}).get("generic_name", []):
        toks = set(re.findall(r"[a-z0-9]+", g.lower()))
        if want <= toks and not (toks - want - SALTS) and "liposom" not in g.lower():
            return True
    return False


def pick(results, drug):
    """Prefer: plain generic name > prescription > has warnings > most recent."""
    def key(r):
        rx = any("PRESCRIPTION" in t for t in r.get("openfda", {}).get("product_type", []))
        return (plain(r, drug), rx, any(r.get(f) for f in WARN_FIELDS), r.get("effective_time", ""))
    best = max(results, key=key) if results else {}
    return best if best and plain(best, drug) else {}


def query(drug):
    """openFDA search with retries: 404 = no label; 5xx and timeouts are
    retried with backoff, then with a smaller result page."""
    import time
    q = urllib.parse.quote(f'openfda.generic_name:"{drug}"')
    last = None
    for limit in (100, 25):
        for attempt in range(3):
            try:
                with urllib.request.urlopen(f"{API}?search={q}&limit={limit}",
                                            timeout=60) as r:
                    return json.load(r).get("results", [])
            except urllib.error.HTTPError as e:
                if e.code == 404:
                    return []
                if e.code < 500 and e.code != 429:
                    sys.exit(f"ERROR: openFDA query failed for {drug}: {e}")
                last = e
            except (urllib.error.URLError, TimeoutError) as e:
                last = e
            time.sleep(2 ** attempt)
        print(f"  openFDA: {drug}: {last}; retrying with limit 25", file=sys.stderr)
    sys.exit(f"ERROR: openFDA query failed for {drug} after retries: {last}")


def fetch(drug, cache, offline):
    path = os.path.join(cache, re.sub(r"\W+", "_", drug.lower()) + ".json")
    if os.path.exists(path):
        return json.load(open(path))
    if offline:
        return None
    results = query(drug)
    label = pick(results, drug)
    os.makedirs(cache, exist_ok=True)
    json.dump(label, open(path, "w"))
    return label


def check(cand, c_label, chemo, d_label):
    flags = []
    c_name = re.compile(rf"\b{re.escape(chemo)}\b", re.I)
    d_name = re.compile(rf"\b{re.escape(cand)}\b", re.I)
    a, b = find(c_label, ALL_FIELDS, c_name), find(d_label, ALL_FIELDS, d_name)
    if a or b:
        flags.append(("named", "", a, b))
    sub = find(d_label, ALL_FIELDS, about(chemo, CYP_SUB))
    mod = find(c_label, ALL_FIELDS, about(cand, CYP_MOD))
    if sub and mod:
        flags.append(("cyp3a4", "", mod, sub))
    for organ, pat in TOXICITY.items():
        a, b = find(c_label, WARN_FIELDS, pat), find(d_label, WARN_FIELDS, pat)
        if a and b:
            flags.append(("toxicity", organ, a, b))
    return flags


def label_info(drug, label):
    if not label:
        return {"drug": drug, "status": "no US drug label (e.g. dietary supplement)",
                "brand": "", "set_id": "", "version": "", "effective_time": ""}
    o = label.get("openfda", {})
    return {"drug": drug, "status": "label",
            "brand": ", ".join(o.get("brand_name", [])[:2]),
            "set_id": label.get("set_id", ""), "version": label.get("version", ""),
            "effective_time": label.get("effective_time", "")}


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--config", required=True, help="resolved case or gene file")
    a.add_argument("--outdir", required=True)
    a.add_argument("--labels", default="resources/labels")
    a.add_argument("--offline", action="store_true", help="use cached labels only")
    args = a.parse_args()

    onc = caseconfig.load(args.config).get("oncology") or {}
    regimens, cands = onc.get("regimens") or {}, onc.get("candidates") or []
    if not regimens or not cands:
        sys.exit(f"{args.config}: no oncology.regimens / oncology.candidates")
    chemo_regs = {}
    for reg, drugs in regimens.items():
        for d in drugs:
            chemo_regs.setdefault(d, []).append(reg)

    labels = {d: fetch(d, args.labels, args.offline) for d in list(chemo_regs) + cands}
    missing = [d for d, l in labels.items() if l is None]
    if missing:
        sys.exit(f"--offline but no cached label for: {', '.join(missing)}")
    nolabel_chemo = [d for d in chemo_regs if not labels[d]]
    if nolabel_chemo:
        sys.exit(f"no openFDA label found for chemotherapy drug(s) {nolabel_chemo}: "
                 f"check the generic names in the config")

    os.makedirs(args.outdir, exist_ok=True)
    with open(os.path.join(args.outdir, "labels.tsv"), "w") as fh:
        cols = ["drug", "status", "brand", "set_id", "version", "effective_time"]
        fh.write("\t".join(cols) + "\n")
        for d, l in labels.items():
            info = label_info(d, l)
            fh.write("\t".join(str(info[c]) for c in cols) + "\n")

    rows = []
    for c in cands:
        for d, regs in chemo_regs.items():
            if not labels[c]:
                continue
            for kind, organ, ev_c, ev_d in check(c, labels[c], d, labels[d]):
                rows.append([c, d, ",".join(regs), kind, organ, ev_c, ev_d])
    with open(os.path.join(args.outdir, "compatibility.tsv"), "w") as fh:
        fh.write("candidate\tchemotherapy\tregimens\tcheck\torgan\t"
                 "candidate_label\tchemotherapy_label\n")
        for r in rows:
            fh.write("\t".join(x.replace("\t", " ").replace("\n", " ") for x in r) + "\n")

    # summary: one line per candidate x chemotherapy
    print("candidate x chemotherapy: flags (see compatibility.tsv for the label text)",
          file=sys.stderr)
    for c in cands:
        if not labels[c]:
            print(f"  {c:<24} no US drug label: not checked", file=sys.stderr)
            continue
        for d in chemo_regs:
            f = [r for r in rows if r[0] == c and r[1] == d]
            s = ", ".join(r[3] + (f"({r[4]})" if r[4] else "") for r in f) or "-"
            print(f"  {c:<24} {d:<18} {s}", file=sys.stderr)
    print(f"wrote {args.outdir}/compatibility.tsv, labels.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()
