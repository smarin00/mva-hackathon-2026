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
    "liver": r"hepatotoxicity|hepatic failure|liver injury|liver failure|liver dysfunction|"
             r"hepatic dysfunction|(?<!pulmonary )veno-occlusive|hepatic injury|"
             r"sinusoidal obstruction|transaminase",
    "heart/QT": r"qt prolongation|qtc prolongation|prolong the qt|cardiomyopathy|"
                r"heart failure|cardiotoxicity|cardiac toxicity",
    "nerves": r"peripheral neuropathy|neurotoxicity|neuropathy",
    "pancreas": r"pancreatitis",
    "kidney": r"nephrotoxicity|renal failure|renal toxicity|acute kidney injury",
}
# CYP3A: the sentence must state it of the drug itself ("X is a substrate of /
# is metabolized by / is an inhibitor of CYP3A"), never a negation ("X is not an
# inhibitor of ...") and never a heading or advice about other drugs.
NEGATED = re.compile(r"\b(not|no|neither|nor)\b[^.]{0,40}(inhibit|induc|substrate|metaboli)", re.I)


def cyp_rx(drug, role):
    d = re.escape(drug)
    if role == "substrate":
        verb = r"(is|are)[^.]{0,40}(substrate|metaboli[sz]ed)[^.]{0,60}(cyp\s?3a|cytochrome p450 3a)"
    else:
        verb = r"(is|acts as)[^.]{0,30}(inhibitor|inducer)[^.]{0,30}(cyp\s?3a|cytochrome p450 3a)"
    return re.compile(rf"\b{d}\b[^.]{{0,60}}?{verb}", re.I)


def sentences(text):
    return [s.strip() for s in re.split(r"(?<=[.;])\s+", text) if s.strip()]


def section_text(label, fields):
    return " ".join(" ".join(label.get(f, [])) for f in fields)


# A keyword preceded by these is a risk factor or a symptom list, not a
# toxicity of the drug ("in patients with renal failure").
RISK_CONTEXT = re.compile(r"(patients? with|history of|in patients|pre-?existing|"
                          r"suggestive|symptoms|suspected|risk factors?).{0,300}$", re.I)
# (the text searched is already a single sentence, so '.' here is a decimal point)


def find(label, fields, pattern, toxicity=False):
    """First label sentence (in `fields`) matching pattern, shortened."""
    rx = pattern if hasattr(pattern, "search") else re.compile(pattern, re.I)
    for s in sentences(section_text(label, fields)):
        m = rx.search(s)
        if not m:
            continue
        if toxicity and RISK_CONTEXT.search(s[:m.start()]):
            continue
        if not toxicity and NEGATED.search(s):
            continue
        return s if len(s) <= 300 else s[:297] + "..."
    return ""


SALTS = {"hydrochloride", "hcl", "sulfate", "sodium", "citrate", "acetate", "phosphate",
         "mesylate", "maleate", "tartrate", "monohydrate", "anhydrous", "extended", "release",
         "for", "injection", "usp", "and", "calcium", "potassium"}


LABEL_TEXT = ("description", "boxed_warning", "package_label_principal_display_panel",
              "spl_product_data_elements", "indications_and_usage")


def plain(label, drug):
    """Generic name is the drug itself (plus salt words): not a liposomal
    formulation (checked in the label text too: some liposomal generics are
    named like the plain drug) and not a combination product."""
    if "liposom" in " ".join(" ".join(label.get(k, [])) for k in LABEL_TEXT).lower():
        return False
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
    rx = any("PRESCRIPTION" in t for t in best.get("openfda", {}).get("product_type", []))
    # only prescription labels count: approved medicines, not OTC/homeopathic listings
    return best if best and plain(best, drug) and rx else {}


def query(drug):
    """openFDA search with retries: 404 = no label; 5xx and timeouts are
    retried with backoff, then with a smaller result page."""
    import time
    # Same form as openFDA's documented examples: literal ':' , quoted phrase,
    # '+' for spaces (percent-encoding the ':' breaks the query).
    q = f'openfda.generic_name:%22{urllib.parse.quote_plus(drug)}%22'
    last = None
    for limit in (100, 25):
        for attempt in range(3):
            req = urllib.request.Request(f"{API}?search={q}&limit={limit}",
                                         headers={"User-Agent": "mva-hackathon-pipeline"})
            try:
                with urllib.request.urlopen(req, timeout=60) as r:
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
    if label:                      # never cache a miss: it may be transient
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
    sub = find(d_label, ALL_FIELDS, cyp_rx(chemo, "substrate"))
    mod = find(c_label, ALL_FIELDS, cyp_rx(cand, "modulator"))
    if sub and mod:
        flags.append(("cyp3a4", "", mod, sub))
    for organ, pat in TOXICITY.items():
        a = find(c_label, WARN_FIELDS, pat, toxicity=True)
        b = find(d_label, WARN_FIELDS, pat, toxicity=True)
        if a and b:
            flags.append(("toxicity", organ, a, b))
    return flags


def label_info(drug, label):
    if not label:
        return {"drug": drug, "status": "no US prescription label (supplement or OTC)",
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
            print(f"  {c:<24} no US prescription label: not checked", file=sys.stderr)
            continue
        for d in chemo_regs:
            f = [r for r in rows if r[0] == c and r[1] == d]
            s = ", ".join(r[3] + (f"({r[4]})" if r[4] else "") for r in f) or "-"
            print(f"  {c:<24} {d:<18} {s}", file=sys.stderr)
    print(f"wrote {args.outdir}/compatibility.tsv, labels.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()
