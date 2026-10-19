#!/usr/bin/env python3
"""
Automated candidate layer for Track 2: map a BubR1-dosage pathway panel to
known drugs in the Open Targets Platform, and classify every drug by whether
its action moves the target in the therapeutically desired direction.

The biology is encoded explicitly in PANEL below: each gene carries the
direction we want (raise BubR1 dose, or do not weaken the checkpoint) and the
literature it rests on. Open Targets supplies the drugs and their action type
(ChEMBL mechanism annotations). The cross of the two gives:

  candidate     drug acts in the desired direction
  do_no_harm    drug acts against it (e.g. any HSP90 inhibitor)
  review        action type ambiguous or not annotated

This is the reproducible, automated part of candidate identification; the
curated shortlist in the report is built on top of its output, not instead
of it.

Named drugs (NAD+ precursors, arimoclomol, senolytics) are also looked up
directly, because several of them act on a pathway rather than through a
single annotated target and would be missed by the target query.

Standard library only. Network access to api.platform.opentargets.org.

Usage
  opentargets_candidates.py --outdir findings/track2/opentargets
  opentargets_candidates.py --introspect Target   # print schema fields
"""

import argparse
import json
import sys
import time
import urllib.error
import urllib.request

API = "https://api.platform.opentargets.org/api/v4/graphql"

# ---------------------------------------------------------------------------
# The biology. desired: 'up' | 'down' | 'no_inhibit'
# ---------------------------------------------------------------------------
PANEL = [
    # (symbol, desired, tier, rationale)
    ("BUB1B",    "up",         "1", "The deficient protein; MVA1 is a BubR1 dosage disorder (Suijkerbuijk 2010)."),
    ("SIRT2",    "up",         "1", "Deacetylates BubR1 K668 and raises BubR1 abundance; extends lifespan of BubR1 hypomorphic mice (North 2014)."),
    ("NAMPT",    "up",         "1", "Rate-limiting NAD+ salvage enzyme; NAD+ supports SIRT2 activity on BubR1 (North 2014)."),
    ("CREBBP",   "down",       "1", "CBP acetylates BubR1 K668, opposing SIRT2 (North 2014). Caution: CREBBP loss causes Rubinstein-Taybi syndrome."),
    ("HSF1",     "up",         "1", "Master regulator of the heat-shock response; chaperone support for a misfolding-prone client."),
    ("HSPA1A",   "up",         "1", "HSP70; folding support for destabilised clients."),
    ("HSP90AA1", "no_inhibit", "1", "Mutant BubR1 is an HSP90 client; HSP90 inhibition depletes it (Suijkerbuijk 2010)."),
    ("HSP90AB1", "no_inhibit", "1", "HSP90 beta; as above."),
    ("TTK",      "no_inhibit", "1", "MPS1 kinase; inhibition weakens the spindle assembly checkpoint further."),
    ("AURKB",    "no_inhibit", "1", "Aurora B; inhibition impairs error correction and the checkpoint."),
    ("BCL2L1",   "down",       "2", "BCL-xL; senolytic target. Senescent-cell clearance delayed ageing phenotypes in BubR1 progeroid mice (Baker 2011)."),
]

NAMED_DRUGS = [
    "arimoclomol",
    "nicotinamide riboside",
    "nicotinamide mononucleotide",
    "nicotinamide",
    "dasatinib",
    "quercetin",
    "fisetin",
    "navitoclax",
]

UP = {"AGONIST", "ACTIVATOR", "POSITIVE ALLOSTERIC MODULATOR", "PARTIAL AGONIST",
      "OPENER", "STABILISER", "POSITIVE MODULATOR", "INDUCER"}
DOWN = {"INHIBITOR", "ANTAGONIST", "NEGATIVE ALLOSTERIC MODULATOR", "BLOCKER",
        "DEGRADER", "NEGATIVE MODULATOR", "INVERSE AGONIST", "ANTISENSE INHIBITOR",
        "RNAI INHIBITOR", "DISRUPTING AGENT", "HYDROLYTIC ENZYME"}


def direction(action_type):
    a = (action_type or "").upper().strip()
    if a in UP:
        return "up"
    if a in DOWN:
        return "down"
    return "unclear"


def verdict(desired, drug_dir):
    if drug_dir == "unclear":
        return "review"
    if desired == "no_inhibit":
        return "do_no_harm" if drug_dir == "down" else "review"
    return "candidate" if drug_dir == desired else "do_no_harm"


# ---------------------------------------------------------------------------
# API
# ---------------------------------------------------------------------------
class GQLError(RuntimeError):
    pass


def gql(query, variables=None, retries=3):
    body = json.dumps({"query": query, "variables": variables or {}}).encode()
    req = urllib.request.Request(API, data=body,
                                 headers={"Content-Type": "application/json"})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                out = json.load(r)
            if out.get("errors"):
                raise GQLError("; ".join(e.get("message", str(e)) for e in out["errors"]))
            return out["data"]
        except urllib.error.HTTPError as e:
            msg = e.read().decode(errors="replace")[:500]
            if e.code < 500 or attempt == retries - 1:
                raise GQLError(f"HTTP {e.code}: {msg}")
        except urllib.error.URLError as e:
            if attempt == retries - 1:
                raise GQLError(f"network: {e.reason}")
        time.sleep(2 * (attempt + 1))


Q_META = "query { meta { name apiVersion { x y z } dataVersion { year month } } }"

Q_SEARCH = """
query($q: String!, $e: [String!]!) {
  search(queryString: $q, entityNames: $e, page: {index: 0, size: 10}) {
    hits { id name entity }
  }
}"""

# Open Targets 26.x: Target.knownDrugs was replaced by
# Target.drugAndClinicalCandidates (rows carry the drug's MoA and the maximum
# clinical stage reached for that target).
Q_CANDIDATES = """
query($id: String!) {
  target(ensemblId: $id) {
    approvedSymbol
    drugAndClinicalCandidates {
      count
      rows {
        maxClinicalStage
        drug {
          id
          name
          drugType
          mechanismsOfAction {
            rows { mechanismOfAction actionType targets { id approvedSymbol } }
          }
        }
        diseases { disease { name } }
      }
    }
  }
}"""

# Drug-level status fields have moved between releases, so the named-drug
# query is assembled at run time from whichever of these the schema offers.
DRUG_OPTIONAL = ["isApproved", "maximumClinicalStage", "maximumClinicalTrialPhase",
                 "hasBeenWithdrawn", "drugType"]


def drug_query(available):
    opt = "\n    ".join(f for f in DRUG_OPTIONAL if f in available)
    return """
query($id: String!) {
  drug(chemblId: $id) {
    id
    name
    %s
    mechanismsOfAction {
      rows { mechanismOfAction actionType targets { id approvedSymbol } }
    }
  }
}""" % opt


def type_fields(type_name):
    q = "query($n: String!) { __type(name: $n) { fields { name } } }"
    t = gql(q, {"n": type_name})["__type"]
    return {f["name"] for f in (t or {}).get("fields", [])}


def stage_rank(stage):
    """Sortable rank for a clinical stage given as a number or a label."""
    if stage is None:
        return 0.0
    u = str(stage).upper()
    if "APPROV" in u:
        return 4.5
    digits = [int(c) for c in u if c.isdigit()]
    return float(max(digits)) if digits else 0.0


def resolve(name, entity):
    hits = gql(Q_SEARCH, {"q": name, "e": [entity]})["search"]["hits"]
    for h in hits:
        if h["name"].lower() == name.lower():
            return h["id"], h["name"]
    return (hits[0]["id"], hits[0]["name"]) if hits else (None, None)


def candidates_for(ensembl_id):
    t = gql(Q_CANDIDATES, {"id": ensembl_id})["target"] or {}
    block = t.get("drugAndClinicalCandidates") or {}
    rows = block.get("rows") or []
    if block.get("count") and block["count"] > len(rows):
        print(f"    note: {block['count']} candidates reported, {len(rows)} returned",
              file=sys.stderr)
    return rows


def action_on(drug, target_id, symbol):
    """actionType(s) of this drug on this target, from its MoA records."""
    acts = set()
    for m in ((drug or {}).get("mechanismsOfAction") or {}).get("rows", []) or []:
        for t in m.get("targets") or []:
            if t.get("id") == target_id or t.get("approvedSymbol") == symbol:
                acts.add(m.get("actionType") or "")
    return acts


def introspect(type_name):
    q = """query($n: String!) { __type(name: $n) { name fields { name
           args { name } type { name kind ofType { name kind } } } } }"""
    t = gql(q, {"n": type_name})["__type"]
    if not t:
        sys.exit(f"no type {type_name}")
    for f in t["fields"]:
        ty = f["type"]["name"] or (f["type"]["ofType"] or {}).get("name")
        args = ",".join(a["name"] for a in f["args"])
        print(f"{f['name']:<32} {ty or '':<28} ({args})")


# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", default="findings/track2/opentargets")
    p.add_argument("--introspect", metavar="TYPE",
                   help="print the fields of a GraphQL type and exit (debug)")
    a = p.parse_args()

    if a.introspect:
        introspect(a.introspect)
        return 0

    import os
    os.makedirs(a.outdir, exist_ok=True)

    try:
        m = gql(Q_META)["meta"]
        av, dv = m["apiVersion"], m["dataVersion"]
        release = f"API {av['x']}.{av['y']}.{av['z']}, data {dv['year']}.{dv['month']:0>2}"
    except GQLError as e:
        release = f"unknown ({e})"
    print(f"Open Targets Platform: {release}", file=sys.stderr)

    # ---- target -> drugs -------------------------------------------------
    out_rows = []
    for sym, desired, tier, why in PANEL:
        tid, tname = resolve(sym, "target")
        if not tid:
            print(f"  {sym}: not found", file=sys.stderr)
            continue
        rows = candidates_for(tid)
        seen = {}
        for r in rows:
            drug = r.get("drug") or {}
            did = drug.get("id")
            if not did:
                continue
            d = seen.setdefault(did, {"drug": drug, "diseases": set(), "stage": None})
            for dz in r.get("diseases") or []:
                d["diseases"].add(((dz or {}).get("disease") or {}).get("name") or "")
            if stage_rank(r.get("maxClinicalStage")) >= stage_rank(d["stage"]):
                d["stage"] = r.get("maxClinicalStage")
        for did, d in seen.items():
            drug = d["drug"]
            acts = action_on(drug, tid, sym)
            dirs = {direction(x) for x in acts} or {"unclear"}
            ddir = dirs.pop() if len(dirs) == 1 else "unclear"
            moas = ((drug.get("mechanismsOfAction") or {}).get("rows") or [])
            out_rows.append({
                "target": sym, "tier": tier, "desired": desired,
                "drug": drug.get("name") or did,
                "chembl_id": did,
                "action_type": ";".join(sorted(acts)) or "",
                "drug_direction": ddir,
                "verdict": verdict(desired, ddir),
                "approved": stage_rank(d["stage"]) >= 4,
                "max_stage": d["stage"],
                "mechanism": "; ".join(sorted({m.get("mechanismOfAction") or ""
                                               for m in moas}))[:200],
                "indications": "; ".join(sorted(x for x in d["diseases"] if x))[:300],
                "rationale": why,
            })
        print(f"  {sym:<9} {tid}  {len(seen):>3} drugs", file=sys.stderr)

    cols = ["target", "tier", "desired", "drug", "chembl_id", "action_type",
            "drug_direction", "verdict", "approved", "max_stage",
            "mechanism", "indications", "rationale"]
    path = os.path.join(a.outdir, "ot_target_drugs.tsv")
    with open(path, "w") as fh:
        fh.write("# Open Targets Platform " + release + "\n")
        fh.write("\t".join(cols) + "\n")
        for r in sorted(out_rows, key=lambda r: (r["verdict"], r["target"],
                                                 -stage_rank(r["max_stage"]))):
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    # ---- named drugs -----------------------------------------------------
    named = []
    dq = drug_query(type_fields("Drug"))
    for n in NAMED_DRUGS:
        did, dname = resolve(n, "drug")
        if not did:
            named.append({"query": n, "found": "no"})
            continue
        info = gql(dq, {"id": did})["drug"] or {}
        moas = ((info.get("mechanismsOfAction") or {}).get("rows") or [])
        stage = info.get("maximumClinicalStage", info.get("maximumClinicalTrialPhase"))
        named.append({
            "query": n, "found": "yes", "drug": info.get("name") or dname,
            "chembl_id": did,
            "approved": info.get("isApproved", stage_rank(stage) >= 4),
            "max_stage": stage,
            "drug_type": info.get("drugType"),
            "mechanisms": " | ".join(
                f"{m.get('mechanismOfAction')} [{m.get('actionType')}] -> "
                + ",".join(t.get("approvedSymbol") or t.get("id") or ""
                           for t in m.get("targets") or [])
                for m in moas) or "not annotated",
        })
        time.sleep(0.1)
    ncols = ["query", "found", "drug", "chembl_id", "approved", "max_stage",
             "drug_type", "mechanisms"]
    npath = os.path.join(a.outdir, "ot_named_drugs.tsv")
    with open(npath, "w") as fh:
        fh.write("# Open Targets Platform " + release + "\n")
        fh.write("\t".join(ncols) + "\n")
        for r in named:
            fh.write("\t".join(str(r.get(c, "")) for c in ncols) + "\n")

    # ---- summary ---------------------------------------------------------
    def show(title, rows):
        print(f"\n{title} ({len(rows)})", file=sys.stderr)
        for r in rows[:25]:
            flag = "APPROVED" if r["approved"] else f"{r['max_stage']}"
            print(f"  {r['target']:<9} {r['drug'][:34]:<34} {r['action_type'][:22]:<22} "
                  f"{flag}", file=sys.stderr)
        if len(rows) > 25:
            print(f"  ... {len(rows) - 25} more in {path}", file=sys.stderr)

    by = lambda v: sorted([r for r in out_rows if r["verdict"] == v],
                          key=lambda r: (not r["approved"], -stage_rank(r["max_stage"])))
    show("CANDIDATES: act in the desired direction", by("candidate"))
    show("DO NO HARM: act against the desired direction", by("do_no_harm"))
    print(f"\nreview (ambiguous action): {len(by('review'))} rows", file=sys.stderr)

    print("\nNamed drugs:", file=sys.stderr)
    for r in named:
        if r["found"] == "no":
            print(f"  {r['query']:<28} not in Open Targets", file=sys.stderr)
        else:
            flag = "APPROVED" if r["approved"] else f"{r['max_stage']}"
            print(f"  {r['query']:<28} {flag:<10} {r['mechanisms'][:70]}",
                  file=sys.stderr)
    print(f"\nwrote {path} and {npath}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except GQLError as e:
        sys.exit(f"Open Targets query failed: {e}\n"
                 f"If this names a missing field, run with --introspect Target "
                 f"(or Drug) and send the output.")
