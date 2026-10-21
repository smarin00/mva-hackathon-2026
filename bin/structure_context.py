#!/usr/bin/env python3
"""
Structural context of BubR1 missense positions in the AlphaFold model,
calibrated against MVA alleles with known experimental stability classes.

For each position: AlphaFold confidence (pLDDT), relative solvent
accessibility (RSA, Shrake-Rupley on the full-length model, normalised by the
Tien et al. 2013 theoretical maxima), side-chain contacts, and side-chain
hydrogen-bond partners -- the interactions a substitution would remove.

Why calibrate: "buried" alone is weak evidence. If the experimentally
verified destabilisers (I909T, L1012P, ...) are buried and the stable
variants are exposed, burial discriminates on this protein, and the patient
residue's burial becomes evidence rather than decoration.

Caveats the output keeps in view:
  * The model is a monomer. A residue exposed here may be buried in a
    complex (BubR1 works in the MCC with Bub3, Mad2 and Cdc20).
  * RSA is only interpretable where pLDDT is high (>70); disordered regions
    are flagged, not interpreted.

Setup: uv pip install biopython
Usage:
  structure_context.py --config config/bub1b.toml --uniprot O60566 \\
      --out findings/track2/structure_context.tsv
  structure_context.py --pdb AF-O60566-F1-model.pdb ...   # offline
"""

import argparse
import json
import sys
import urllib.request

# Tien et al. 2013, theoretical maximum ASA (A^2)
MAX_ASA = {"ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167,
           "GLN": 225, "GLU": 223, "GLY": 104, "HIS": 224, "ILE": 197,
           "LEU": 201, "LYS": 236, "MET": 224, "PHE": 240, "PRO": 159,
           "SER": 155, "THR": 172, "TRP": 285, "TYR": 263, "VAL": 174}
THREE = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
         "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
         "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
         "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL"}
BACKBONE = {"N", "CA", "C", "O"}

# Calibration alleles come from the case configuration (config/<gene>.toml).

BURIED_RSA = 0.20          # conventional cut-off
CONTACT_A = 4.5            # heavy-atom contact distance
HBOND_A = 3.5              # donor/acceptor N/O distance


def fetch_model(uniprot, dest):
    api = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot}"
    with urllib.request.urlopen(api, timeout=60) as r:
        entries = json.load(r)
    url = entries[0]["pdbUrl"]
    version = entries[0].get("latestVersion") or entries[0].get("modelCreatedDate")
    urllib.request.urlretrieve(url, dest)
    return url, version


def parse(v):
    return v[0], int(v[1:-1]), v[-1]


def analyse(structure, pos):
    from Bio.PDB import NeighborSearch
    chain = next(structure[0].get_chains())
    res = chain[pos]
    atoms = [a for a in structure[0].get_atoms() if a.element != "H"]
    ns = NeighborSearch(atoms)

    rname = res.get_resname()
    rsa = res.sasa / MAX_ASA[rname]
    plddt = sum(a.get_bfactor() for a in res) / len(res)

    side = [a for a in res if a.get_id() not in BACKBONE] or list(res)
    contacts, hbonds = set(), set()
    for a in side:
        for b in ns.search(a.coord, CONTACT_A):
            other = b.get_parent()
            if other is res or abs(other.id[1] - pos) < 2:
                continue
            tag = f"{other.get_resname()}{other.id[1]}"
            contacts.add(tag)
            if (a.element in "NO" and b.element in "NO"
                    and (a - b) <= HBOND_A):
                hbonds.add(f"{tag}:{b.get_id()}")
    return {
        "resname": rname, "plddt": plddt, "rsa": rsa,
        "n_contacts": len(contacts),
        "contacts": ",".join(sorted(contacts, key=lambda s: int(s[3:]))),
        "hbonds": ",".join(sorted(hbonds)),
    }


def main():
    p = argparse.ArgumentParser()
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--uniprot")
    g.add_argument("--pdb")
    p.add_argument("--config", required=True, help="config/<gene>.toml")
    p.add_argument("--patient", help="default: patient.missense in the config")
    p.add_argument("--out", required=True)
    a = p.parse_args()

    import os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import caseconfig
    cfg = caseconfig.load(a.config)
    CALIBRATION = caseconfig.calibration(cfg)
    a.patient = a.patient or cfg["patient"].get("missense")
    if not a.patient:
        sys.exit("needs a patient missense (--patient or patient.missense in the config)")

    from Bio.PDB import PDBParser
    from Bio.PDB.SASA import ShrakeRupley

    source = a.pdb
    if a.uniprot:
        source = f"AF-{a.uniprot}-F1-model.pdb"
        url, ver = fetch_model(a.uniprot, source)
        print(f"AlphaFold model: {url} ({ver})", file=sys.stderr)

    s = PDBParser(QUIET=True).get_structure("m", source)
    ShrakeRupley().compute(s[0], level="R")
    chain = next(s[0].get_chains())

    rows = []
    for v, cls in [(a.patient, "PATIENT")] + CALIBRATION:
        wt, pos, mt = parse(v)
        if pos not in chain or chain[pos].get_resname() != THREE[wt]:
            got = chain[pos].get_resname() if pos in chain else "absent"
            sys.exit(f"{v}: model has {got} at {pos}, expected {THREE[wt]}. "
                     f"Stopping rather than analysing the wrong residue.")
        r = analyse(s, pos)
        r.update(variant=v, cls=cls)
        rows.append(r)

    cols = ["variant", "cls", "resname", "plddt", "rsa", "n_contacts",
            "hbonds", "contacts"]
    with open(a.out, "w") as fh:
        fh.write(f"# model: {source}\n")
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(f"{r[c]:.3f}" if isinstance(r[c], float)
                               else str(r[c]) for c in cols) + "\n")

    print(f"\n{'variant':<8} {'class':<14} {'pLDDT':>6} {'RSA':>6}  "
          f"{'burial':<9} {'contacts':>8}  side-chain H-bonds", file=sys.stderr)
    for r in rows:
        if r["plddt"] < 70:
            burial = "low-conf"
        else:
            burial = "buried" if r["rsa"] < BURIED_RSA else "exposed"
        print(f"{r['variant']:<8} {r['cls']:<14} {r['plddt']:6.1f} "
              f"{r['rsa']:6.2f}  {burial:<9} {r['n_contacts']:8d}  "
              f"{r['hbonds'] or '-'}", file=sys.stderr)

    conf = [r for r in rows if r["cls"] != "PATIENT" and r["plddt"] >= 70]
    d = [r["rsa"] for r in conf if r["cls"] == "destabilising"]
    st = [r["rsa"] for r in conf if r["cls"] == "stable"]
    if d and st:
        sep = max(d) < min(st)
        print(f"\ncalibration (pLDDT>=70 only): destabilising RSA max "
              f"{max(d):.2f}, stable RSA min {min(st):.2f} -> "
              f"{'SEPARATE' if sep else 'OVERLAP'}", file=sys.stderr)
    print("note: monomer model -- interface residues of the MCC may appear "
          "exposed here.", file=sys.stderr)
    print(f"wrote {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
