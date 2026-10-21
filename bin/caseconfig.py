#!/usr/bin/env python3
"""
Load a case configuration (config/*.toml). Standard library only (tomllib),
so the Makefile can read it before the Python environment exists.

As a module:   cfg = caseconfig.load("config/bub1b.toml")
From make:     caseconfig.py config/bub1b.toml gene            -> BUB1B
               caseconfig.py config/bub1b.toml patient.ptc_pos  -> 40209701
Missing keys print nothing (exit 0), so optional values stay empty in make.
"""
import sys
import tomllib


def load(path):
    with open(path, "rb") as fh:
        cfg = tomllib.load(fh)
    for key in ("gene", "uniprot"):
        if key not in cfg:
            sys.exit(f"{path}: required key '{key}' missing")
    cfg.setdefault("calibration", [])
    cfg.setdefault("pathway", [])
    cfg.setdefault("named_drugs", [])
    cfg.setdefault("stability_features", [])
    cfg.setdefault("patient", {})
    cfg.setdefault("domains", "")
    return cfg


def calibration(cfg):
    """[(variant, class), ...] in the order given."""
    return [(c["variant"], c["class"]) for c in cfg["calibration"]]


def panel(cfg):
    """[(symbol, desired, tier, rationale), ...]"""
    return [(p["symbol"], p["desired"], str(p.get("tier", "1")), p.get("rationale", ""))
            for p in cfg["pathway"]]


def get(cfg, dotted):
    node = cfg
    for part in dotted.split("."):
        if not isinstance(node, dict) or part not in node:
            return None
        node = node[part]
    return node


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: caseconfig.py CONFIG KEY")
    v = get(load(sys.argv[1]), sys.argv[2])
    if v is not None:
        print(v if not isinstance(v, list) else " ".join(map(str, v)))
