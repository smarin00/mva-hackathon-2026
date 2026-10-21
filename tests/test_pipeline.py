"""
Tests for the Track 1 submission rules and the Track 2 scripts, on small
synthetic data in tests/data/ (a made-up 20 kb chr15 carrying a BUB1B-like
gene, a random test protein, invented ClinVar and Exomiser records). No
patient data, no network.

    pytest -q            (from the repository root)
"""
import csv
import json
import os
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "tests", "data")
BIN = os.path.join(ROOT, "bin")
sys.path.insert(0, BIN)


def d(name):
    return os.path.join(DATA, name)


def run(script, *args):
    r = subprocess.run([sys.executable, os.path.join(BIN, script), *map(str, args)],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-2000:]
    return r


def tsv(path):
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


# ---------------------------------------------------------------- config
def test_gene_file_loads_with_defaults():
    import caseconfig
    cfg = caseconfig.load(os.path.join(ROOT, "config/genes/cep57.yaml"))
    assert cfg["gene"] == "CEP57" and cfg["calibration"] == [] and cfg["pathway"] == []
    bub = caseconfig.load(os.path.join(ROOT, "config/genes/bub1b.yaml"))
    assert len(caseconfig.calibration(bub)) == 6 and bub["oncology"]["candidates"]


def test_resolved_case_merges_patient(tmp_path):
    import caseconfig
    rc = tmp_path / "case.yaml"
    rc.write_text("gene_config: config/genes/bub1b.yaml\n"
                  "patient: {chrom: chr15, ptc_pos: 1, missense: N1002K}\n")
    cfg = caseconfig.load(str(rc))
    assert cfg["patient"]["missense"] == "N1002K" and cfg["uniprot"] == "O60566"
    mk = caseconfig.make_vars(str(rc))
    assert "GENE := BUB1B" in mk and "PANEL := yes" in mk and "ONCOLOGY := yes" in mk


# ---------------------------------------------------------------- Track 2 case
def test_resolve_case_derives_consequences(tmp_path):
    case = tmp_path / "params.yaml"
    case.write_text("track2:\n  gene: BUB1B\n  variants:\n"
                    "    - {chrom: chr15, pos: 11861, ref: C, alt: A}\n"
                    "    - {chrom: chr15, pos: 15155, ref: A, alt: C}\n")
    out = tmp_path / "case.yaml"
    run("resolve_case.py", "--case", case, "--gtf", d("test.gtf"),
        "--fasta", d("test_genome.fa"), "--out", out)
    import yaml
    r = yaml.safe_load(out.read_text())
    assert [(v["protein"], v["type"]) for v in r["variants"]] == [("S737*", "ptc"), ("N1002H", "missense")]
    assert r["patient"]["ptc_codon"] == 737 and r["patient"]["missense"] == "N1002H"


def test_resolve_case_rejects_wrong_reference_base(tmp_path):
    case = tmp_path / "params.yaml"
    case.write_text("track2:\n  gene: BUB1B\n  variants:\n    - {chrom: chr15, pos: 11861, ref: G, alt: A}\n")
    r = subprocess.run([sys.executable, os.path.join(BIN, "resolve_case.py"), "--case", str(case),
                        "--gtf", d("test.gtf"), "--fasta", d("test_genome.fa"),
                        "--out", str(tmp_path / "c.yaml")], cwd=ROOT, capture_output=True, text=True)
    assert r.returncode != 0 and "does not match the reference" in r.stderr


# ---------------------------------------------------------------- atlas / ClinVar
def test_variant_atlas_matches_expected(tmp_path):
    run("variant_atlas.py", "--config", "config/genes/bub1b.yaml", "--fasta", d("test_protein.fasta"),
        "--gtf", d("test.gtf"), "--esm-cache", d("test_esm_scan.npz"), "--outdir", tmp_path)
    got = (tmp_path / "atlas_truncating.tsv").read_text()
    assert got == open(d("expected_atlas_truncating.tsv")).read()
    outcomes = [r["outcome"] for r in tsv(tmp_path / "atlas_truncating.tsv")]
    assert outcomes.count("null (NMD)") == 969 and len(outcomes) == 1050
    assert (tmp_path / "atlas_figure.png").stat().st_size > 10_000


def test_clinvar_validation_matches_expected(tmp_path):
    atlas = tmp_path / "atlas"
    run("variant_atlas.py", "--config", "config/genes/bub1b.yaml", "--fasta", d("test_protein.fasta"),
        "--gtf", d("test.gtf"), "--esm-cache", d("test_esm_scan.npz"), "--outdir", atlas)
    out = tmp_path / "cv.tsv"
    run("clinvar_validate.py", "--config", "config/genes/bub1b.yaml", "--clinvar", d("test_clinvar.vcf.gz"),
        "--fasta", d("test_genome.fa"), "--gtf", d("test.gtf"), "--protein", d("test_protein.fasta"),
        "--atlas", atlas / "atlas_missense.tsv", "--out", out)
    assert out.read_text() == open(d("expected_clinvar_validation.tsv")).read()


def test_alphamissense_isoform_guard_and_agreement(tmp_path):
    import pysam
    import alphamissense_compare as ac
    atlas = tmp_path / "atlas"
    run("variant_atlas.py", "--config", "config/genes/bub1b.yaml", "--fasta", d("test_protein.fasta"),
        "--gtf", d("test.gtf"), "--esm-cache", d("test_esm_scan.npz"), "--outdir", atlas)
    chrom, lo, hi = ac.gene_span(d("test.gtf"), "BUB1B")
    rows = tsv(atlas / "atlas_missense.tsv")[::3]
    for name, wt_shift in (("ok", False), ("bad", True)):
        p = tmp_path / f"am_{name}.tsv"
        with open(p, "w") as fh:
            fh.write("# test\n#CHROM\tPOS\tREF\tALT\tgenome\tuniprot_id\ttranscript_id\t"
                     "protein_variant\tam_pathogenicity\tam_class\n")
            for i, r in enumerate(rows):
                s = min(1.0, max(0.0, -float(r["esm_llr"]) / 12))
                wt = "W" if wt_shift else r["wt"]
                fh.write(f"{chrom}\t{lo + i % (hi - lo)}\tA\tC\thg38\tO60566\tT\t"
                         f"{wt}{r['pos']}{r['mt']}\t{s:.3f}\tambiguous\n")
        lines = sorted(open(p).read().splitlines()[2:], key=lambda l: int(l.split("\t")[1]))
        p.write_text("# test\n#CHROM\tPOS\tREF\tALT\tgenome\tuniprot_id\ttranscript_id\t"
                     "protein_variant\tam_pathogenicity\tam_class\n" + "\n".join(lines) + "\n")
        pysam.tabix_index(str(p), seq_col=0, start_col=1, end_col=1, meta_char="#", force=True)
    case = tmp_path / "case.yaml"
    case.write_text("gene_config: config/genes/bub1b.yaml\npatient: {missense: N1002K}\n")
    base = ["--config", case, "--gtf", d("test.gtf"), "--atlas", atlas / "atlas_missense.tsv"]
    run("alphamissense_compare.py", "--am", tmp_path / "am_ok.tsv.gz", *base, "--outdir", tmp_path / "o")
    assert "rho = -1.00" in (tmp_path / "o" / "summary.txt").read_text()
    r = subprocess.run([sys.executable, os.path.join(BIN, "alphamissense_compare.py"),
                        "--am", str(tmp_path / "am_bad.tsv.gz"), *map(str, base), "--outdir",
                        str(tmp_path / "o2")], cwd=ROOT, capture_output=True, text=True)
    assert r.returncode != 0 and "isoform mismatch" in r.stderr


# ---------------------------------------------------------------- related genes
def test_related_genes_finds_mva_series_and_orphanet(tmp_path):
    r = run("related_genes.py", "--gene", "BUB1B", "--g2d", d("hpo_genes_to_disease.txt"),
            "--hpoa", d("hpo_phenotype.hpoa"), "--out", tmp_path / "rel.tsv")
    assert r.stdout.split() == ["BUB1", "BUB3", "CENATAC", "CEP57", "MAD1L1", "TRIP13"]
    assert run("related_genes.py", "--gene", "PAH", "--g2d", d("hpo_genes_to_disease.txt"),
               "--hpoa", d("hpo_phenotype.hpoa")).stdout.strip() == ""


# ---------------------------------------------------------------- Track 1 submission
def test_submission_rules(tmp_path):
    out = tmp_path / "s.csv"
    r = run("make_submission.py", "--ranked", d("test_exomiser.variants.tsv"),
            "--gene-diseases", d("test_genes_to_disease.txt"), "--out", out)
    rows = list(csv.DictReader(open(out)))
    genes = [row["notes"].split(":")[0] for row in rows]
    assert genes[0] == "BUB1B" and rows[0]["chrom_2"] == "chr15"     # pair on one row, first
    assert "HK1" not in genes                                         # too common for AD
    assert "CFTR" in genes                                            # common recessive allele kept
    assert genes[-1] == "MTCH2"                                       # no disease: demoted to last
    assert "too common for AD, kept as AR" in r.stderr
    epcr = [float(row["epcr"]) for row in rows]
    assert epcr == sorted(epcr, reverse=True) and len(set(epcr)) == len(epcr)


def test_submission_curation_applies(tmp_path):
    cur = tmp_path / "cur.tsv"
    cur.write_text("MTCH2\texclude\tno disease\nLZTR1\tsecondary\tSecondary finding\n")
    out = tmp_path / "s.csv"
    run("make_submission.py", "--ranked", d("test_exomiser.variants.tsv"), "--curation", cur, "--out", out)
    rows = list(csv.DictReader(open(out)))
    assert not any(r["notes"].startswith("MTCH2") for r in rows)
    assert rows[-1]["finding_type"] == "secondary" and rows[-1]["notes"].startswith("LZTR1")


# ---------------------------------------------------------------- oncology rules
def test_oncology_matching_on_real_label_sentences():
    import oncology_compatibility as o
    T = o.TOXICITY

    def f(text, pattern, tox=False):
        return bool(o.find({"x": [text]}, ["x"], pattern, toxicity=tox))
    assert not f("Arimoclomol is not an inhibitor or inducer of CYP1A2, CYP2C9 or CYP3A4/5.",
                 o.cyp_rx("arimoclomol", "modulator"))
    assert not f("7.1 Effect of Other Drugs on Dasatinib Strong CYP3A4 Inhibitors The coadministration "
                 "with strong CYP3A inhibitors may increase dasatinib concentrations.",
                 o.cyp_rx("dasatinib", "modulator"))
    assert f("Dasatinib is a time-dependent inhibitor of CYP3A4.", o.cyp_rx("dasatinib", "modulator"))
    assert f("Vincristine sulfate is a substrate for cytochrome P450 3A isozymes (CYP3A).",
             o.cyp_rx("vincristine", "substrate"))
    assert not f("The risk is increased in patients with diabetes, renal failure, or hypothyroidism;",
                 T["kidney"], True)
    assert not f("not administered to patients with serum bilirubin > 2.0 mg/dL, or transaminases > 3 "
                 "times ULN.", T["liver"], True)
    assert not f("Pulmonary veno-occlusive disease leading to respiratory failure may occur.",
                 T["liver"], True)
    assert f("Veno-occlusive Disease: Can cause severe or fatal VOD.", T["liver"], True)


def test_oncology_label_selection():
    import oncology_compatibility as o

    def lab(g, pt="HUMAN PRESCRIPTION DRUG", **kw):
        return {"openfda": {"generic_name": [g], "product_type": [pt]}, "effective_time": "2025",
                "warnings": ["x"], **{k: [v] for k, v in kw.items()}}
    assert o.pick([lab("DOXORUBICIN HYDROCHLORIDE", boxed_warning="liposome injection"),
                   lab("DOXORUBICIN HYDROCHLORIDE")], "doxorubicin").get("boxed_warning") is None
    assert o.pick([lab("LOVASTATIN AND NIACIN"), lab("NIACIN")], "niacin")["openfda"]["generic_name"] == ["NIACIN"]
    assert o.pick([lab("QUERCETIN", pt="HUMAN OTC DRUG")], "quercetin") == {}


def test_oncology_end_to_end_offline(tmp_path):
    labels = tmp_path / "labels"
    labels.mkdir()

    def write(drug, **fields):
        lab = {k: [v] for k, v in fields.items()}
        lab.update(openfda={"generic_name": [drug.upper()], "product_type": ["HUMAN PRESCRIPTION DRUG"]},
                   set_id="s", version="1", effective_time="20250101")
        (labels / (drug.replace(" ", "_") + ".json")).write_text(json.dumps(lab))
    write("vincristine", warnings_and_cautions="Myelosuppression may occur.",
          clinical_pharmacology="Vincristine is a substrate of CYP3A4.")
    for drugname in ("dactinomycin", "cyclophosphamide", "irinotecan", "doxorubicin"):
        write(drugname, warnings_and_cautions="Severe myelosuppression.")
    write("dasatinib", warnings_and_cautions="Myelosuppression may occur.",
          clinical_pharmacology="Dasatinib is a time-dependent inhibitor of CYP3A4.")
    for drugname in ("niacin", "arimoclomol", "valproic acid"):
        write(drugname, warnings_and_cautions="Hypersensitivity.")
    for drugname in ("quercetin", "nicotinamide riboside"):
        (labels / (drugname.replace(" ", "_") + ".json")).write_text("{}")
    run("oncology_compatibility.py", "--config", "config/genes/bub1b.yaml", "--outdir", tmp_path / "o",
        "--labels", labels, "--offline")
    rows = tsv(tmp_path / "o" / "compatibility.tsv")
    flags = {(r["candidate"], r["chemotherapy"], r["check"]) for r in rows}
    assert ("dasatinib", "vincristine", "cyp3a4") in flags
    assert not any(r["candidate"] == "arimoclomol" for r in rows)


# ---------------------------------------------------------------- Makefile
@pytest.mark.parametrize("target", ["help", "-n clean"])
def test_makefile_parses_without_case(target):
    r = subprocess.run(["make", *target.split(), "CASE=/nonexistent.yaml"], cwd=ROOT,
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
