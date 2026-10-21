# MVA Hackathon 2026: from a child's genome to drug hypotheses

Code for the [Rare Disease, Real Kid: MVA Hackathon 2026](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026)
(Sage Bionetworks, MVA Society, Hugging Face, BEACON).

- **Track 1** finds the genetic cause: it takes one patient's whole-genome
  sequencing, from raw reads or from the provided variant file, to a ranked,
  validated list of candidate variants.
- **Track 2** explains the mechanism and proposes approved medicines: it works
  out what the variants do to the protein and screens drug databases for
  medicines that act in the right direction, and for medicines to avoid.

> **No patient data is in this repository.** The dataset is access-controlled
> and was analysed only on a dedicated server under the hackathon's data-use
> terms. See [Data handling](#data-handling).

## Results at a glance

| | Result | Details |
|---|---|---|
| Track 1 | Both causal BUB1B variants ranked first as a compound-heterozygous pair (rank points 100/100, F-max 1.000) | `docs/smarin00_track1_report.pdf` |
| Mechanism | p.Leu737\* is degraded by nonsense-mediated decay (null allele); p.Asn1002Lys most likely destabilises BubR1, so the disease is a matter of protein *quantity* | `docs/track2_conclusions.md` |
| Candidates | Niacin (NAD⁺ → SIRT2 → BubR1 stability), arimoclomol (chaperone support), dasatinib + quercetin (senescence); do-no-harm list incl. HSP90, NAMPT, Aurora B and MPS1 inhibitors | `docs/track2_candidates.md` |
| Resource | Atlas of every possible BubR1 missense and premature-stop variant, each with predicted mechanism and therapy direction | `findings/track2/atlas/` |

These are hypotheses for experimental follow-up, not evidence that any
medicine works.

## Quick start

On a fresh Ubuntu server (16 CPU cores, 32 GB RAM, ~700 GB disk):

```bash
git clone https://github.com/smarin00/mva-hackathon-2026 && cd mva-hackathon-2026
make install     # system tools                        (~15 min, needs sudo)
make refs        # reference genome and databases      (1-3 h, ~100 GB)
make setup       # Python environment                  (~5 min)
make data        # patient dataset (approved access only, see below)
make inputs      # prepare inputs; creates params.yaml -> add the HPO terms
make all         # Track 1, then Track 2 on its top candidate, then the
                 # atlases for related genes                (~4 h VCF route)
```

`make` on its own prints the menu; each part also runs on its own
(`make track1-vcf`, `make track2`, `make mva-genes`). Every step can be
re-run safely: finished work is reused, not repeated.

**One file describes the case.** `params.yaml` holds everything about one
patient: data paths, HPO terms, the Track 1 route (`mode: vcf` or `full`) and
the Track 2 variants. Left at `auto`, Track 2 takes the gene and variants from
Track 1's top-ranked candidate and works out each variant's consequence (stop,
missense) from the reference genome. Knowledge about a gene (domains,
calibration alleles, drug-screen panel) lives in `config/genes/<gene>.yaml`;
for a gene without one, a minimal file is generated from UniProt.

**Track 2 runs without the patient data.** It needs only public resources (the
Ensembl annotation, the protein sequence and structure, drug databases) and
the two variants, which are in `config/cases/mva1_proband.yaml`. On any Linux
machine with [uv](https://docs.astral.sh/uv/):
`make setup gtf track2 CASE=config/cases/mva1_proband.yaml REFS=$HOME/refs`.

## What you need

| | |
|---|---|
| Computer | Linux (tested on Ubuntu 26.04), 16 cores, 32 GB RAM. No GPU needed. |
| Disk | ~100 GB for references; +~500 GB for a full read-level run |
| Time | VCF route ~2 h; full read-level route ~15 h; Track 2 ~1.5 h |
| Patient data | Approved access to `SageBio/mva-hackathon-2026-data` on Hugging Face, then `hf auth login` |
| Phenotype | HPO terms for the patient in `params.yaml` (template: `params.example.yaml`) |

## The steps, explained

**1. Tools and references** (`make install`, `make refs`). Installs aligners,
variant callers and annotation tools at fixed versions, and downloads the
reference genome, Ensembl annotation, ClinVar, AlphaMissense and the Exomiser
databases. Every version and URL is listed in `docs/reference_data.md`.

**2. Track 1: find the cause.** Two routes, both ending in a validated
submission file:

```
FASTQ (4 lanes)                                provided VCF
   | fastp -> bwa mem -> sort / markdup -> CRAM    |
   v                                               |
DeepVariant re-calling (primary chromosomes) --+   |
Mutect2 (mosaic, candidate panel)              |   |
Delly (structural variants)                    |   |
ExpansionHunter (repeat expansions)            |   |
mosdepth + B-allele frequency (aneuploidy)     |   |
                                               v   v
                   bcftools: PASS, primary contigs, split, left-align
                                     |
                 +-------------------+--------------------+
                 v                                        v
       Exomiser 15.1 (HPO-driven ranking)      VEP 116 (+ AlphaMissense, ClinVar)
                 |
                 v
       WhatsHap read-backed phasing (candidate panel)
                 |
                 v
       make_submission.py -> validated Track 1 CSV
```

- `make track1-vcf`: fast route from the provided VCF.
- `make track1`: full route from raw reads; adds re-calling, mosaic, structural,
  repeat and aneuploidy analyses.
- `make smoke`: runs the full route on a 2 M-read subset of lane 1, to check
  an installation before committing to the ~15 h full run.

How the submission file is built (`bin/make_submission.py`):

- Compound heterozygotes go on one row, since only a row holding both alleles
  scores a full match; the single alleles follow as hedges.
- Variants are paired only if Exomiser scored them under a recessive mode and
  they are at least 50 bp apart.
- Excluded: ClinVar-benign variants, hypervariable loci (HLA, KIR, MUC,
  olfactory receptors, Ig/TCR), lone heterozygous hits in recessive-only genes,
  and genes below a minimum Exomiser score. Gene-level curation with reasons is
  in `docs/curation.tsv`.
- The output is checked against the scorer's rules, and every variant is
  confirmed to exist exactly as written in the source VCF.

**3. Track 2: mechanism and medicines** (`make track2`).

| Step | Script | Question it answers |
|---|---|---|
| `make nmd` | `bin/nmd_check.py` | Is the truncating allele's mRNA degraded, or does it make protein? |
| `make esm` | `bin/esm_variant_score.py` | Does the missense look like the known MVA stability mutants? (ESM-2, calibrated on six experimentally characterised alleles) |
| `make structure` | `bin/structure_context.py` | Where does the residue sit in the AlphaFold structure, and which contacts does the change break? |
| `make opentargets` | `bin/opentargets_candidates.py` | Which drugs push the pathway in the desired direction, and which push it the wrong way? |
| `make case` | `bin/resolve_case.py` | Which gene and variants, and what does each do to the protein (stop, missense)? |
| `make atlas` | `bin/variant_atlas.py` | The same questions for *every* possible BubR1 variant, for other families |
| `make clinvar` | `bin/clinvar_validate.py` | Does the atlas agree with ClinVar variants it never saw? |
| `make mva-genes` | `bin/variant_atlas.py` | Constraint and NMD maps for the other MVA genes, CEP57 (MVA2) and TRIP13 (MVA3) |
| `make clinvar` | `bin/clinvar_validate.py` | Does the atlas agree with ClinVar's pathogenic and benign BUB1B variants, which it never saw? |

`make aneuploidy-sensitivity` re-runs the aneuploidy analysis under three GC
correction models; chromosomes whose depth shift follows GC content are
flagged as unconfirmed rather than called.

## What comes out

| Folder | Contents | In the repository? |
|---|---|---|
| `results_vcf/`, `results_full/` | Alignments, VCFs, rankings, QC | **Never**: patient data |
| `findings/` | Submission files, figures, per-chromosome summaries, Track 2 tables | Yes: collected by `make findings`, reviewed before commit |
| `docs/` | Reports, decisions, curation, data handling | Yes |

## Reusing it

- **Another patient:** a new case file (copy `params.example.yaml`: data
  paths, HPO terms), then `make all CASE=other.yaml`. Track 2 follows Track 1's
  top candidate, in whatever gene it is.
- **Another BubR1 variant:** look it up in `findings/track2/atlas/`, or list
  it under `track2.variants` in a case file and run `make track2 CASE=...`.
- **Another gene:** `make gene-atlas G=SYMBOL` writes `config/genes/<gene>.yaml`
  from UniProt and builds the atlas. Adding knowledge to that file switches
  on more analyses: a curated domain map replaces UniProt's annotations,
  calibration alleles turn constraint into mechanism calls, and a `pathway`
  panel switches on the drug screen. That curation, not the code, is the
  real work for a new gene.

## Repository layout

```
Makefile                single entry point (run `make`)
main.nf                 Track 1 pipeline (Nextflow)
nextflow.config         resources, containers, error handling
params.example.yaml     case file template -> params.yaml (not committed)
bin/                    all analysis scripts (Track 1 and Track 2)
config/genes/           gene knowledge: domains, calibration, drug panel (one YAML each)
config/cases/           committed Track 2 case for reproducing without patient data
scripts/                install_tools.sh, setup_refs.sh, prepare_inputs.sh
resources/              chr_map.txt (contig renaming); downloaded files land here
findings/               shareable results
docs/
  reference_data.md     tools and references, versions and sources
  decisions.md          analysis decisions and rationale
  curation.tsv          gene-level curation of submitted candidates
  bub1b_domains.tsv     BubR1 functional map with literature sources
  track2_*.md           mechanism conclusions and drug candidates
  llm_usage.md          AI-assistant use declaration
  DELETION.md           data deletion plan
requirements.in/.lock   Python environment
```

## Data handling

- Data obtained under the hackathon's gated-access terms (WCG IRB protocol
  #20252010); never redistributed.
- Analysis on a dedicated, key-only cloud server. Patient data never enters
  this repository: `.gitignore` excludes every data format and all results
  folders, and the git history has been checked.
- All data and genotype-level derived files are deleted by 23 November 2026;
  plan and checklist in `docs/DELETION.md`, and `make deletion-check` lists
  anything still on disk.
- AI-assistant use is declared in `docs/llm_usage.md`.

## Acknowledgement

This work was made possible through the Hackathon, organized by Sage
Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The
Benchmarking, Evaluation, and Assessment Consortium for Science), with prize
sponsorship from AWS and Anthropic. We are deeply grateful to the child and
their family who generously contributed their data and their story to advance
research into this rare disease. We acknowledge their trust in making this
Hackathon possible.

## Licence

Code and documentation are released under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), in line with the
hackathon's submission terms. The underlying dataset is not covered by this
licence and is not distributed here.

## Author

Santiago Marin Martinez ([@smarin00](https://github.com/smarin00))
