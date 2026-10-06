# MVA Hackathon 2026 — singleton WGS analysis pipeline

Reproducible Nextflow pipeline for the
[Rare Disease, Real Kid: MVA Hackathon 2026](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026)
(Sage Bionetworks, MVA Society, Hugging Face, BEACON). It takes one patient's whole-genome
sequencing from raw reads or a variant call file to a ranked list of candidate causal variants,
with read-level analyses aimed at Mosaic Variegated Aneuploidy (MVA).

> **No patient data is stored in this repository.** The dataset is access-controlled and is
> used only on a dedicated server under the hackathon's data-use terms. See
> [Data handling](#data-handling).

## What the pipeline does

```
FASTQ (4 lanes)                                provided VCF
   │ fastp → bwa mem → sort / markdup → CRAM       │
   ▼                                               │
DeepVariant re-calling ──────────────┐             │
Mutect2 (mosaic, candidate panel)    │             │
Delly (structural variants)          │             │
ExpansionHunter (repeat expansions)  │             │
CNVpytor + mosdepth (aneuploidy)     │             │
                                     ▼             ▼
                          bcftools normalisation (PASS, split, left-align)
                                     │
                  ┌──────────────────┴───────────────────┐
                  ▼                                      ▼
        Exomiser 15.1 (HPO-driven ranking)     VEP 116 (+ AlphaMissense, ClinVar)
                  │
                  ▼
     WhatsHap read-backed phasing (candidate panel)
                  │
                  ▼
     make_submission.py → validated Track 1 CSV
```

Two run modes:

| Mode | Input | Use |
|---|---|---|
| `--mode vcf` | Provided VCF | Fast baseline: normalisation, Exomiser, VEP, submission (~2 h) |
| `--mode full` (default) | Raw FASTQ | Re-alignment, re-calling, mosaic / SV / repeat / aneuploidy analyses |

## Submission logic

`bin/make_submission.py` turns Exomiser's ranking into the Track 1 format. It doesn't just
copy the top rows:

- **Compound heterozygotes on one row.** Only a row holding both alleles can score a full
  match, so recessive pairs are submitted together, with the single alleles below as hedges.
- **Biologically valid pairs only.** Variants are paired only if Exomiser scored them under a
  recessive mode, and only if they are at least 50 bp apart; closer variants are usually one
  complex event on one haplotype.
- **Exclusions:** ClinVar-benign variants, hypervariable loci (HLA, KIR, MUC, olfactory
  receptors, Ig/TCR), lone heterozygous hits in recessive-only genes, and genes below a
  minimum Exomiser score.
- **Documented curation:** gene-level decisions with reasons, in `docs/curation.tsv`.
- **Self-validation:** the output is checked against the scorer's rules, and every variant is
  confirmed to exist exactly as written in the source VCF.

## Requirements

- Linux, Java ≥ 17, [Nextflow](https://www.nextflow.io/) 26.04, [Apptainer](https://apptainer.org/)
- [uv](https://docs.astral.sh/uv/) for the Python environment (`requirements.lock`)
- Native tools: samtools, bcftools, tabix, bwa, fastp, mosdepth, delly, ExpansionHunter
- Containers (pulled automatically or pinned locally): DeepVariant, GATK, Ensembl VEP 116
- Hardware used: 16 vCPU, 32 GB RAM, 600 GB NVMe (Hetzner CPX62)

Exact versions: `docs/reference_data.md` and `results/pipeline_info/versions.txt`.

## Reference data

| Resource | Version |
|---|---|
| Genome | GRCh38 no-alt analysis set + hs38d1 decoys (NCBI), with BWA index |
| Exomiser data | 2602 (hg38 + phenotype) |
| VEP cache / plugins | 116 / VEP_plugins release/116 |
| AlphaMissense | hg38 precomputed scores |
| ClinVar | NCBI GRCh38 VCF, renamed to chr contigs |
| ExpansionHunter catalog | GRCh38 catalog bundled with v5.0.0 |

Download sources and dates: `docs/reference_data.md`.

## Usage

```bash
# Python environment
uv venv --python 3.12 && source .venv/bin/activate
uv pip sync requirements.lock

# Configure (paths, HPO terms, optional resources)
cp params.example.yaml params.yaml      # params.yaml is gitignored

# Validate
nextflow lint main.nf nextflow.config

# Run
nextflow run main.nf -params-file params.yaml --mode vcf    # baseline from the VCF
nextflow run main.nf -params-file params.yaml               # full run from FASTQ
```

Add `-resume` to reuse finished steps. Note that `-resume` does not detect edits to scripts in
`bin/`: change the process command in `main.nf`, or clear that task's work directory, to force
a rerun.

## Repository layout

```
main.nf                 pipeline
nextflow.config         parameters, resources, containers
params.example.yaml     configuration template
bin/                    make_submission.py, submission_table.py, read_phenotype.py
docs/
  reference_data.md     tools and reference data, with versions and sources
  decisions.md          analysis decisions and rationale
  curation.tsv          gene-level curation of submitted candidates
  llm_usage.md          AI-assistant use declaration
  DELETION.md           data deletion plan
  *_report.md           Track 1 / Track 2 reports
requirements.in/.lock   Python environment
```

## Outputs

| Folder | Content | In repo? |
|---|---|---|
| `results/submission/` | Track 1 CSV | No |
| `results/ranking/` | Exomiser gene and variant rankings | No |
| `results/annotation/`, `results/variants/`, `results/alignment/` | VCFs, CRAM | No (genotype-level) |
| `results/pipeline_info/` | Nextflow report, timeline, trace, versions | Yes, as evidence |

## Data handling

- Data obtained under the hackathon's gated-access terms (WCG IRB protocol #20252010); never
  redistributed.
- Analysis on a dedicated, key-only cloud server. Patient data never enters this repository;
  `.gitignore` excludes all data formats, and the git history has been checked.
- All data and genotype-level derived files are deleted within 30 days of the hackathon close,
  following the organisers' guidance; plan and checklist in `docs/DELETION.md`.
- AI-assistant use is declared in `docs/llm_usage.md` and in each report.

## Acknowledgement

This work was made possible through the Hackathon, organized by Sage Bionetworks in
partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and
Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are
deeply grateful to the child and their family who generously contributed their data and their
story to advance research into this rare disease. We acknowledge their trust in making this
Hackathon possible.

## Licence

Code and documentation are released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/),
in line with the hackathon's submission terms. The underlying dataset is not covered by this
licence and is not distributed here.

## Author

Santiago Marin Martinez ([@smarin00](https://github.com/smarin00))
