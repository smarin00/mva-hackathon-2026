# findings/

Shareable results: everything here is a *finding* in the sense of the
hackathon's data terms ("a handful of named variants is a finding; a
genome-wide table of genotypes is the dataset in a different format").

| Path | Contents | Made by |
|---|---|---|
| `track1/` | Submission CSVs, Nextflow run reports | `make findings` |
| `aneuploidy/` | Karyotype figure and per-chromosome summary tables (three GC models) | `make findings` |
| `track2/case.yaml` | Resolved Track 2 case: gene, transcript, the two variants and their consequences | `make case` |
| `track2/nmd_*.txt` | NMD geometry of the patient's stop-gain | `make nmd` |
| `track2/esm_scores.tsv` | ESM-2 score of the patient missense and six calibration alleles | `make esm` |
| `track2/structure_context.tsv` | AlphaFold structural context, same seven residues | `make structure` |
| `track2/spliceai_two_alleles.tsv` | SpliceAI delta scores for the two causal alleles | `make spliceai` |
| `track2/opentargets/` | Direction-aware drug screen and named-drug lookup | `make opentargets` |
| `track2/clinvar_validation.tsv` | ClinVar BUB1B missense variants placed in the atlas (independent check) | `make clinvar` |
| `track2/atlas/` | All possible BubR1 missense and premature-stop variants; ESM cache | `make atlas` |
| `track2/alphamissense/` | AlphaMissense against ESM-2: calibration alleles, agreement, ClinVar, VUS triage | `make alphamissense` |
| `track2/oncology/` | Candidates against chemotherapy regimens, with the FDA label sentences and label versions | `make oncology` |
| `track2/related_genes.tsv` | Genes linked to MVA in HPO (OMIM series, Orphanet group) | `make related-genes` |
| `genes/<GENE>/` | Constraint and NMD atlas for each related gene (no patient data) | `make related-genes` |
| `versions.txt` | Tool versions and container digests | `make versions` |

Track 2 outputs are computed from public data plus the two named variants
(UniProt, AlphaFold DB, Open Targets, Ensembl, ClinVar, AlphaMissense, HPO,
openFDA). The figures in `docs/figures/` are drawn from these tables by
`scripts/report_figures.py`.

Before committing anything new here, check it is a summary or a named-variant
list and not a per-variant genome-wide table.
