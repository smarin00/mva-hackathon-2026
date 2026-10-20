# findings/

Shareable results: everything here is a *finding* in the sense of the
hackathon's data terms ("a handful of named variants is a finding; a
genome-wide table of genotypes is the dataset in a different format").

| Path | Contents | Made by |
|---|---|---|
| `track1/` | Submission CSVs, Nextflow run reports | `make findings` |
| `aneuploidy/` | Karyotype figure and per-chromosome summary tables (three GC models) | `make findings` |
| `track2/nmd_*.txt` | NMD geometry of the patient's stop-gain | `make nmd` |
| `track2/esm_scores.tsv` | ESM-2 score of the patient missense and six calibration alleles | `make esm` |
| `track2/structure_context.tsv` | AlphaFold structural context, same seven residues | `make structure` |
| `track2/spliceai_two_alleles.tsv` | SpliceAI delta scores for the two causal alleles | `make spliceai` |
| `track2/opentargets/` | Direction-aware drug screen and named-drug lookup | `make opentargets` |
| `track2/clinvar_validation.tsv` | ClinVar BUB1B missense variants placed in the atlas (independent check) | `make clinvar` |
| `track2/atlas/` | All possible BubR1 missense and premature-stop variants; ESM cache | `make atlas` |
| `versions.txt` | Tool versions and container digests | `make versions` |

Track 2 outputs other than the NMD file are computed from public data only
(UniProt, AlphaFold DB, Open Targets, Ensembl).

Before committing anything new here, check it is a summary or a named-variant
list and not a per-variant genome-wide table.
