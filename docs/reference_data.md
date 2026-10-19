# Reference data and tool versions

Everything here is installed by `make install` (tools) and `make refs`
(references and containers); the exact commands are in
`scripts/install_tools.sh` and `scripts/setup_refs.sh`. `make versions`
records the versions actually present, with container sha256 digests, in
`findings/versions.txt`.

## Reference data

| Resource | Version | Source | Downloaded |
|---|---|---|---|
| Genome | GRCh38 no-alt analysis set + hs38d1 decoys, prebuilt BWA index | NCBI `GCA_000001405.15 ... seqs_for_alignment_pipelines.ucsc_ids` | 2026-10-04 |
| Gene annotation | Ensembl 116 GTF | ftp.ensembl.org release-116 | 2026-10-05 |
| VEP cache | homo_sapiens_vep_116_GRCh38 (indexed) | ftp.ensembl.org release-116 | 2026-10-05 |
| VEP plugins | Ensembl/VEP_plugins, branch `release/116` | GitHub | 2026-10-05 |
| AlphaMissense | AlphaMissense_hg38.tsv.gz (Cheng et al. 2023), tabix-indexed | storage.googleapis.com/dm_alphamissense | 2026-10-05 |
| ClinVar | NCBI GRCh38 VCF, fileDate 2026-10-04, renamed to chr contigs | ftp.ncbi.nlm.nih.gov/pub/clinvar | 2026-10-05 |
| Exomiser data | 2602_hg38 + 2602_phenotype (ClinVar 2026-02-08, HPO 2026-01-08); dbNSFP 4.5a pathogenicity scores | Exomiser discussion #640 | 2026-10-05 |
| ExpansionHunter catalog | GRCh38 catalog bundled with v5.0.0 | GitHub release | 2026-10-05 |
| Candidate panel | 132 genes (MVA, SAC, kinetochore, centromere, cohesion, APC/C, microcephaly and CIN differentials) -> 129 merged regions, 12.20 Mb, 10 kb flanks | built from Ensembl 116 GTF (`bin/build_panel.py`, `docs/panel_genes.txt`) | 2026-10-10 |
| BubR1 sequence | UniProt O60566 (1050 aa) | rest.uniprot.org | 2026-10-19 |
| BubR1 structure | AlphaFold DB model AF-O60566-F1 | alphafold.ebi.ac.uk | 2026-10-19 |
| Drug-target data | Open Targets Platform, API 26.9.0, data 26.09 | api.platform.opentargets.org | 2026-10-19 |
| Protein language model | ESM-2 650M (`esm2_t33_650M_UR50D`) | fair-esm | 2026-10-19 |
| Splice prediction | SpliceAI 1.3.1 (built-in GRCh38 annotation), two patient alleles only | PyPI, separate Python 3.11 env (TensorFlow 2.15 CPU) | 2026-10-19 |

## Tools

| Tool | Version | Installed as |
|---|---|---|
| Nextflow | 26.04.6 | get.nextflow.io |
| Java | OpenJDK 21 | apt |
| Apptainer | from ppa:apptainer/ppa | apt |
| bwa | 0.7.19-r1273 | apt |
| fastp | 1.0.1 | apt |
| samtools, bcftools, tabix | Ubuntu 26.04 packages | apt |
| mosdepth | 0.3.14 | static binary, GitHub release |
| Delly | 2.7.0 (HTSlib 1.21) | static binary, GitHub release |
| ExpansionHunter | 5.0.0 | binary, GitHub release |
| DeepVariant | 1.10.0 | container `google/deepvariant:1.10.0` |
| GATK | 4.7.0.0 | container `broadinstitute/gatk:4.7.0.0` |
| Ensembl VEP | 116.2 | container `ensemblorg/ensembl-vep:release_116.2` |
| Exomiser | CLI 15.1.0 (Java 21) | GitHub release |
| Python | 3.12, packages pinned in `requirements.lock`; PyTorch CPU build | uv |

## Not used

- SpliceAI genome-wide precomputed scores: the BaseSpace share is US-region
  only and the account is EU. Coding and splice-site SNV scores come from
  Exomiser 2602 (dbNSFP 4.5a); the two causal alleles were scored with SpliceAI
  directly (`make spliceai`, `docs/decisions.md`).
- Pangolin: not run (see `docs/decisions.md`).
