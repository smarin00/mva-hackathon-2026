# Data deletion plan

Deadline: 30 days after hackathon close (close 2026-10-24) -> delete by **2026-11-23**
Confirm by email: MVAHackathon2026@synapse.org (address agreed on the dataset access form); cc RarediseaserealkidMVAhackathon2026@synapse.org (address on the Rules tab); cc cnb@sagebase.org (address given in community notes)
Rule (organiser, Discussion #2): delete anything that carries the child's genome at
genotype scale; keep findings. "A handful of named variants is a finding; a genome-wide
table of genotypes is the dataset in a different format."

## DELETE — server (Hetzner CPX62, project mva-hackathon)
- [ ] ~/mva/data/          FASTQ, provided VCF + .tbi, phenotype .docx, README
- [ ] ~/mva/derived/       vendor.chr.vcf.gz + .tbi, phenotype.txt (if created)
- [ ] ~/mva/testdata/      2 M-read FASTQ subset of lane 1 (smoke test): PATIENT READS
- [ ] ~/mva/work/          all Nextflow intermediates (BAM/CRAM, VCFs, caches)
- [ ] ~/mva/results/alignment/   CRAM + index
- [ ] ~/mva/results/variants/    normalised / DeepVariant / mosaic / SV / repeat VCFs
- [ ] ~/mva/results/annotation/  VEP-annotated VCFs
- [ ] ~/mva/results/ranking/*.variants.tsv, *.html, *.json   (variant-level Exomiser output)
- [ ] ~/mva/results/phasing/, results/qc/ (per-sample QC)
- [ ] ~/mva/results_full/   full read-level run: CRAM, all VCFs, phasing, QC, mosdepth
- [ ] ~/mva/results_vcf/, ~/mva/results_smoke/   other run outputs
- [ ] Check nothing is left:  make deletion-check   (should print nothing)
- [ ] ~/.cache/huggingface/  (xet download cache)
- [ ] Then: delete the server instance itself + any snapshots/backups in Hetzner console

## DELETE — elsewhere
- [ ] Laptop: downloaded submission CSVs and any other files from the server (Downloads folder)
- [ ] GitHub: confirm no data files ever committed:
      git log --all --name-only --format='' | grep -E '\.(vcf|bam|cram|fastq|tbi|csi)' || echo clean
- [ ] claude.ai: delete any chats where variant blocks or phenotype text were pasted

## KEEP (findings; must contain no items from the delete list)
- Submission CSVs (a handful of named variants)
- HPO term list
- Exomiser *.genes.tsv (gene rankings), mechanism and pathway analyses, drug candidates
- Aneuploidy figure + per-chromosome summary table (aggregate, not genotypes)
- Code, docs, report, pitch, leaderboard entry
- results/pipeline_info/ (run reports; no genotypes)
- Everything in findings/ (collected by `make findings`, reviewed before each commit):
  Track 2 tables (ESM scores, structural context, Open Targets screen, variant atlas)

## Attestation (planned; send after deletion is complete)
"I will delete all MVA Hackathon data and derived genotype-level files from all systems
under my control (cloud server and its storage, laptop, code repository) no later than
23 November 2026, 30 days after the hackathon closed on 24 October 2026.
Third-party LLM used: Anthropic Claude via claude.ai (Pro - Opus 5.5), model training disabled,
no outputs rated; no genome-scale data was shared with it."

