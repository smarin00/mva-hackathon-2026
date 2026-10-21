# MVA Hackathon 2026 -- single entry point.  Run `make` for the menu.
#
# One case file (params.yaml, from params.example.yaml) describes a patient
# for both tracks; gene knowledge lives in config/genes/<gene>.yaml and is
# generated from UniProt when missing. After the one-time setup,
#     make all
# runs Track 1, hands the top-ranked candidate to Track 2, and builds the
# atlases for genes that cause related diseases (looked up in HPO).

# A recipe that fails deletes its half-written target, so a rerun redoes it.
.DELETE_ON_ERROR:

# ---- case -----------------------------------------------------------------
CASE ?= params.yaml
PYC   = $(if $(wildcard .venv/bin/python),.venv/bin/python,python3)
CFG   = $(PYC) bin/caseconfig.py
# value from the case file (evaluated only where used)
case  = $(shell $(CFG) $(CASE) $(1) 2>/dev/null)
ROUTE   = $(or $(call case,mode),full)
# auto | list of symbols | empty (mva_genes: older name of the same key)
RELATED = $(or $(call case,related_genes),$(call case,mva_genes))
# GENE, UNIPROT, DOMAINS, PANEL, PATIENT_*: from the resolved case, written
# by `make case` as build/case.mk and included below once it exists.

# ---- locations ------------------------------------------------------------
REFS       ?= /opt/refs
CONTAINERS ?= /opt/containers
GTF        ?= $(REFS)/Homo_sapiens.GRCh38.116.gtf.gz
FASTA      ?= $(REFS)/GRCh38.fa
OUT_VCF    ?= results_vcf
OUT_FULL   ?= results_full
FINDINGS   ?= findings
PY         ?= .venv/bin/python
NXF        ?= nextflow run main.nf -params-file $(CASE)

T2    := $(FINDINGS)/track2
# internal state (not results): generated make variables, gene lists, temp files
BUILD := build
RCASE := $(T2)/case.yaml
-include $(wildcard $(BUILD)/case.mk)
OUT    = $(if $(filter vcf,$(ROUTE)),$(OUT_VCF),$(OUT_FULL))
SUBMISSION = $(OUT)/submission/submission.csv


.PHONY: help all install refs gtf setup data inputs inputs-run lint smoke track1-vcf track1 \
	    case gene-config gene-atlas related-genes related-genes-run mva-genes hpo \
	    track2 track2-steps \
	    nmd esm structure opentargets atlas spliceai clinvar aneuploidy-sensitivity \
	    findings versions deletion-check clean

help:
	@echo "MVA Hackathon 2026 pipeline. Typical order:"
	@echo ""
	@echo "  make install      system tools (sudo): Nextflow, Apptainer, bwa, samtools, ..."
	@echo "  make refs         reference genome, VEP, Exomiser, ClinVar, containers (~100 GB)"
	@echo "  make gtf          only the Ensembl annotation (enough for Track 2, 60 MB)"
	@echo "  make setup        Python environment (.venv)"
	@echo "  make data         download the patient dataset (needs approved access)"
	@echo "  make inputs       prepare VCF and gene panel; creates params.yaml"
	@echo ""
	@echo "  make all          everything below, as set in params.yaml (CASE=...)"
	@echo ""
	@echo "  make track1-vcf   Track 1 from the provided VCF       (~2 h)  -> $(OUT_VCF)/"
	@echo "  make track1       Track 1 from raw reads (FASTQ)      (~15 h) -> $(OUT_FULL)/"
	@echo "  make track2       Track 2 mechanism + drug analyses   (~1.5 h) -> $(T2)/"
	@echo "                    variants from params.yaml, or Track 1's top candidate"
	@echo "  make related-genes  atlas for genes causing related diseases (HPO;"
	@echo "                    related_genes: auto in params.yaml, or a list)"
	@echo "  make gene-atlas G=CEP57   atlas for any gene -> $(FINDINGS)/genes/"
	@echo ""
	@echo "  make smoke        quick pipeline test on a 2 M-read subset"
	@echo "  make findings     copy shareable results into $(FINDINGS)/ (safe to commit)"
	@echo "  make versions     record tool and container versions"
	@echo "  make deletion-check   list patient-derived files still on disk"

# ---- setup ------------------------------------------------------------------
install:
	bash scripts/install_tools.sh

refs:
	REFS=$(REFS) CONTAINERS=$(CONTAINERS) bash scripts/setup_refs.sh
	$(MAKE) --no-print-directory hpo

# HPO gene-disease annotations (OMIM + Orphanet), pinned release
HPO_REL ?= v2026-09-01
HPO      = $(REFS)/hpo/$(HPO_REL)
hpo: $(HPO)/genes_to_disease.txt $(HPO)/phenotype.hpoa
$(HPO)/%:
	mkdir -p $(@D)
	wget -q -O $@ https://github.com/obophenotype/human-phenotype-ontology/releases/download/$(HPO_REL)/$*

# Track 2 needs only the gene annotation, not the full reference set
gtf: $(GTF)
$(GTF):
	mkdir -p $(dir $(GTF))
	wget -q -O $@ https://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz

setup: .venv/.done
.venv/.done: requirements.lock
	uv venv --python 3.12
	uv pip sync requirements.lock
	uv pip install torch==2.14.1 --index-url https://download.pytorch.org/whl/cpu
	touch $@

data:
	DATASET=$(call case,dataset) bash scripts/prepare_inputs.sh data

inputs:
	@test -f $(CASE) || { cp params.example.yaml $(CASE); \
	    echo "Created $(CASE) from the template: add the HPO terms before running."; }
	$(MAKE) --no-print-directory inputs-run
inputs-run:
	REFS=$(REFS) PANEL_GENES=$(or $(call case,panel_genes),docs/panel_genes.txt) \
	PANEL_BED=$(or $(call case,panel_bed),resources/mva_panel.bed) \
	    bash scripts/prepare_inputs.sh vcf panel

# ---- everything -------------------------------------------------------------
all:
	$(MAKE) --no-print-directory $(if $(filter vcf,$(ROUTE)),track1-vcf,track1 aneuploidy-sensitivity)
	$(MAKE) --no-print-directory track2
	$(MAKE) --no-print-directory related-genes
	$(MAKE) --no-print-directory findings versions

lint:
	nextflow lint main.nf nextflow.config

# ---- Track 1 ----------------------------------------------------------------
smoke:
	bash scripts/prepare_inputs.sh smoke
	$(NXF) --fastq_dir testdata --outdir results_smoke --sample_id TESTRUN -resume

track1-vcf:
	$(NXF) --mode vcf --outdir $(OUT_VCF) -resume
	@echo "Submission: $(OUT_VCF)/submission/submission.csv"

track1:
	$(NXF) --outdir $(OUT_FULL) -resume
	@echo "Submission: $(OUT_FULL)/submission/submission.csv"

# ---- Track 2 ----------------------------------------------------------------
# The case is resolved first (variants -> consequences, gene file generated if
# missing); the steps then run in a second make, which reads the result.
# Steps whose input is absent are skipped: no stop -> no NMD check, no
# missense -> no ESM/structure, no pathway panel -> no drug screen.
track2: case
	$(MAKE) --no-print-directory track2-steps
	@echo "Track 2 outputs in $(T2)/"

# case.mk deliberately has no rule of its own: make would otherwise try to
# rebuild it (and the whole case) before every command, even `make help`.
case: $(RCASE)
	mkdir -p $(BUILD)
	$(PY) bin/caseconfig.py --make $(RCASE) > $(BUILD)/case.mk
$(RCASE): $(CASE) $(GTF) $(wildcard $(SUBMISSION)) $(wildcard config/genes/*.yaml)
	$(PY) bin/resolve_case.py --case $(CASE) --gtf $(GTF) --fasta $(FASTA) \
	    --submission $(SUBMISSION) --out $@

track2-steps: $(if $(PATIENT_PTC_POS),nmd) $(if $(PATIENT_MISSENSE),esm structure) \
	    $(if $(PANEL),opentargets) $(if $(PATIENT_PTC_POS)$(PATIENT_MISSENSE),spliceai) \
	    atlas clinvar

nmd: $(T2)/nmd_$(GENE)_$(PATIENT_PTC_POS).txt
$(T2)/nmd_$(GENE)_$(PATIENT_PTC_POS).txt: $(GTF)
	mkdir -p $(T2)
	$(PY) bin/nmd_check.py --gtf $(GTF) --gene $(GENE) --pos $(PATIENT_PTC_POS) \
	    $(if $(DOMAINS),--domains $(DOMAINS)) | tee $@

esm: $(T2)/esm_scores.tsv
$(T2)/esm_scores.tsv: resources/$(UNIPROT).fasta
	mkdir -p $(T2)
	$(PY) bin/esm_variant_score.py --config $(RCASE) --fasta $< --out $@

structure: $(T2)/structure_context.tsv
$(T2)/structure_context.tsv: resources/AF-$(UNIPROT)-F1-model.pdb
	mkdir -p $(T2)
	$(PY) bin/structure_context.py --config $(RCASE) --pdb $< --out $@

opentargets: $(T2)/opentargets/ot_target_drugs.tsv
$(T2)/opentargets/ot_target_drugs.tsv:
	$(PY) bin/opentargets_candidates.py --config $(RCASE) --outdir $(T2)/opentargets

atlas: $(T2)/atlas/atlas_figure.pdf
$(T2)/atlas/atlas_figure.pdf: $(GTF) resources/$(UNIPROT).fasta resources/AF-$(UNIPROT)-F1-model.pdb
	$(PY) bin/variant_atlas.py --config $(RCASE) --fasta resources/$(UNIPROT).fasta \
	    --gtf $(GTF) --pdb resources/AF-$(UNIPROT)-F1-model.pdb --outdir $(T2)/atlas

# Atlas for any gene, without patient data: make gene-atlas G=CEP57.
# The gene file is generated from UniProt if it does not exist yet.
G   ?=
GLC  = $(shell echo $(G) | tr A-Z a-z)
GCFG = config/genes/$(GLC).yaml
GUNI = $(if $(wildcard $(GCFG)),$(shell $(CFG) $(GCFG) uniprot))
GA   = $(FINDINGS)/genes/$(shell echo $(G) | tr a-z A-Z)

gene-config:
	@test -n "$(G)" || { echo "usage: make gene-config G=SYMBOL"; exit 1; }
	$(MAKE) --no-print-directory $(GCFG)
config/genes/%.yaml:
	$(PY) bin/make_gene_config.py --gene $* --out $@

gene-atlas: gene-config
	$(MAKE) --no-print-directory $(GA)/atlas_figure.pdf G=$(G)
$(FINDINGS)/genes/%/atlas_figure.pdf: $(GTF) resources/$(GUNI).fasta resources/AF-$(GUNI)-F1-model.pdb
	$(PY) bin/variant_atlas.py --config $(GCFG) --fasta resources/$(GUNI).fasta \
	    --gtf $(GTF) --pdb resources/AF-$(GUNI)-F1-model.pdb --outdir $(@D)

# Genes causing other forms of the same disease: the OMIM series and the
# Orphanet group of the Track 2 gene (bin/related_genes.py). One gene's atlas
# failing (e.g. not in the annotation) is reported and does not stop the rest.
related-genes: $(if $(filter auto,$(RELATED)),case)
	$(MAKE) --no-print-directory related-genes-run
related-genes-run: $(if $(filter auto,$(RELATED)),$(BUILD)/related_genes.txt)
	@failed=""; \
	for g in $(if $(filter auto,$(RELATED)),$$(cat $(BUILD)/related_genes.txt),$(RELATED)); do \
	    $(MAKE) --no-print-directory gene-atlas G=$$g || failed="$$failed $$g"; done; \
	if [ -n "$$failed" ]; then echo "WARNING: no atlas for:$$failed (see messages above)"; fi
mva-genes: related-genes

$(BUILD)/related_genes.txt: $(RCASE) $(HPO)/genes_to_disease.txt $(HPO)/phenotype.hpoa
	mkdir -p $(@D)
	$(PY) bin/related_genes.py --gene $(GENE) --g2d $(HPO)/genes_to_disease.txt \
	    --hpoa $(HPO)/phenotype.hpoa --out $(T2)/related_genes.tsv > $@.tmp
	mv $@.tmp $@

# Public protein inputs for any UniProt accession
resources/%.fasta resources/AF-%-F1-model.pdb:
	UNIPROT=$* bash scripts/prepare_inputs.sh public

# Independent check of the atlas: where do ClinVar's BUB1B missense variants
# (pathogenic, benign, uncertain) fall? Needs the ClinVar file from `make refs`.
CLINVAR ?= $(REFS)/clinvar/clinvar.chr.vcf.gz
clinvar: $(T2)/clinvar_validation.tsv
$(T2)/clinvar_validation.tsv: $(T2)/atlas/atlas_figure.pdf
	$(PY) bin/clinvar_validate.py --clinvar $(CLINVAR) --fasta $(FASTA) --gtf $(GTF) \
	    --config $(RCASE) --protein resources/$(UNIPROT).fasta \
	    --atlas $(T2)/atlas/atlas_missense.tsv --out $@

# SpliceAI on the two causal alleles only: rules out a splicing mechanism.
# Own environment: SpliceAI needs TensorFlow 2.15, which needs Python 3.11.
SPLICEAI_PY := .venv-spliceai/bin/python
spliceai: $(T2)/spliceai_two_alleles.tsv
$(T2)/spliceai_two_alleles.tsv: $(SPLICEAI_PY)
	mkdir -p $(T2) $(BUILD)
	len=$$(awk -v c=$(PATIENT_CHROM) '$$1==c{print $$2}' $(FASTA).fai); \
	printf '##fileformat=VCFv4.2\n##contig=<ID=%s,length=%s>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' \
	    $(PATIENT_CHROM) $$len > $(BUILD)/two.vcf; \
	for a in "$(PATIENT_PTC_ALLELE)" "$(PATIENT_MIS_ALLELE)"; do \
	    set -- $$a; printf '%s\t%s\t.\t%s\t%s\t.\t.\t.\n' $(PATIENT_CHROM) $$1 $$2 $$3 >> $(BUILD)/two.vcf; \
	done
	.venv-spliceai/bin/spliceai -I $(BUILD)/two.vcf -O $(BUILD)/two.spliceai.vcf \
	    -R $(FASTA) -A grch38 -D 4999
	{ printf 'chrom\tpos\tref\talt\tgene\tDS_AG\tDS_AL\tDS_DG\tDS_DL\n'; \
	  grep -v '^#' $(BUILD)/two.spliceai.vcf | awk -F'\t' '{ \
	    split($$8, kv, "SpliceAI="); n = split(kv[2], ann, ","); \
	    for (i = 1; i <= n; i++) { split(ann[i], f, "|"); \
	      print $$1"\t"$$2"\t"$$4"\t"$$5"\t"f[2]"\t"f[3]"\t"f[4]"\t"f[5]"\t"f[6] } }'; \
	} | tee $@
	rm -f $(BUILD)/two.vcf $(BUILD)/two.spliceai.vcf

$(SPLICEAI_PY):
	uv venv .venv-spliceai --python 3.11
	uv pip install --python $@ spliceai "tensorflow-cpu==2.15.1" "setuptools<81"

# Aneuploidy under three GC-correction models (needs a finished `make track1`)
aneuploidy-sensitivity:
	mkdir -p $(OUT_FULL)/aneuploidy_sensitivity
	for m in local 2d poly; do \
	    $(PY) bin/aneuploidy_plot.py \
	        --bins $$(ls $(OUT_FULL)/qc/mosdepth/*.regions.bed.gz | head -1) \
	        --fasta $(FASTA) --vcf $$(ls $(OUT_FULL)/variants/*.dv.vcf.gz | head -1) \
	        --gc-model $$m \
	        --out-figure $(OUT_FULL)/aneuploidy_sensitivity/karyotype_$$m.pdf \
	        --out-table $(OUT_FULL)/aneuploidy_sensitivity/mosaic_$$m.tsv; \
	done

# ---- sharing and compliance -------------------------------------------------
# Copies ONLY findings (named variants, summaries, figures) out of the results
# folders. Nothing at genotype scale is copied.
findings:
	mkdir -p $(FINDINGS)/track1 $(FINDINGS)/aneuploidy
	-cp $(OUT_VCF)/submission/submission.csv $(FINDINGS)/track1/submission_vcf.csv
	-cp $(OUT_FULL)/submission/submission.csv $(FINDINGS)/track1/submission_full.csv
	-cp $(OUT_FULL)/aneuploidy/*.pdf $(OUT_FULL)/aneuploidy/*.png \
	    $(OUT_FULL)/aneuploidy/*.tsv* $(FINDINGS)/aneuploidy/
	-cp $(OUT_FULL)/aneuploidy_sensitivity/mosaic_*.tsv* $(FINDINGS)/aneuploidy/
	-cp $(OUT_FULL)/pipeline_info/*.html $(OUT_FULL)/pipeline_info/trace.tsv \
	    $(FINDINGS)/track1/
	@echo "Review $(FINDINGS)/ before committing."

versions:
	@mkdir -p $(FINDINGS)
	@{ echo "# Tool and container versions, $$(date -u +%F)"; \
	   nextflow -version 2>&1 | grep -m1 version; \
	   for t in samtools bcftools fastp mosdepth delly apptainer java; do \
	       printf '%s: ' $$t; ($$t --version 2>&1 || $$t 2>&1) | grep -m1 -iE 'version|[0-9]+\.[0-9]+' ; done; \
	   printf 'bwa: '; bwa 2>&1 | grep -m1 '^Version'; \
	   printf 'ExpansionHunter: '; ExpansionHunter --version 2>&1 | grep -m1 -oE 'v[0-9]+\.[0-9]+\.[0-9]+'; \
	   echo "## containers (sha256)"; sha256sum $(CONTAINERS)/*.sif; \
	   echo "## python"; $(PY) -m pip freeze 2>/dev/null || uv pip freeze; \
	} | tee $(FINDINGS)/versions.txt

deletion-check:
	@echo "Patient-derived files still on this machine (see docs/DELETION.md):"
	@for p in data derived testdata work results results_* ~/.cache/huggingface; do \
	    [ -e "$$p" ] && du -sh "$$p" 2>/dev/null; done; true

clean:
	rm -rf work .nextflow* results_smoke $(BUILD)
