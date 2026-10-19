# MVA Hackathon 2026 -- single entry point.  Run `make` for the menu.
#
# Everything case-specific is in the block below; override on the command
# line, e.g.  make track2 PATIENT_MISSENSE=R550Q

# ---- case -----------------------------------------------------------------
GENE              ?= BUB1B
UNIPROT           ?= O60566
# GRCh38 position (VCF POS) of the stop-gain, and its codon
PATIENT_PTC_POS   ?= 40209701
PATIENT_PTC_CODON ?= 737
PATIENT_MISSENSE  ?= N1002K
# Genomic alleles (GRCh38) of the two causal variants, for splice scoring
PATIENT_CHROM     ?= chr15
PATIENT_PTC_ALLELE ?= 40209701 T G
PATIENT_MIS_ALLELE ?= 40220612 T G

# ---- locations ------------------------------------------------------------
REFS       ?= /opt/refs
CONTAINERS ?= /opt/containers
GTF        ?= $(REFS)/Homo_sapiens.GRCh38.116.gtf.gz
FASTA      ?= $(REFS)/GRCh38.fa
OUT_VCF    ?= results_vcf
OUT_FULL   ?= results_full
FINDINGS   ?= findings
PY         ?= .venv/bin/python
NXF        ?= nextflow run main.nf -params-file params.yaml

T2 := $(FINDINGS)/track2


.PHONY: help install refs gtf setup data inputs lint smoke track1-vcf track1 \
	    track2 nmd esm structure opentargets atlas spliceai aneuploidy-sensitivity \
	    findings versions deletion-check clean

help:
	@echo "MVA Hackathon 2026 pipeline. Typical order:"
	@echo ""
	@echo "  make install      system tools (sudo): Nextflow, Apptainer, bwa, samtools, ..."
	@echo "  make refs         reference genome, VEP, Exomiser, ClinVar, containers (~100 GB)"
	@echo "  make gtf          only the Ensembl annotation (enough for Track 2, 60 MB)"
	@echo "  make setup        Python environment (.venv)"
	@echo "  make data         download the patient dataset (needs approved access)"
	@echo "  make inputs       prepare VCF, gene panel and public protein files"
	@echo ""
	@echo "  make track1-vcf   Track 1 from the provided VCF       (~2 h)  -> $(OUT_VCF)/"
	@echo "  make track1       Track 1 from raw reads (FASTQ)      (~15 h) -> $(OUT_FULL)/"
	@echo "  make track2       Track 2 mechanism + drug analyses   (~1.5 h) -> $(T2)/"
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

# Track 2 needs only the gene annotation, not the full reference set
gtf: $(GTF)
$(GTF):
	mkdir -p $(dir $(GTF))
	wget -q -O $@ https://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz

setup: .venv/.done
.venv/.done: requirements.lock
	uv venv --python 3.12
	uv pip sync requirements.lock
	uv pip install torch --index-url https://download.pytorch.org/whl/cpu
	touch $@

data:
	bash scripts/prepare_inputs.sh data

inputs:
	REFS=$(REFS) UNIPROT=$(UNIPROT) bash scripts/prepare_inputs.sh vcf panel public
	@test -f params.yaml || { cp params.example.yaml params.yaml; \
	    echo "Created params.yaml from the template: add the HPO terms before running."; }

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
track2: nmd esm structure opentargets atlas spliceai
	@echo "Track 2 outputs in $(T2)/"

nmd: $(T2)/nmd_$(GENE)_$(PATIENT_PTC_POS).txt
$(T2)/nmd_$(GENE)_$(PATIENT_PTC_POS).txt: $(GTF)
	mkdir -p $(T2)
	$(PY) bin/nmd_check.py --gtf $(GTF) --gene $(GENE) --pos $(PATIENT_PTC_POS) \
	    --domains docs/bub1b_domains.tsv | tee $@

esm: $(T2)/esm_scores.tsv
$(T2)/esm_scores.tsv: resources/$(UNIPROT).fasta
	mkdir -p $(T2)
	$(PY) bin/esm_variant_score.py --fasta $< --patient $(PATIENT_MISSENSE) --out $@

structure: $(T2)/structure_context.tsv
$(T2)/structure_context.tsv: resources/AF-$(UNIPROT)-F1-model.pdb
	mkdir -p $(T2)
	$(PY) bin/structure_context.py --pdb $< --patient $(PATIENT_MISSENSE) --out $@

opentargets: $(T2)/opentargets/ot_target_drugs.tsv
$(T2)/opentargets/ot_target_drugs.tsv:
	$(PY) bin/opentargets_candidates.py --outdir $(T2)/opentargets

atlas: $(T2)/atlas/atlas_figure.pdf
$(T2)/atlas/atlas_figure.pdf: $(GTF) resources/$(UNIPROT).fasta resources/AF-$(UNIPROT)-F1-model.pdb
	$(PY) bin/bub1b_atlas.py --fasta resources/$(UNIPROT).fasta --gtf $(GTF) \
	    --gene $(GENE) --domains docs/bub1b_domains.tsv \
	    --pdb resources/AF-$(UNIPROT)-F1-model.pdb \
	    --patient-missense $(PATIENT_MISSENSE) --patient-ptc $(PATIENT_PTC_CODON) \
	    --outdir $(T2)/atlas

resources/$(UNIPROT).fasta resources/AF-$(UNIPROT)-F1-model.pdb:
	UNIPROT=$(UNIPROT) bash scripts/prepare_inputs.sh public

# SpliceAI on the two causal alleles only: rules out a splicing mechanism.
# Own environment: SpliceAI needs TensorFlow 2.15, which needs Python 3.11.
SPLICEAI_PY := .venv-spliceai/bin/python
spliceai: $(T2)/spliceai_two_alleles.tsv
$(T2)/spliceai_two_alleles.tsv: $(SPLICEAI_PY)
	mkdir -p $(T2)
	len=$$(awk -v c=$(PATIENT_CHROM) '$$1==c{print $$2}' $(FASTA).fai); \
	printf '##fileformat=VCFv4.2\n##contig=<ID=%s,length=%s>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' \
	    $(PATIENT_CHROM) $$len > $(T2)/.two.vcf; \
	for a in "$(PATIENT_PTC_ALLELE)" "$(PATIENT_MIS_ALLELE)"; do \
	    set -- $$a; printf '%s\t%s\t.\t%s\t%s\t.\t.\t.\n' $(PATIENT_CHROM) $$1 $$2 $$3 >> $(T2)/.two.vcf; \
	done
	.venv-spliceai/bin/spliceai -I $(T2)/.two.vcf -O $(T2)/.two.spliceai.vcf \
	    -R $(FASTA) -A grch38 -D 4999
	{ printf 'chrom\tpos\tref\talt\tgene\tDS_AG\tDS_AL\tDS_DG\tDS_DL\n'; \
	  grep -v '^#' $(T2)/.two.spliceai.vcf | awk -F'\t' '{ \
	    split($$8, kv, "SpliceAI="); n = split(kv[2], ann, ","); \
	    for (i = 1; i <= n; i++) { split(ann[i], f, "|"); \
	      print $$1"\t"$$2"\t"$$4"\t"$$5"\t"f[2]"\t"f[3]"\t"f[4]"\t"f[5]"\t"f[6] } }'; \
	} | tee $@
	rm -f $(T2)/.two.vcf $(T2)/.two.spliceai.vcf

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
	rm -rf work .nextflow* results_smoke
