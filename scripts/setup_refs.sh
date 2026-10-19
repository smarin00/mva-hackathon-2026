#!/usr/bin/env bash
# Download and prepare every reference resource, at the releases used for the
# submission. Reconstructed from the analysis server's setup history.
#
# Space: ~100 GB in $REFS (genome + BWA index ~9 GB, VEP cache ~25 GB,
# Exomiser data ~60 GB). Time: 1-3 h, dominated by downloads.
# Each step is skipped if its output already exists, so the script can be
# re-run after an interruption.
#
#   REFS=/opt/refs CONTAINERS=/opt/containers bash scripts/setup_refs.sh
set -euo pipefail

REFS=${REFS:-/opt/refs}
CONTAINERS=${CONTAINERS:-/opt/containers}
REPO=$(cd "$(dirname "$0")/.." && pwd)
CHR_MAP="$REPO/resources/chr_map.txt"

ENSEMBL=116
EXOMISER_VER=15.1.0
EXOMISER_DATA=2602
EXOMISER_DATA_URL=https://g-879a9f.f5dc97.75bc.dn.glob.us/data   # Exomiser discussion #640
NCBI=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids
GENOME=GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set

for d in "$REFS" "$CONTAINERS"; do
    [ -d "$d" ] || { sudo mkdir -p "$d" && sudo chown "$USER" "$d"; }
done
step() { echo; echo "== $*"; }

# --- 1. Genome: GRCh38 no-alt analysis set + hs38d1 decoys, prebuilt BWA index
step "GRCh38 analysis set"
cd "$REFS"
if [ ! -s "$GENOME.fna" ]; then
    wget -q "$NCBI/$GENOME.fna.gz" && gunzip "$GENOME.fna.gz"
fi
[ -s "$GENOME.fna.fai" ] || wget -q "$NCBI/$GENOME.fna.fai"
if [ ! -s "$GENOME.fna.bwt" ]; then
    wget -q "$NCBI/$GENOME.fna.bwa_index.tar.gz"
    tar -xzf "$GENOME.fna.bwa_index.tar.gz" && rm "$GENOME.fna.bwa_index.tar.gz"
fi
[ -s "$GENOME.dict" ] || samtools dict "$GENOME.fna" -o "$GENOME.dict"
# Stable short names used by nextflow.config (params.fasta = $REFS/GRCh38.fa)
ln -sf "$GENOME.dict" GRCh38.dict
for ext in "" .fai .amb .ann .bwt .pac .sa; do
    ln -sf "$GENOME.fna$ext" "GRCh38.fa$ext"
done

# --- 2. Gene annotation (panel building, NMD geometry, variant atlas)
step "Ensembl $ENSEMBL GTF"
[ -s "Homo_sapiens.GRCh38.$ENSEMBL.gtf.gz" ] || wget -q \
    "https://ftp.ensembl.org/pub/release-$ENSEMBL/gtf/homo_sapiens/Homo_sapiens.GRCh38.$ENSEMBL.gtf.gz"

# --- 3. VEP: container, offline cache, plugins
step "VEP $ENSEMBL (container, cache, plugins)"
[ -s "$CONTAINERS/vep_$ENSEMBL.sif" ] || apptainer pull \
    "$CONTAINERS/vep_$ENSEMBL.sif" "docker://ensemblorg/ensembl-vep:release_$ENSEMBL.2"
mkdir -p "$REFS/vep" && cd "$REFS/vep"
if [ ! -d "homo_sapiens/${ENSEMBL}_GRCh38" ]; then
    wget -q "https://ftp.ensembl.org/pub/release-$ENSEMBL/variation/indexed_vep_cache/homo_sapiens_vep_${ENSEMBL}_GRCh38.tar.gz"
    tar -xzf "homo_sapiens_vep_${ENSEMBL}_GRCh38.tar.gz" && rm "homo_sapiens_vep_${ENSEMBL}_GRCh38.tar.gz"
fi
[ -d Plugins ] || git clone -q --depth 1 -b "release/$ENSEMBL" \
    https://github.com/Ensembl/VEP_plugins.git Plugins

# --- 4. AlphaMissense (VEP plugin input)
step "AlphaMissense"
mkdir -p "$REFS/alphamissense" && cd "$REFS/alphamissense"
[ -s AlphaMissense_hg38.tsv.gz ] || wget -q \
    https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
[ -s AlphaMissense_hg38.tsv.gz.tbi ] || tabix -s 1 -b 2 -e 2 -f -S 1 AlphaMissense_hg38.tsv.gz

# --- 5. ClinVar, renamed to chr contigs. NCBI serves only the current weekly
#        release; the fileDate actually used is recorded in docs/reference_data.md.
step "ClinVar (current release)"
mkdir -p "$REFS/clinvar" && cd "$REFS/clinvar"
if [ ! -s clinvar.chr.vcf.gz ]; then
    wget -q https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
    tabix -p vcf clinvar.vcf.gz
    bcftools view -r "$(cut -f1 "$CHR_MAP" | paste -sd,)" clinvar.vcf.gz \
        | bcftools annotate --rename-chrs "$CHR_MAP" -O z -o clinvar.chr.vcf.gz
    tabix -p vcf clinvar.chr.vcf.gz
fi
echo "ClinVar $(bcftools view -h clinvar.chr.vcf.gz | grep -m1 '^##fileDate')"

# --- 6. Exomiser CLI + data release
step "Exomiser $EXOMISER_VER, data $EXOMISER_DATA"
mkdir -p "$REFS/exomiser" "$REFS/exomiser-data"
cd "$REFS/exomiser"
if [ ! -d "exomiser-cli-$EXOMISER_VER" ]; then
    wget -q "https://github.com/exomiser/Exomiser/releases/download/$EXOMISER_VER/exomiser-cli-$EXOMISER_VER-distribution.zip"
    unzip -q "exomiser-cli-$EXOMISER_VER-distribution.zip"
fi
cd "$REFS/exomiser-data"
for part in hg38 phenotype; do
    if [ ! -d "${EXOMISER_DATA}_$part" ]; then
        wget -q "$EXOMISER_DATA_URL/${EXOMISER_DATA}_$part.zip"
        unzip -q "${EXOMISER_DATA}_$part.zip" && rm "${EXOMISER_DATA}_$part.zip"
    fi
done
cat > "$REFS/exomiser/exomiser-cli-$EXOMISER_VER/application.properties" <<PROPS
exomiser.data-directory=$REFS/exomiser-data
remm.version=0.4
cadd.version=1.7
exomiser.hg38.data-version=$EXOMISER_DATA
exomiser.phenotype.data-version=$EXOMISER_DATA
logging.level.com.zaxxer.hikari=ERROR
PROPS

# --- 7. ExpansionHunter catalog ships with the binary (scripts/install_tools.sh)

# --- 8. Containers
step "Containers (DeepVariant 1.10.0, GATK 4.7.0.0)"
[ -s "$CONTAINERS/deepvariant_1.10.0.sif" ] || apptainer pull \
    "$CONTAINERS/deepvariant_1.10.0.sif" docker://google/deepvariant:1.10.0
[ -s "$CONTAINERS/gatk_4.7.0.0.sif" ] || apptainer pull \
    "$CONTAINERS/gatk_4.7.0.0.sif" docker://broadinstitute/gatk:4.7.0.0

echo
echo "References ready in $REFS, containers in $CONTAINERS."
echo "Record container digests for provenance:  make versions"
