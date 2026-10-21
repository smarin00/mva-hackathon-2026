#!/usr/bin/env bash
# Prepare pipeline inputs. Three independent parts, selected by argument:
#
#   data     download the patient dataset from Hugging Face (needs approved
#            access to SageBio/mva-hackathon-2026-data and `hf auth login`)
#   vcf      rename the provided VCF's contigs to chr-style -> derived/
#   panel    build the candidate-gene panel BED from the Ensembl GTF
#   public   fetch the public protein inputs for Track 2 (UniProt, AlphaFold)
#   smoke    cut a small FASTQ subset into testdata/ for a quick pipeline test
#
#   bash scripts/prepare_inputs.sh data vcf panel public
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
cd "$REPO"
REFS=${REFS:-/opt/refs}
DATASET=SageBio/mva-hackathon-2026-data
VENDOR_VCF=data/WGS_EX2312012_HGWCNDSX7.vcf.gz
UNIPROT=${UNIPROT:-O60566}

for part in "$@"; do
case "$part" in
data)
    echo "== downloading $DATASET into data/ (FASTQ ~100 GB)"
    hf download "$DATASET" --repo-type dataset --local-dir ./data \
        --include "*.fastq.gz" "*.vcf.gz" "*.vcf.gz.tbi" "*.docx" "README.md"
    ;;
vcf)
    echo "== provided VCF -> chr-named contigs (primary chromosomes only)"
    mkdir -p derived
    bcftools view -r "$(cut -f1 resources/chr_map.txt | paste -sd,)" "$VENDOR_VCF" \
        | bcftools annotate --rename-chrs resources/chr_map.txt \
            -O z -o derived/vendor.chr.vcf.gz
    tabix -f -p vcf derived/vendor.chr.vcf.gz
    echo "   $(bcftools index -n derived/vendor.chr.vcf.gz) records"
    ;;
panel)
    echo "== candidate panel from Ensembl GTF"
    python3 bin/build_panel.py --genes docs/panel_genes.txt \
        --gtf "$REFS/Homo_sapiens.GRCh38.116.gtf.gz" --pad 10000 \
        --out resources/mva_panel.bed
    ;;
public)
    echo "== UniProt $UNIPROT sequence and AlphaFold model (public)"
    # fetch URL OUT: fail loudly with the URL; never leave a partial file
    fetch() {
        curl -sSfL --retry 3 "$1" -o "$2.part" \
            || { rm -f "$2.part"; echo "ERROR: download failed: $1" >&2; exit 1; }
        mv "$2.part" "$2"
    }
    fetch "https://rest.uniprot.org/uniprotkb/$UNIPROT.fasta" "resources/$UNIPROT.fasta"
    fetch "https://alphafold.ebi.ac.uk/api/prediction/$UNIPROT" "resources/AF-$UNIPROT.json"
    url=$(python3 -c 'import json,sys
d = json.load(open(sys.argv[1]))
if not d: sys.exit("ERROR: AlphaFold has no prediction for " + sys.argv[2])
print(d[0]["pdbUrl"])' "resources/AF-$UNIPROT.json" "$UNIPROT")
    rm -f "resources/AF-$UNIPROT.json"
    fetch "$url" "resources/AF-$UNIPROT-F1-model.pdb"
    echo "   $url"
    ;;
smoke)
    echo "== testdata/: first 2 M read pairs of lane 1 (PATIENT DATA: never commit)"
    mkdir -p testdata
    for R in 1 2; do
        zcat "data/WGS_EX2312012_HGWCNDSX7_S16_L001_R${R}_001.fastq.gz" | head -n 8000000 \
            | gzip > "testdata/TEST_S16_L001_R${R}_001.fastq.gz" || true
    done
    ;;
*)  echo "unknown part: $part (data | vcf | panel | public | smoke)"; exit 1 ;;
esac
done
