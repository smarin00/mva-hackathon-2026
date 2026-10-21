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
    # try URL OUT: download to OUT.part, rename on success; returns curl's status
    try() {
        curl -sSfL --retry 3 -A "mva-hackathon-pipeline (curl)" "$1" -o "$2.part" \
            && mv "$2.part" "$2" || { rm -f "$2.part"; return 1; }
    }
    fetch() { try "$1" "$2" || { echo "ERROR: download failed: $1" >&2; exit 1; }; }
    PY=python3; [ -x .venv/bin/python ] && PY=.venv/bin/python

    fetch "https://rest.uniprot.org/uniprotkb/$UNIPROT.fasta" "resources/$UNIPROT.fasta"
    fetch "https://alphafold.ebi.ac.uk/api/prediction/$UNIPROT" "resources/AF-$UNIPROT.json"
    urls=$($PY -c 'import json,sys
d = json.load(open(sys.argv[1]))
if not d: sys.exit("ERROR: AlphaFold has no prediction for " + sys.argv[2])
print(d[0]["pdbUrl"], d[0]["cifUrl"])' "resources/AF-$UNIPROT.json" "$UNIPROT")         || { rm -f "resources/AF-$UNIPROT.json"; exit 1; }
    url=${urls% *}; cif=${urls#* }
    rm -f "resources/AF-$UNIPROT.json"
    pdb="resources/AF-$UNIPROT-F1-model.pdb"
    if ! try "$url" "$pdb"; then
        # Some entries are served as mmCIF only: fetch it and convert.
        echo "   PDB file not available; using mmCIF and converting" >&2
        url=$cif
        fetch "$url" "resources/AF-$UNIPROT-F1-model.cif"
        $PY - "resources/AF-$UNIPROT-F1-model.cif" "$pdb" <<'EOF'
import sys
from Bio.PDB import MMCIFParser, PDBIO
io = PDBIO()
io.set_structure(MMCIFParser(QUIET=True).get_structure("m", sys.argv[1]))
io.save(sys.argv[2])
EOF
        rm -f "resources/AF-$UNIPROT-F1-model.cif"
    fi
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
