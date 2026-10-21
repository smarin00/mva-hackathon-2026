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
# Case-specific values come from the case file via make (see params.example.yaml)
DATASET=${DATASET:-}
PANEL_GENES=${PANEL_GENES:-docs/panel_genes.txt}
PANEL_BED=${PANEL_BED:-resources/mva_panel.bed}
# the dataset's single VCF and its lane-1 reads, whatever the sample is called
VENDOR_VCF=${VENDOR_VCF:-$(ls data/*.vcf.gz 2>/dev/null | head -1 || true)}
UNIPROT=${UNIPROT:-}

for part in "$@"; do
case "$part" in
data)
    [ -n "$DATASET" ] || { echo "ERROR: set dataset: in the case file (Hugging Face dataset id)"; exit 1; }
    echo "== downloading $DATASET into data/ (FASTQ ~100 GB)"
    hf download "$DATASET" --repo-type dataset --local-dir ./data \
        --include "*.fastq.gz" "*.vcf.gz" "*.vcf.gz.tbi" "*.docx" "README.md"
    ;;
vcf)
    echo "== provided VCF -> chr-named contigs (primary chromosomes only)"
    [ -n "$VENDOR_VCF" ] || { echo "ERROR: no data/*.vcf.gz (run: make data)"; exit 1; }
    mkdir -p derived
    bcftools view -r "$(cut -f1 resources/chr_map.txt | paste -sd,)" "$VENDOR_VCF" \
        | bcftools annotate --rename-chrs resources/chr_map.txt \
            -O z -o derived/vendor.chr.vcf.gz
    tabix -f -p vcf derived/vendor.chr.vcf.gz
    echo "   $(bcftools index -n derived/vendor.chr.vcf.gz) records"
    ;;
panel)
    echo "== candidate panel from Ensembl GTF"
    echo "   genes: $PANEL_GENES -> $PANEL_BED"
    python3 bin/build_panel.py --genes "$PANEL_GENES" \
        --gtf "$REFS/Homo_sapiens.GRCh38.116.gtf.gz" --pad 10000 \
        --out "$PANEL_BED"
    ;;
public)
    [ -n "$UNIPROT" ] || { echo "ERROR: set UNIPROT=<accession> (make does this per gene)"; exit 1; }
    echo "== UniProt $UNIPROT sequence and AlphaFold model (public)"
    # Downloads use Python's urllib, the client that fetched the BUB1B model
    # in structure_context.py. EBI's AlphaFold file server refused curl (403)
    # for the same URLs.
    PY=python3; [ -x .venv/bin/python ] && PY=.venv/bin/python
    $PY - "$UNIPROT" <<'EOF'
import json, os, sys, urllib.request
acc = sys.argv[1]

def get(url, dest=None):
    """Download url (to dest via dest.part, or return bytes); exit with the URL on failure."""
    try:
        with urllib.request.urlopen(url, timeout=120) as r:
            data = r.read()
    except Exception as e:
        sys.exit(f"ERROR: download failed: {url} ({e})")
    if dest is None:
        return data
    with open(dest + ".part", "wb") as fh:
        fh.write(data)
    os.replace(dest + ".part", dest)

get(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", f"resources/{acc}.fasta")
entries = json.loads(get(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}"))
if not entries:
    sys.exit(f"ERROR: AlphaFold has no prediction for {acc}")
e = entries[0]
pdb = f"resources/AF-{acc}-F1-model.pdb"
try:
    get(e["pdbUrl"], pdb)
    print("   " + e["pdbUrl"])
except SystemExit as err:
    # Some entries are served as mmCIF only: fetch it and convert (pLDDT is
    # kept in the B-factor column).
    print(f"   {err}; trying mmCIF", file=sys.stderr)
    cif = pdb[:-4] + ".cif"
    get(e["cifUrl"], cif)
    from Bio.PDB import MMCIFParser, PDBIO
    io = PDBIO()
    io.set_structure(MMCIFParser(QUIET=True).get_structure("m", cif))
    io.save(pdb)
    os.remove(cif)
    print("   " + e["cifUrl"] + " (converted to PDB)")
EOF
    ;;
smoke)
    echo "== testdata/: first 2 M read pairs of lane 1 (PATIENT DATA: never commit)"
    mkdir -p testdata
    for R in 1 2; do
        fq=$(ls data/*_L001_R${R}_001.fastq.gz 2>/dev/null | head -1 || true)
        [ -n "$fq" ] || { echo "ERROR: no data/*_L001_R${R}_001.fastq.gz (run: make data)"; exit 1; }
        zcat "$fq" | head -n 8000000 \
            | gzip > "testdata/TEST_S16_L001_R${R}_001.fastq.gz" || true
    done
    ;;
*)  echo "unknown part: $part (data | vcf | panel | public | smoke)"; exit 1 ;;
esac
done
