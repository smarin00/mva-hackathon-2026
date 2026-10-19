#!/usr/bin/env bash
# Install the command-line tools the pipeline needs, at the versions used for
# the hackathon submission. Reconstructed from the setup shell history of the
# analysis server (Ubuntu 26.04, Hetzner CPX62).
#
# Needs sudo for apt and Apptainer. Everything else goes to $TOOLS (default
# /opt/tools), which this script puts on PATH via ~/.bashrc.
#
#   bash scripts/install_tools.sh
set -euo pipefail

TOOLS=${TOOLS:-/opt/tools}
NXF_VER=${NXF_VER:-26.04.6}
MOSDEPTH_VER=0.3.14
DELLY_VER=2.7.0
EH_VER=5.0.0

echo "== apt packages (samtools, bcftools, tabix, bwa, fastp, Java 21, build deps)"
sudo apt-get update
sudo apt-get install -y build-essential git tmux curl wget unzip pkg-config \
    zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev \
    libncurses-dev libdeflate-dev openjdk-21-jre-headless \
    samtools bcftools tabix bwa fastp

echo "== Apptainer (containers for DeepVariant, GATK, VEP)"
if ! command -v apptainer >/dev/null; then
    sudo add-apt-repository -y ppa:apptainer/ppa
    sudo apt-get update && sudo apt-get install -y apptainer
fi

sudo mkdir -p "$TOOLS" && sudo chown "$USER" "$TOOLS"
cd "$TOOLS"

echo "== Nextflow $NXF_VER"
if [ ! -x nextflow ]; then
    curl -s https://get.nextflow.io | NXF_VER=$NXF_VER bash
fi

echo "== mosdepth $MOSDEPTH_VER"
[ -x mosdepth ] || { wget -q -O mosdepth \
    "https://github.com/brentp/mosdepth/releases/download/v${MOSDEPTH_VER}/mosdepth" && chmod +x mosdepth; }

echo "== delly $DELLY_VER"
[ -x delly ] || { wget -q -O delly \
    "https://github.com/dellytools/delly/releases/download/v${DELLY_VER}/delly-v${DELLY_VER}-linux-amd64" && chmod +x delly; }

echo "== ExpansionHunter $EH_VER"
if [ ! -d "ExpansionHunter-v${EH_VER}-linux_x86_64" ]; then
    wget -q "https://github.com/Illumina/ExpansionHunter/releases/download/v${EH_VER}/ExpansionHunter-v${EH_VER}-linux_x86_64.tar.gz"
    tar -xzf "ExpansionHunter-v${EH_VER}-linux_x86_64.tar.gz"
fi
ln -sf "$TOOLS/ExpansionHunter-v${EH_VER}-linux_x86_64/bin/ExpansionHunter" ExpansionHunter

echo "== uv (Python environment manager)"
command -v uv >/dev/null || curl -LsSf https://astral.sh/uv/install.sh | sh

line="export PATH=\"$TOOLS:\$PATH\""
grep -qxF "$line" ~/.bashrc || echo "$line" >> ~/.bashrc

echo
echo "Done. Open a new shell (or: source ~/.bashrc), then check versions with:"
echo "  make versions"
