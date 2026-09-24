#!/bin/bash
### General options
#BSUB -q hpc
#BSUB -J m2_fetch_lotus
### -- download ~170MB, build a 554 Mbp bowtie2 index: light --
#BSUB -n 8
#BSUB -R "span[hosts=1] rusage[mem=4GB]"
#BSUB -M 4500MB
#BSUB -W 3:00
#BSUB -u josne@dtu.dk
#BSUB -B
#BSUB -N
#BSUB -o m2_fetch_lotus_%J.out
#BSUB -e m2_fetch_lotus_%J.err

# Fixes a failure in submit_host_fraction.sh: the NCBI FTP folder for
# GCA_012489685 (LjGifu_v1.2) holds no genome FASTA, so curl saved a 990-byte
# HTML error page and the Lotus index build failed silently. ENA serves the
# assembly (chromosomes 1-6 + MT + Pltd). This version fails loudly.
set -uo pipefail

WORK="/work3/josne/github/27221-Microbiome-Engineering_2027/Module2_MetaG/coassembly56/host_diag"
GEN="${WORK}/genomes/lotus_Gifu_v1.2.fna.gz"
mkdir -p "${WORK}/genomes" "${WORK}/index"
source /work3/josne/miniconda3/etc/profile.d/conda.sh
conda activate anvio-9

echo "Job started $(date), ID ${LSB_JOBID:-none}, node $(hostname)"
rm -f "${GEN}"
curl -fSL --retry 3 --retry-delay 15 -o "${GEN}" \
    "https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_012489685?download=true&gzip=true" \
    || { echo "ERROR: download failed"; exit 1; }
gzip -t "${GEN}" || { echo "ERROR: not a valid gzip"; exit 1; }
ls -la "${GEN}"

zcat "${GEN}" | grep '>' > "${WORK}/genomes/lotus_headers.txt"
echo "sequences: $(wc -l < ${WORK}/genomes/lotus_headers.txt)"
cat "${WORK}/genomes/lotus_headers.txt"
nbp=$(zcat "${GEN}" | grep -v '>' | tr -d '\n' | wc -c)
echo "total bp: ${nbp}"
[ "${nbp}" -gt 400000000 ] || { echo "ERROR: genome smaller than expected (${nbp} bp)"; exit 1; }

# The assembly already contains Pltd/MT, so no separate chloroplast file.
rm -f "${WORK}"/index/lotus.*
bowtie2-build --threads 8 "${GEN}" "${WORK}/index/lotus" > "${WORK}/index/lotus.build.log" 2>&1 \
    || { echo "ERROR: bowtie2-build failed"; tail -5 "${WORK}/index/lotus.build.log"; exit 1; }
ls -la "${WORK}"/index/lotus.*
echo "Job finished $(date)"
