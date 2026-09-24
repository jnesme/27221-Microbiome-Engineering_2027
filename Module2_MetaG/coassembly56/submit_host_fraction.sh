#!/bin/bash
### General options
### -- specify queue --
#BSUB -q hpc
### -- set the job Name --
#BSUB -J m2_host_fraction
### -- barley index build (~4.4 Gbp) is the heavy step; 16 cores, 96GB total --
#BSUB -n 16
#BSUB -R "span[hosts=1] rusage[mem=6GB]"
#BSUB -M 6500MB
### -- set walltime limit: hh:mm --
#BSUB -W 12:00
### -- set the email address --
#BSUB -u josne@dtu.dk
### -- send notification at start --
#BSUB -B
### -- send notification at completion --
#BSUB -N
### -- Specify the output and error file. %J is the job-id --
#BSUB -o m2_host_fraction_%J.out
#BSUB -e m2_host_fraction_%J.err

# Estimate the plant-host DNA fraction in representative reconstitution samples
# (Lotus / barley), to decide whether to deplete host reads before the
# 56-sample co-assembly. Subsamples ~500K read pairs per run and maps them to
# the host nuclear genome + chloroplast (organelle reads are host too), and to
# the Phase 1 bacterial reference for context.
set -uo pipefail

BASE="/work3/josne/github/27221-Microbiome-Engineering_2027/Module2_MetaG"
DIAG="${BASE}/coassembly56"
WORK="${DIAG}/host_diag"
RAW="${BASE}/raw_reads"
PHASE1_REF="${BASE}/session2_prep/REFERENCE-contigs-exported.fa"
NPAIRS=500000
THREADS=16

mkdir -p "${WORK}"/{genomes,index,subsample}
source /work3/josne/miniconda3/etc/profile.d/conda.sh
conda activate anvio-9

echo "=========================================="
echo "Module 2 - host DNA fraction diagnostic"
echo "Job started on $(date), ID: ${LSB_JOBID:-none}, node: $(hostname)"
echo "=========================================="

fetch_ncbi_genome() {   # accession outfile
    local acc=$1 out=$2
    [ -s "$out" ] && return 0
    local p="${acc:4:3}/${acc:7:3}/${acc:10:3}"
    local base="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/${p}/"
    local dir; dir=$(curl -s -m 60 "$base" | grep -oE 'GCA_[^"/]+' | head -1)
    [ -n "$dir" ] || { echo "ERROR: no NCBI dir for $acc"; return 1; }
    curl -sS --retry 3 -o "$out" "${base}${dir}/${dir}_genomic.fna.gz"
}
fetch_ena_seq() {       # accession outfile
    [ -s "$2" ] || curl -sS --retry 3 -o "$2" "https://www.ebi.ac.uk/ena/browser/api/fasta/$1"
}

echo "--- fetching host references ---"
fetch_ncbi_genome GCA_012489685 "${WORK}/genomes/lotus_Gifu_v1.2.fna.gz" || exit 1   # Lotus japonicus Gifu v1.2
fetch_ncbi_genome GCA_058768555 "${WORK}/genomes/barley_conrad_v1.0.fna.gz" || exit 1 # Hordeum vulgare (chromosome-level)
fetch_ena_seq AP002983 "${WORK}/genomes/lotus_chloroplast.fa"   # Lotus japonicus chloroplast
fetch_ena_seq EF115541 "${WORK}/genomes/barley_chloroplast.fa"  # barley (Morex) chloroplast
ls -la "${WORK}/genomes"

echo "--- building bowtie2 indexes ---"
[ -s "${WORK}/index/lotus.1.bt2" ]  || bowtie2-build --threads ${THREADS} \
    "${WORK}/genomes/lotus_Gifu_v1.2.fna.gz,${WORK}/genomes/lotus_chloroplast.fa" "${WORK}/index/lotus" > "${WORK}/index/lotus.build.log" 2>&1
[ -s "${WORK}/index/phase1.1.bt2" ] || bowtie2-build --threads ${THREADS} \
    "${PHASE1_REF}" "${WORK}/index/phase1" > "${WORK}/index/phase1.build.log" 2>&1
[ -s "${WORK}/index/barley.1.bt2l" ] || [ -s "${WORK}/index/barley.1.bt2" ] || bowtie2-build --large-index --threads ${THREADS} \
    "${WORK}/genomes/barley_conrad_v1.0.fna.gz,${WORK}/genomes/barley_chloroplast.fa" "${WORK}/index/barley" > "${WORK}/index/barley.build.log" 2>&1
ls "${WORK}/index" | head -30

echo "--- subsample and map ---"
OUT="${DIAG}/host_fraction.tsv"
printf "run_accession\tsample_alias\thost\tpairs\thost_pct\tphase1_ref_pct\n" > "${OUT}"
tail -n +2 "${DIAG}/host_diag_runs.tsv" | while IFS=$'\t' read -r run alias host; do
    sub1="${WORK}/subsample/${run}_1.fq"; sub2="${WORK}/subsample/${run}_2.fq"
    if [ ! -s "$sub1" ]; then
        zcat "${RAW}/${run}_1.fastq.gz" 2>/dev/null | head -n $((NPAIRS*4)) > "$sub1" || true
        zcat "${RAW}/${run}_2.fastq.gz" 2>/dev/null | head -n $((NPAIRS*4)) > "$sub2" || true
    fi
    pairs=$(( $(wc -l < "$sub1") / 4 ))
    if [ "$host" = "Lj" ]; then hidx="lotus"; else hidx="barley"; fi
    hpct=$(bowtie2 -p ${THREADS} -x "${WORK}/index/${hidx}" -1 "$sub1" -2 "$sub2" --no-unal -S /dev/null 2>&1 | grep "overall alignment rate" | awk '{print $1}')
    ppct=$(bowtie2 -p ${THREADS} -x "${WORK}/index/phase1" -1 "$sub1" -2 "$sub2" --no-unal -S /dev/null 2>&1 | grep "overall alignment rate" | awk '{print $1}')
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$run" "$alias" "$host" "$pairs" "$hpct" "$ppct" | tee -a "${OUT}"
done

echo "=========================================="
echo "Job finished on $(date)"
echo "Result table: ${OUT}"
cat "${OUT}"
echo "=========================================="
