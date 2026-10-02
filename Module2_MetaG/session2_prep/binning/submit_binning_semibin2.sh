#!/bin/bash
### General options
### -- specify queue --
#BSUB -q hpc
### -- set the job Name --
#BSUB -J m2s2_binning_semibin2
### -- Full SemiBin2 run (single_easy_bin: "single or co-assembly" mode -
### COASM56 is one co-assembly, not per-sample assemblies, so this is the
### right subcommand, not multi_easy_bin), default 20 epochs, using the
### exported contigs FASTA + all 56 BAMs for cross-sample coverage. CPU-only
### PyTorch in the semibin2 env (no CUDA) - training runs on CPU.
### Sized from the calibration job (29541293, 1 epoch, exit 0): log timestamps
### split the 24m28s wall time into stages: ~6m24s fixed setup (56-BAM
### depth/composition generation, independent of epoch count), ~3m31s for the
### 1 training epoch, and ~14m21s binning+reclustering (mostly the naive ORF
### finder - the single biggest stage, also ~epoch-independent). So full run
### (20 epochs) ~= 6m24s + 20*3m31s + 14m21s =~ 1h31m; walltime below has >3x
### headroom. Peak RSS was 18.6 GB and shouldn't grow much with epoch count
### (same data/model size, just more iterations) - memory below has ~1.7x
### headroom. Only ~4.5 of the 16 requested cores were used on average
### (455% CPU) - several stages (ORF finder, reclustering) are not
### highly parallel - kept at 16 anyway since the BAM-processing stage does
### use them and the headroom costs little.
### Runs concurrently with submit_binning_metabat2.sh: this script never
### touches the anvi'o profile/contigs DBs, so there is no write conflict
### with that job (or, once reinstalled, submit_binning_concoct.sh). --
#BSUB -n 16
#BSUB -R "span[hosts=1] rusage[mem=2GB]"
#BSUB -M 2048MB
### -- set walltime limit: hh:mm - measured full-run estimate ~1h31m, >3x
### headroom --
#BSUB -W 5:00
### -- set the email address --
#BSUB -u josne@dtu.dk
### -- send notification at start --
#BSUB -B
### -- send notification at completion --
#BSUB -N
### -- Specify the output and error file. %J is the job-id --
#BSUB -o m2s2_binning_semibin2_%J.out
#BSUB -e m2s2_binning_semibin2_%J.err

PREP_DIR="/work3/josne/github/27221-Microbiome-Engineering_2027/Module2_MetaG/session2_prep"
BIN_DIR="${PREP_DIR}/binning"
MAP_DIR="${PREP_DIR}/04_MAPPING-coasm56/COASM56"
OUT="${BIN_DIR}/semibin2_full"

echo "=========================================="
echo "Module 2 Session 2 binning: SemiBin2 (full run, single_easy_bin)"
echo "Job started on $(date)"
echo "Job ID: $LSB_JOBID"
echo "Running on node: $(hostname)"
echo "=========================================="

source /work3/josne/miniconda3/etc/profile.d/conda.sh
conda activate semibin2

mkdir -p "${OUT}"
BAMS=$(ls "${MAP_DIR}"/*.bam)
echo "Using $(echo "${BAMS}" | wc -w) BAM files from ${MAP_DIR}"

/usr/bin/time -v -o "${BIN_DIR}/semibin2_full_time.txt" \
    SemiBin2 single_easy_bin \
        -i "${PREP_DIR}/COASM56-contigs-exported.fa" \
        -b ${BAMS} \
        -o "${OUT}" \
        -p 16 \
        > "${BIN_DIR}/semibin2_full.log" 2>&1

EXIT_CODE=$?
CONVERT_EXIT=0

if [ ${EXIT_CODE} -eq 0 ]; then
    # SemiBin2 writes its own contig->bin assignment table directly
    # (<out>/contig_bins.tsv: "contig\tbin", bin "-1" = unbinned) - no need to
    # parse the bin FASTAs (that earlier approach used Biopython, which is
    # NOT installed in the semibin2 env and silently failed with exit 1 on
    # the first full run; this file is both correct and simpler).
    awk -F'\t' 'NR>1 && $2!="-1" {print $1"\tSemiBin2_"$2}' \
        "${OUT}/contig_bins.tsv" > "${BIN_DIR}/semibin2_contig2bin.tsv"
    CONVERT_EXIT=$?
    echo "contig2bin export exit: ${CONVERT_EXIT}"
    wc -l "${BIN_DIR}/semibin2_contig2bin.tsv"
fi

echo "=========================================="
echo "Job finished on $(date)"
echo "SemiBin2 exit code: ${EXIT_CODE}, conversion exit code: ${CONVERT_EXIT}"
echo "--- time -v summary ---"
grep -E "Elapsed|Maximum resident|Percent of CPU" "${BIN_DIR}/semibin2_full_time.txt" 2>/dev/null
echo "=========================================="

# Previously this masked a failed conversion step behind SemiBin2's own
# (successful) exit code - report the real combined status.
if [ ${EXIT_CODE} -ne 0 ]; then
    exit ${EXIT_CODE}
fi
exit ${CONVERT_EXIT}
