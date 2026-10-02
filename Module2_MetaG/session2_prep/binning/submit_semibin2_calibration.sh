#!/bin/bash
### General options
### -- specify queue --
#BSUB -q hpc
### -- set the job Name --
#BSUB -J m2s2_semibin2_calib
### -- calibration run only: --epochs 1 instead of the default 20, to measure
### real wall time and peak memory on this dataset before sizing the full
### SemiBin2 job (same approach used for megahit's memory, see
### coassembly56/megahit_memory_probe.sh). 16 threads, modest memory: SemiBin2's
### footprint scales with contigs/samples, not raw read volume, so this should
### be well under 100 GB; this run itself is the test of that. --
#BSUB -n 16
#BSUB -R "span[hosts=1] rusage[mem=6GB]"
#BSUB -M 6144MB
### -- set walltime limit: hh:mm (one epoch should be a small fraction of a
### full 20-epoch run; generous ceiling since this is unmeasured) --
#BSUB -W 4:00
### -- set the email address --
#BSUB -u josne@dtu.dk
### -- send notification at start --
#BSUB -B
### -- send notification at completion --
#BSUB -N
### -- Specify the output and error file. %J is the job-id --
#BSUB -o m2s2_semibin2_calib_%J.out
#BSUB -e m2s2_semibin2_calib_%J.err

PREP_DIR="/work3/josne/github/27221-Microbiome-Engineering_2027/Module2_MetaG/session2_prep"
BIN_DIR="${PREP_DIR}/binning"
MAP_DIR="${PREP_DIR}/04_MAPPING-coasm56/COASM56"
OUT="${BIN_DIR}/semibin2_calib"

echo "=========================================="
echo "Module 2 Session 2 binning: SemiBin2 calibration (1 epoch)"
echo "Job started on $(date)"
echo "Job ID: $LSB_JOBID"
echo "Running on node: $(hostname)"
echo "=========================================="

source /work3/josne/miniconda3/etc/profile.d/conda.sh
conda activate semibin2

mkdir -p "${OUT}"
BAMS=$(ls "${MAP_DIR}"/*.bam)
echo "Using $(echo "${BAMS}" | wc -w) BAM files from ${MAP_DIR}"

/usr/bin/time -v -o "${BIN_DIR}/semibin2_calib_time.txt" \
    SemiBin2 single_easy_bin \
        -i "${PREP_DIR}/COASM56-contigs-exported.fa" \
        -b ${BAMS} \
        -o "${OUT}" \
        -p 16 \
        --epochs 1 \
        > "${BIN_DIR}/semibin2_calib.log" 2>&1

EXIT_CODE=$?

echo "=========================================="
echo "Job finished on $(date)"
echo "Exit code: ${EXIT_CODE}"
echo "--- time -v summary ---"
grep -E "Elapsed|Maximum resident|Percent of CPU" "${BIN_DIR}/semibin2_calib_time.txt" 2>/dev/null
echo "=========================================="

exit ${EXIT_CODE}
