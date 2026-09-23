#!/bin/bash
### General options
### -- specify queue --
#BSUB -q hpc
### -- set the job Name --
#BSUB -J m2_download
### -- download is network/IO-bound, not CPU-bound: a handful of concurrent
### curl transfers is plenty, no point reserving a full node --
#BSUB -n 4
#BSUB -R "span[hosts=1] rusage[mem=1GB]"
#BSUB -M 1200MB
### -- set walltime limit: hh:mm (22GB over the network, generous padding) --
#BSUB -W 2:00
### -- set the email address --
#BSUB -u josne@dtu.dk
### -- send notification at start --
#BSUB -B
### -- send notification at completion --
#BSUB -N
### -- Specify the output and error file. %J is the job-id --
#BSUB -o m2_download_%J.out
#BSUB -e m2_download_%J.err

MODULE2_DIR="/work3/josne/github/27221-Microbiome-Engineering_2027/Module2_MetaG"

echo "=========================================="
echo "Module 2 - Download sample subset"
echo "Job started on $(date)"
echo "Job ID: $LSB_JOBID"
echo "Running on node: $(hostname)"
echo "Manifest: ${MODULE2_DIR}/M2_sample_selection.tsv"
echo "Output:   ${MODULE2_DIR}/raw_reads"
echo "=========================================="

cd "${MODULE2_DIR}" || exit 1
/work3/josne/miniconda3/bin/python3 download_selection.py

EXIT_CODE=$?

echo "=========================================="
echo "Job finished on $(date)"
echo "Exit code: ${EXIT_CODE}"
echo "=========================================="

exit ${EXIT_CODE}
