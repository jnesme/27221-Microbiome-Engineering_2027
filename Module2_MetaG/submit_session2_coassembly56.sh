#!/bin/bash
### General options
### -- specify queue --
#BSUB -q hpc
### -- set the job Name --
#BSUB -J m2s2_coasm56
### -- co-assembly of ALL 56 host-depleted samples, then mapping, profiling and
### merging of all 56 against it. One node: megahit is pinned to 600 GB
### (--memory in config_coasm56.json), 32 slots x 22 GB = 704 GB reserved. --
#BSUB -n 32
#BSUB -R "span[hosts=1] rusage[mem=22GB]"
#BSUB -M 22528MB
### -- 72h is the queue maximum. Submit this script TWICE: the second time with
### -w "ended(<first job id>)" so that a walltime kill is resumed by Snakemake
### (it is a no-op if the first run finished). --
#BSUB -W 72:00
### -- set the email address --
#BSUB -u josne@dtu.dk
### -- send notification at start --
#BSUB -B
### -- send notification at completion --
#BSUB -N
### -- Specify the output and error file. %J is the job-id --
#BSUB -o m2s2_coasm56_%J.out
#BSUB -e m2s2_coasm56_%J.err

PREP_DIR="/work3/josne/github/27221-Microbiome-Engineering_2027/Module2_MetaG/session2_prep"

echo "=========================================="
echo "Module 2 Session 2 prep: co-assembly of all 56 host-depleted samples"
echo "Job started on $(date)"
echo "Job ID: $LSB_JOBID"
echo "Running on node: $(hostname)"
echo "=========================================="

source /work3/josne/miniconda3/etc/profile.d/conda.sh
conda activate anvio-9

cd "${PREP_DIR}" || exit 1

# A previous run killed at the walltime leaves Snakemake locks behind
if [ -n "$(ls -A .snakemake/locks 2>/dev/null)" ]; then
    echo "Stale Snakemake locks found, unlocking"
    anvi-run-workflow -w metagenomics -c config_coasm56.json --skip-dry-run -A --unlock --jobs 1
fi

anvi-run-workflow -w metagenomics -c config_coasm56.json --skip-dry-run -A --jobs 32

EXIT_CODE=$?
if [ ${EXIT_CODE} -ne 0 ]; then
    echo "Workflow failed or was incomplete (exit ${EXIT_CODE}); not exporting contigs."
    exit ${EXIT_CODE}
fi

# The workflow's reformatted-contigs fasta is a Snakemake temp() file: export a
# persistent copy and record assembly statistics.
anvi-export-contigs -c 03_CONTIGS-coasm56/COASM56-contigs.db -o "${PREP_DIR}/COASM56-contigs-exported.fa"
EXPORT_EXIT=$?
seqkit stats -a "${PREP_DIR}/COASM56-contigs-exported.fa"

echo "=========================================="
echo "Job finished on $(date)"
echo "Workflow exit code: ${EXIT_CODE}, export exit code: ${EXPORT_EXIT}"
echo "=========================================="
exit ${EXPORT_EXIT}
