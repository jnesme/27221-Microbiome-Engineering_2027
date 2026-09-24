#!/bin/bash
### General options
### -- specify queue --
#BSUB -q hpc
### -- set the job Name --
#BSUB -J m2s2_phase1_refassembly
### -- co-assembly of the 16 Input (SynCom-only, no host) samples --
#BSUB -n 14
#BSUB -R "span[hosts=1] rusage[mem=6GB]"
#BSUB -M 6500MB
### -- set walltime limit: hh:mm (generous ceiling, not an estimate - this
### hasn't been run before at this scale) --
#BSUB -W 24:00
### -- set the email address --
#BSUB -u josne@dtu.dk
### -- send notification at start --
#BSUB -B
### -- send notification at completion --
#BSUB -N
### -- Specify the output and error file. %J is the job-id --
#BSUB -o m2s2_phase1_%J.out
#BSUB -e m2s2_phase1_%J.err

PREP_DIR="/work3/josne/github/27221-Microbiome-Engineering_2027/Module2_MetaG/session2_prep"

echo "=========================================="
echo "Module 2 Session 2 prep - Phase 1: reference co-assembly"
echo "Job started on $(date)"
echo "Job ID: $LSB_JOBID"
echo "Running on node: $(hostname)"
echo "=========================================="

source /work3/josne/miniconda3/etc/profile.d/conda.sh
conda activate anvio-9

cd "${PREP_DIR}" || exit 1

anvi-run-workflow -w metagenomics -c config_phase1.json --skip-dry-run -A --jobs 14

EXIT_CODE=$?
if [ ${EXIT_CODE} -ne 0 ]; then
    echo "Phase 1 workflow failed with exit code ${EXIT_CODE} - not exporting contigs."
    exit ${EXIT_CODE}
fi

# The workflow's own reformatted-contigs fasta is a Snakemake temp() file
# that gets cleaned up - export a persistent copy from the final contigs
# database instead, for phase 2's fasta_txt to reference.
echo "=========================================="
echo "Exporting reference contigs fasta for phase 2"
echo "=========================================="
anvi-export-contigs \
    -c 03_CONTIGS/REFERENCE-contigs.db \
    -o "${PREP_DIR}/REFERENCE-contigs-exported.fa"

EXPORT_EXIT=$?

echo "=========================================="
echo "Job finished on $(date)"
echo "Workflow exit code: ${EXIT_CODE}, export exit code: ${EXPORT_EXIT}"
echo "=========================================="

exit ${EXPORT_EXIT}
