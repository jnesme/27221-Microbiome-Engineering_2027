#!/bin/bash
### General options
### -- specify queue --
#BSUB -q hpc
### -- set the job Name --
#BSUB -J m2s2_binning_dastool
### -- Import SemiBin2's contig-level result as a collection, then reconcile
### CONCOCT + METABAT2 + SEMIBIN2 with DAS_Tool (the dastool driver pulls those
### three collections straight from the profile DB, exports contig2bin files,
### and shells out to the real DAS_Tool binary; see anvio/drivers/dastool.py).
### usearch (DAS_Tool's default single-copy-gene search engine) is not
### installed here, so --search-engine diamond is required. Submit this AFTER
### both submit_binning_semibin2.sh and submit_binning_concoct_metabat2.sh
### finish: bsub -w "ended(<semibin2 job>) && ended(<concoct/metabat2 job>)" --
#BSUB -n 8
#BSUB -R "span[hosts=1] rusage[mem=8GB]"
#BSUB -M 8192MB
### -- set walltime limit: hh:mm (DAS_Tool's single-copy-gene search over
### ~160k splits x 3 collections; generous ceiling, unmeasured) --
#BSUB -W 8:00
### -- set the email address --
#BSUB -u josne@dtu.dk
### -- send notification at start --
#BSUB -B
### -- send notification at completion --
#BSUB -N
### -- Specify the output and error file. %J is the job-id --
#BSUB -o m2s2_binning_dastool_%J.out
#BSUB -e m2s2_binning_dastool_%J.err

PREP_DIR="/work3/josne/github/27221-Microbiome-Engineering_2027/Module2_MetaG/session2_prep"
BIN_DIR="${PREP_DIR}/binning"
CONTIGS_DB="${PREP_DIR}/03_CONTIGS-coasm56/COASM56-contigs.db"
PROFILE_DB="${PREP_DIR}/06_MERGED-coasm56/COASM56/PROFILE.db"
SEMIBIN2_TSV="${BIN_DIR}/semibin2_contig2bin.tsv"

echo "=========================================="
echo "Module 2 Session 2 binning: import SemiBin2 + DAS_Tool consensus"
echo "Job started on $(date)"
echo "Job ID: $LSB_JOBID"
echo "Running on node: $(hostname)"
echo "=========================================="

source /work3/josne/miniconda3/etc/profile.d/conda.sh
conda activate anvio-9

if [ ! -s "${SEMIBIN2_TSV}" ]; then
    echo "ERROR: ${SEMIBIN2_TSV} missing or empty - did submit_binning_semibin2.sh finish successfully?"
    exit 1
fi

echo "--- import SemiBin2 collection (contigs-mode) ---"
anvi-import-collection "${SEMIBIN2_TSV}" -c "${CONTIGS_DB}" -p "${PROFILE_DB}" \
    -C SEMIBIN2 --contigs-mode
IMPORT_EXIT=$?
echo "import exit: ${IMPORT_EXIT}"
if [ ${IMPORT_EXIT} -ne 0 ]; then
    exit ${IMPORT_EXIT}
fi

echo "--- DAS_Tool consensus over CONCOCT, METABAT2, SEMIBIN2 ---"
anvi-cluster-contigs -p "${PROFILE_DB}" -c "${CONTIGS_DB}" \
    -C DASTOOL_CONSENSUS --driver dastool \
    -S "CONCOCT,METABAT2,SEMIBIN2" --search-engine diamond -T 8 --just-do-it
DASTOOL_EXIT=$?
echo "DAS_Tool exit: ${DASTOOL_EXIT}"

echo "=========================================="
echo "--- all collections now in the profile DB ---"
anvi-show-collections-and-bins -p "${PROFILE_DB}" 2>&1
echo "=========================================="
echo "Job finished on $(date)"
echo "import exit: ${IMPORT_EXIT}, DAS_Tool exit: ${DASTOOL_EXIT}"
echo "=========================================="

exit ${DASTOOL_EXIT}
