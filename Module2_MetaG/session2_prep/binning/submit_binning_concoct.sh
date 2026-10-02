#!/bin/bash
### General options
### -- specify queue --
#BSUB -q hpc
### -- set the job Name --
#BSUB -J m2s2_binning_concoct
### -- CONCOCT via anvi-cluster-contigs, on the merged 56-sample profile
### (160,325 splits). Split out of the original combined CONCOCT+MetaBAT2 job:
### CONCOCT first failed there with a numpy ABI mismatch in the anvio-9 env
### (its compiled vbgmm extension was built against numpy 1.x, but something
### had pip-installed SemiBin 2.5.0 directly into anvio-9 too, which requires
### numpy>=2 and silently bumped it there - unrelated to the dedicated
### semibin2 env used for the real binning job). Fixed by reinstalling
### concoct+biopython via conda (conda-forge/bioconda), which downgraded
### numpy back to 1.26.4; verified `concoct --version` runs clean afterwards.
### Submit this ONLY after submit_binning_metabat2.sh has ended (bsub -w
### "ended(<that job id>)"), never concurrently with it: anvi'o opens the
### profile DB with plain sqlite3.connect() (anvio/db.py), no WAL/busy-timeout
### tuning, so two anvi-cluster-contigs processes writing a new collection
### into the same PROFILE.db at once risk "database is locked". --
#BSUB -n 8
#BSUB -R "span[hosts=1] rusage[mem=8GB]"
#BSUB -M 8192MB
### -- set walltime limit: hh:mm (CONCOCT's VGMM fit is largely
### single-threaded and untested at this split count; generous ceiling,
### not an estimate) --
#BSUB -W 24:00
### -- set the email address --
#BSUB -u josne@dtu.dk
### -- send notification at start --
#BSUB -B
### -- send notification at completion --
#BSUB -N
### -- Specify the output and error file. %J is the job-id --
#BSUB -o m2s2_binning_concoct_%J.out
#BSUB -e m2s2_binning_concoct_%J.err

PREP_DIR="/work3/josne/github/27221-Microbiome-Engineering_2027/Module2_MetaG/session2_prep"
CONTIGS_DB="${PREP_DIR}/03_CONTIGS-coasm56/COASM56-contigs.db"
PROFILE_DB="${PREP_DIR}/06_MERGED-coasm56/COASM56/PROFILE.db"

echo "=========================================="
echo "Module 2 Session 2 binning: CONCOCT"
echo "Job started on $(date)"
echo "Job ID: $LSB_JOBID"
echo "Running on node: $(hostname)"
echo "=========================================="

source /work3/josne/miniconda3/etc/profile.d/conda.sh
conda activate anvio-9

/usr/bin/time -v -o "${PREP_DIR}/binning/concoct_time.txt" \
    anvi-cluster-contigs -p "${PROFILE_DB}" -c "${CONTIGS_DB}" \
        -C CONCOCT --driver concoct -T 8 --just-do-it
EXIT_CODE=$?

echo "=========================================="
echo "Job finished on $(date)"
echo "Exit code: ${EXIT_CODE}"
echo "--- time -v summary ---"
grep -E "Elapsed|Maximum resident|Percent of CPU" "${PREP_DIR}/binning/concoct_time.txt" 2>/dev/null
echo "=========================================="

exit ${EXIT_CODE}
