#!/bin/bash
### General options
### -- specify queue --
#BSUB -q hpc
### -- set the job Name --
#BSUB -J m2s2_binning_metabat2_direct
### -- Run MetaBAT2 directly (bypassing anvi'o's `--driver metabat2`, which is
### broken against this MetaBAT2 version: it hardcodes `-l`/`--onlyLabel`
### output and then parses lines with `l.split("\t")[0]`, but confirmed by a
### controlled toy-data reproduction that `-l` output has NO tab character at
### all - every line is a bare contig name, so the parser returns the whole
### line including its trailing "\n" and the lookup fails for every contig,
### 100% of the time, not data-dependent). This script uses the DEFAULT
### output mode instead (real per-bin FASTA files, no -l).
### Depth source: NOT jgi_summarize_bam_contig_depths (binning/metabat2_depth.txt) -
### confirmed, via direct samtools cross-checks on several contigs, to produce
### wildly wrong (inflated and negative) depth values for this dataset; an
### independent Biostars report corroborates this as a known Bioconda-build
### bug, not present in the author's own build. Using
### binning/metabat2_depth_anvio.txt instead (export_anvio_coverage_for_metabat2.py,
### the SAME per-contig coverage anvi'o's own CONCOCT driver already used
### successfully - verified to match real samtools depth exactly on the
### contigs jgi_summarize got wrong), with --cvExt (a coverage file without
### variance; verified format empirically on toy data: contigName, contigLen,
### one column per sample, no header requirement issue, no totalAvgDepth). --
#BSUB -n 8
#BSUB -R "span[hosts=1] rusage[mem=8GB]"
#BSUB -M 8192MB
### -- set walltime limit: hh:mm (the depth file is already made, the slow
### part; metabat2's own clustering on 158k contigs took a comparable order
### of magnitude inside the failed 18-min job - generous ceiling, not an
### estimate) --
#BSUB -W 6:00
### -- set the email address --
#BSUB -u josne@dtu.dk
### -- send notification at start --
#BSUB -B
### -- send notification at completion --
#BSUB -N
### -- Specify the output and error file. %J is the job-id --
#BSUB -o m2s2_binning_metabat2_direct_%J.out
#BSUB -e m2s2_binning_metabat2_direct_%J.err

PREP_DIR="/work3/josne/github/27221-Microbiome-Engineering_2027/Module2_MetaG/session2_prep"
BIN_DIR="${PREP_DIR}/binning"
OUT="${BIN_DIR}/metabat2_full"
CONTIGS_DB="${PREP_DIR}/03_CONTIGS-coasm56/COASM56-contigs.db"
PROFILE_DB="${PREP_DIR}/06_MERGED-coasm56/COASM56/PROFILE.db"

echo "=========================================="
echo "Module 2 Session 2 binning: MetaBAT2 (direct, default FASTA output)"
echo "Job started on $(date)"
echo "Job ID: $LSB_JOBID"
echo "Running on node: $(hostname)"
echo "=========================================="

source /work3/josne/miniconda3/etc/profile.d/conda.sh
conda activate anvio-9

mkdir -p "${OUT}"

/usr/bin/time -v -o "${BIN_DIR}/metabat2_direct_time.txt" \
    metabat2 -i "${PREP_DIR}/COASM56-contigs-exported.fa" \
             -a "${BIN_DIR}/metabat2_depth_anvio.txt" \
             --cvExt \
             -o "${OUT}/bin" \
             -t 8 \
             > "${BIN_DIR}/metabat2_direct.log" 2>&1
EXIT_CODE=$?
echo "metabat2 exit: ${EXIT_CODE}"

if [ ${EXIT_CODE} -eq 0 ]; then
    python3 - "${OUT}" "${BIN_DIR}/metabat2_contig2bin.tsv" <<'PYEOF'
import sys, glob, os
from Bio import SeqIO
outdir, dest = sys.argv[1], sys.argv[2]
with open(dest, "w") as f:
    for fa in sorted(glob.glob(os.path.join(outdir, "bin.*.fa"))):
        # anvi'o rejects bin names with characters outside [A-Za-z0-9_] - the
        # default "bin.1.fa" naming leaves a literal "." otherwise.
        bin_name = "MetaBAT2_" + os.path.basename(fa).rsplit(".fa", 1)[0].replace(".", "_")
        for rec in SeqIO.parse(fa, "fasta"):
            f.write(f"{rec.id}\t{bin_name}\n")
print(f"wrote {dest}")
PYEOF
    echo "contig2bin export exit: $?"
    wc -l "${BIN_DIR}/metabat2_contig2bin.tsv"

    echo "--- import METABAT2 collection (contigs-mode) ---"
    anvi-import-collection "${BIN_DIR}/metabat2_contig2bin.tsv" -c "${CONTIGS_DB}" -p "${PROFILE_DB}" \
        -C METABAT2 --contigs-mode
    IMPORT_EXIT=$?
    echo "import exit: ${IMPORT_EXIT}"
else
    IMPORT_EXIT=1
fi

echo "=========================================="
echo "Job finished on $(date)"
echo "metabat2 exit: ${EXIT_CODE}, import exit: ${IMPORT_EXIT}"
echo "--- time -v summary ---"
grep -E "Elapsed|Maximum resident|Percent of CPU" "${BIN_DIR}/metabat2_direct_time.txt" 2>/dev/null
echo "=========================================="

exit ${IMPORT_EXIT}
