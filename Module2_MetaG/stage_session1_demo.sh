#!/bin/bash
# Stage the small, single-sample dataset used for Module 2 Session 1's live
# assembly demo. Reproducible from the already-downloaded full subset
# (Module2_MetaG/raw_reads/, see M2_sample_selection.tsv) - just symlinks the
# chosen sample under the R1/R2 naming M2_S1_reads_to_contigsDB.ipynb expects.
#
# Sample choice: Input_LjSC_B2 (SRR29710017) - the LjSC SynCom sequenced
# directly, no plant host. ~375K read pairs, ~52MB gzipped: a pure bacterial
# community, small enough to QC+assemble live in class, with no host-DNA
# noise complicating the first "what is a contigs database" demo. The
# host-selection biology (did the SynCom survive on the plant?) is Session
# 2's job, using the precomputed full profile.
#
# Note: a smaller candidate (Input_LjSC_B3 / SRR29710015, ~133K reads) was
# tried first and validated end-to-end - it ran without error but produced
# ZERO single-copy marker gene hits (too shallow for a ~177-member SynCom).
# SRR29710017 was chosen after confirming it actually assembles something
# real: 388 contigs, 9 Bacteria_71 + 2 16S rRNA HMM hits, ~4 min total
# compute (see Module2_MetaG/instructor_test_run/ for the validation run).
set -euo pipefail

MODULE2_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAMPLE="SRR29710017"

SRC_DIR="${MODULE2_DIR}/raw_reads"
DEST_DIR="${MODULE2_DIR}/raw_reads_subsampled"

if [ ! -f "${SRC_DIR}/${SAMPLE}_1.fastq.gz" ]; then
    echo "ERROR: ${SRC_DIR}/${SAMPLE}_1.fastq.gz not found."
    echo "Run submit_download.sh first to fetch the full M2_sample_selection.tsv subset."
    exit 1
fi

mkdir -p "${DEST_DIR}"
ln -sf "../raw_reads/${SAMPLE}_1.fastq.gz" "${DEST_DIR}/${SAMPLE}_R1.fastq.gz"
ln -sf "../raw_reads/${SAMPLE}_2.fastq.gz" "${DEST_DIR}/${SAMPLE}_R2.fastq.gz"

echo "Staged ${SAMPLE} into ${DEST_DIR}"
ls -la "${DEST_DIR}"
