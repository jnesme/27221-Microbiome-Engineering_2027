#!/bin/bash
# Per-sample host-read depletion, sourced by submit_host_depletion.sh.
#
# Rule: a read PAIR is removed if EITHER mate aligns (bowtie2 end-to-end,
# default sensitivity) to a host NUCLEAR sequence. This is stricter than
# bowtie2 --un-conc, which keeps pairs where only one mate aligns.
# Organelle sequences are deliberately NOT depleted (Lotus index is
# chromosomes-only; alignments to the barley chloroplast contig are ignored),
# so organelle contigs still assemble and can be flagged as a teaching point.
#
# Requires exported: RAW, OUTDIR, TMPD, IDX (index dir), LOG, THREADS, SEQKIT_J
# deplete_one RUN MODE   where MODE is Lj | Hv | Input (Input = both hosts)

deplete_one() {
    local run=$1 mode=$2
    local r1="${RAW}/${run}_1.fastq.gz" r2="${RAW}/${run}_2.fastq.gz"
    local o1="${OUTDIR}/${run}_1.fastq.gz" o2="${OUTDIR}/${run}_2.fastq.gz"
    local names="${TMPD}/${run}.names" done="${TMPD}/${run}.done"

    if [ -e "$done" ] && [ -s "$o1" ] && [ -s "$o2" ]; then
        echo "[$(date +%T)] $run: already done, skipping" >> "${LOG}"; return 0
    fi
    local idxs
    case "$mode" in
        Lj)    idxs="lotus_nuclear" ;;
        Hv)    idxs="barley" ;;
        Input) idxs="lotus_nuclear barley" ;;
        *) echo "$run: unknown mode $mode" >> "${LOG}"; return 1 ;;
    esac

    : > "${names}.raw"
    local idx
    for idx in $idxs; do
        # names of pairs with an aligned mate; ignore the barley chloroplast contig
        bowtie2 -p "${THREADS}" --no-unal -x "${IDX}/${idx}" -1 "$r1" -2 "$r2" 2>> "${TMPD}/${run}.bowtie2.log" \
            | samtools view -F 4 - \
            | awk '$3 !~ /EF115541/ {print $1}' >> "${names}.raw" \
            || { echo "[$(date +%T)] $run: mapping to $idx FAILED" >> "${LOG}"; return 1; }
    done
    sort -u "${names}.raw" > "${names}"; rm -f "${names}.raw"
    local removed; removed=$(wc -l < "${names}")

    seqkit grep -j "${SEQKIT_J}" -v -f "${names}" "$r1" -o "$o1" \
        && seqkit grep -j "${SEQKIT_J}" -v -f "${names}" "$r2" -o "$o2" \
        || { echo "[$(date +%T)] $run: seqkit FAILED" >> "${LOG}"; rm -f "$o1" "$o2"; return 1; }

    echo "${removed}" > "${TMPD}/${run}.removed"
    touch "$done"
    echo "[$(date +%T)] $run ($mode): removed ${removed} pairs" >> "${LOG}"
}
