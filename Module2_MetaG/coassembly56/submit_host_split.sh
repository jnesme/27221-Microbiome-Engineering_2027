#!/bin/bash
### General options
### -- specify queue --
#BSUB -q hpc
### -- set the job Name --
#BSUB -J m2_host_split
### -- runs only after the host-fraction job (29481226) and the Lotus fix job (29481885) --
#BSUB -w "ended(29481226) && ended(29481885)"
### -- mapping ~1M reads x 6 samples: light. 8 cores, 64GB total --
#BSUB -n 8
#BSUB -R "span[hosts=1] rusage[mem=8GB]"
#BSUB -M 8500MB
### -- set walltime limit: hh:mm --
#BSUB -W 4:00
### -- set the email address --
#BSUB -u josne@dtu.dk
### -- send notification at start --
#BSUB -B
### -- send notification at completion --
#BSUB -N
### -- Specify the output and error file. %J is the job-id --
#BSUB -o m2_host_split_%J.out
#BSUB -e m2_host_split_%J.err

# Follow-up to submit_host_fraction.sh. Splits host-mapping reads into
# nuclear vs plastid, and measures how many reads map BOTH to the host and to
# the bacterial Phase 1 reference (the plastid-homology concern): those are the
# reads whose origin is ambiguous. Reuses the indexes and subsamples built by
# job 29481226.
set -uo pipefail

DIAG="/work3/josne/github/27221-Microbiome-Engineering_2027/Module2_MetaG/coassembly56"
WORK="${DIAG}/host_diag"
TMP="${WORK}/split_tmp"
THREADS=8
# organelle sequence names: Lotus assembly's own chloroplast+mitochondrion; barley chloroplast (no barley mito in the index)
PL_RX_LJ='AP022636|AP022637'
PL_RX_HV='EF115541'

mkdir -p "${TMP}"
source /work3/josne/miniconda3/etc/profile.d/conda.sh
conda activate anvio-9

echo "Job started on $(date), ID: ${LSB_JOBID:-none}, node: $(hostname)"

OUT="${DIAG}/host_split.tsv"
printf "run_accession\tsample_alias\thost\treads\tnuclear_pct\torganelle_pct\tphase1ref_pct\torganelle_and_phase1\tnuclear_and_phase1\torganelle_and_phase1_pct_of_organelle\tnuclear_and_phase1_pct_of_nuclear\n" > "${OUT}"

tail -n +2 "${DIAG}/host_diag_runs.tsv" | while IFS=$'\t' read -r run alias host; do
    s1="${WORK}/subsample/${run}_1.fq"; s2="${WORK}/subsample/${run}_2.fq"
    [ -s "$s1" ] || { echo "missing subsample $s1"; continue; }
    reads=$(( $(wc -l < "$s1") / 4 * 2 ))
    if [ "$host" = "Lj" ]; then hidx="lotus"; PL_RX="${PL_RX_LJ}"; else hidx="barley"; PL_RX="${PL_RX_HV}"; fi

    # host mapping -> plastid/other key lists
    bowtie2 -p ${THREADS} -x "${WORK}/index/${hidx}" -1 "$s1" -2 "$s2" --no-unal 2>/dev/null \
        | python3 "${DIAG}/split_sam.py" --prefix "${TMP}/${run}.host" --plastid-regex "${PL_RX}" > "${TMP}/${run}.host.counts"
    # bacterial (Phase 1 reference) mapping -> all keys in 'other'
    bowtie2 -p ${THREADS} -x "${WORK}/index/phase1" -1 "$s1" -2 "$s2" --no-unal 2>/dev/null \
        | python3 "${DIAG}/split_sam.py" --prefix "${TMP}/${run}.bact" > "${TMP}/${run}.bact.counts"

    np=$(awk '$1=="plastid"{print $2}' "${TMP}/${run}.host.counts")
    nn=$(awk '$1=="other"{print $2}'   "${TMP}/${run}.host.counts")
    nb=$(awk '$1=="other"{print $2}'   "${TMP}/${run}.bact.counts")
    sort -u "${TMP}/${run}.host.plastid.keys" > "${TMP}/${run}.p.sorted"
    sort -u "${TMP}/${run}.host.other.keys"   > "${TMP}/${run}.n.sorted"
    sort -u "${TMP}/${run}.bact.other.keys"   > "${TMP}/${run}.b.sorted"
    pb=$(comm -12 "${TMP}/${run}.p.sorted" "${TMP}/${run}.b.sorted" | wc -l)
    nbov=$(comm -12 "${TMP}/${run}.n.sorted" "${TMP}/${run}.b.sorted" | wc -l)

    awk -v r="$run" -v a="$alias" -v h="$host" -v R="$reads" -v np="$np" -v nn="$nn" -v nb="$nb" -v pb="$pb" -v nbov="$nbov" 'BEGIN{
        printf "%s\t%s\t%s\t%d\t%.2f\t%.2f\t%.2f\t%d\t%d\t%.2f\t%.2f\n", r,a,h,R,100*nn/R,100*np/R,100*nb/R,pb,nbov,(np>0?100*pb/np:0),(nn>0?100*nbov/nn:0)}' | tee -a "${OUT}"
done

echo "Job finished on $(date)"
cat "${OUT}"
