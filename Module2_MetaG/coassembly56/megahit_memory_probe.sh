#!/bin/bash
# Measure MEGAHIT peak memory and time on a fraction of the depleted 56-sample
# dataset, to size the real co-assembly. Usage: megahit_memory_probe.sh PERCENT THREADS MEM_BYTES [MEM_FLAG=1]
# Takes the first PERCENT% of pairs of EVERY sample (keeps sample proportions),
# runs megahit with the same options as config_coasm56.json, reports peak RSS.
set -uo pipefail
P=$1; T=$2; MEM=$3; FLAG=${4:-1}
DIAG="/work3/josne/github/27221-Microbiome-Engineering_2027/Module2_MetaG/coassembly56"
WD="${DIAG}/megahit_probe/p${P}_f${FLAG}"
source /work3/josne/miniconda3/etc/profile.d/conda.sh
conda activate anvio-9
rm -rf "${WD}"; mkdir -p "${WD}/reads"
echo "probe ${P}% threads ${T} cap ${MEM} mem-flag ${FLAG} on $(hostname), $(date)"

R1=""; R2=""; total=0
while IFS=$'\t' read -r run alias mode n; do
    out=$(awk -F'\t' -v r="$run" '$1==r{print $7}' "${DIAG}/depletion_summary.tsv")
    K=$(( out * P / 100 ))
    seqkit head -n "$K" "${DIAG}/depleted_reads/${run}_1.fastq.gz" -o "${WD}/reads/${run}_1.fq.gz" 2>/dev/null
    seqkit head -n "$K" "${DIAG}/depleted_reads/${run}_2.fastq.gz" -o "${WD}/reads/${run}_2.fq.gz" 2>/dev/null
    R1="${R1:+$R1,}${WD}/reads/${run}_1.fq.gz"; R2="${R2:+$R2,}${WD}/reads/${run}_2.fq.gz"
    total=$(( total + K ))
done < "${DIAG}/depletion_samples.tsv"
echo "subset: ${total} pairs"

/usr/bin/time -v -o "${WD}/time.txt" \
    megahit -1 "$R1" -2 "$R2" -o "${WD}/asm" -t "$T" --min-contig-len 1000 --memory "$MEM" --mem-flag "$FLAG" > "${WD}/megahit.log" 2>&1
RC=$?
peak_kb=$(awk '/Maximum resident/{print $NF}' "${WD}/time.txt")
wall=$(awk '/Elapsed \(wall/{print $NF}' "${WD}/time.txt")
seqkit stats -T "${WD}/asm/final.contigs.fa" 2>/dev/null | tail -1 | cut -f4,5,13 > "${WD}/contig_stats.txt"
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$P" "$total" "$T" "$((peak_kb/1024/1024))" "$wall" "$(cat ${WD}/contig_stats.txt | tr '\t' ' ')" "$FLAG" >> "${DIAG}/megahit_probe/results.tsv"
echo "exit $RC, peak RSS $((peak_kb/1024)) MB, wall ${wall}"
rm -rf "${WD}/reads" "${WD}/asm/intermediate_contigs"
exit $RC
