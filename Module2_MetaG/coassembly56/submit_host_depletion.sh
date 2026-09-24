#!/bin/bash
### General options
#BSUB -q hpc
#BSUB -J m2_host_depletion
### -- 3 samples in parallel x 8 bowtie2 threads = 24 slots. Memory: each
### barley bowtie2 process holds ~7GB of index; with name sets/sorting ~30GB peak,
### reserved 3GB x 24 = 72GB (2x headroom). --
#BSUB -n 24
#BSUB -R "span[hosts=1] rusage[mem=3GB]"
#BSUB -M 3500MB
### -- ~440M read-mappings (Input samples map to two genomes) plus gzip
### filtering of 22GB: estimated 2-5h; 8h ceiling. --
#BSUB -W 8:00
#BSUB -u josne@dtu.dk
#BSUB -B
#BSUB -N
#BSUB -o m2_host_depletion_%J.out
#BSUB -e m2_host_depletion_%J.err

# Remove plant host NUCLEAR reads (Lotus / barley) from all 56 samples before
# the co-assembly. Reconstitution samples: their own host; Input samples: both
# hosts (they carry trace host DNA). Organelle reads are kept on purpose.
# See host_depletion_lib.sh for the removal rule. Resumable (skips finished runs).
set -uo pipefail

BASE="/work3/josne/github/27221-Microbiome-Engineering_2027/Module2_MetaG"
DIAG="${BASE}/coassembly56"
export RAW="${BASE}/raw_reads"
export OUTDIR="${DIAG}/depleted_reads"
export TMPD="${DIAG}/depletion_tmp"
export IDX="${DIAG}/host_diag/index"
export LOG="${DIAG}/depletion_progress.log"
export THREADS=8
export SEQKIT_J=4
NPAR=3

mkdir -p "${OUTDIR}" "${TMPD}"
source /work3/josne/miniconda3/etc/profile.d/conda.sh
conda activate anvio-9
source "${DIAG}/host_depletion_lib.sh"

echo "Job started $(date), ID ${LSB_JOBID:-none}, node $(hostname)"

# --- Lotus nuclear-only index (chromosomes only; drops the assembly's chloroplast + mitochondrion)
if [ ! -s "${IDX}/lotus_nuclear.1.bt2" ]; then
    zcat "${DIAG}/host_diag/genomes/lotus_Gifu_v1.2.fna.gz" \
        | seqkit grep -n -r -p "chromosome" -o "${DIAG}/host_diag/genomes/lotus_nuclear.fa.gz"
    nseq=$(zcat "${DIAG}/host_diag/genomes/lotus_nuclear.fa.gz" | grep -c '>')
    echo "Lotus nuclear sequences: ${nseq} (expected 6)"
    [ "${nseq}" -eq 6 ] || { echo "ERROR: unexpected sequence count"; exit 1; }
    bowtie2-build --threads 24 "${DIAG}/host_diag/genomes/lotus_nuclear.fa.gz" "${IDX}/lotus_nuclear" \
        > "${IDX}/lotus_nuclear.build.log" 2>&1 || { echo "ERROR: lotus index build failed"; exit 1; }
fi
[ -s "${IDX}/barley.1.bt2l" ] || { echo "ERROR: barley index missing"; exit 1; }

# --- sample list: biggest first for load balance
python3 - <<'EOF'
import csv
base="/work3/josne/github/27221-Microbiome-Engineering_2027/Module2_MetaG"
rows=list(csv.DictReader(open(f"{base}/M2_sample_selection.tsv"),delimiter="\t"))
out=[]
for r in rows:
    mode = "Input" if r["sample_group"].startswith("Input") else r["host"]
    out.append((int(r["read_count"]), r["run_accession"], r["sample_alias"], mode))
out.sort(reverse=True)
with open(f"{base}/coassembly56/depletion_samples.tsv","w") as f:
    for n,run,alias,mode in out: f.write(f"{run}\t{alias}\t{mode}\t{n}\n")
print(len(out),"samples")
EOF

export -f deplete_one
cut -f1,3 "${DIAG}/depletion_samples.tsv" \
    | xargs -P ${NPAR} -L 1 bash -c 'deplete_one "$0" "$1"'
echo "xargs exit: $?"

# --- summary + verification: output pairs must equal input - removed, mates in sync
echo "verifying outputs..."
seqkit stats -T -j 24 "${OUTDIR}"/*_1.fastq.gz "${OUTDIR}"/*_2.fastq.gz > "${DIAG}/depletion_out_stats.tsv" 2>/dev/null
python3 - <<'EOF'
import csv,os,sys
base="/work3/josne/github/27221-Microbiome-Engineering_2027/Module2_MetaG/coassembly56"
stats={r["file"].split("/")[-1]:int(r["num_seqs"].replace(",","")) for r in csv.DictReader(open(f"{base}/depletion_out_stats.tsv"),delimiter="\t")}
bad=0
with open(f"{base}/depletion_summary.tsv","w") as f:
    f.write("run\talias\tmode\tpairs_in\tpairs_removed\tpct_removed\tpairs_out\tstatus\n")
    for line in open(f"{base}/depletion_samples.tsv"):
        run,alias,mode,n=line.rstrip("\n").split("\t"); n=int(n)
        rp=f"{base}/depletion_tmp/{run}.removed"
        if not os.path.exists(rp): f.write(f"{run}\t{alias}\t{mode}\t{n}\t\t\t\tMISSING\n"); bad+=1; continue
        rem=int(open(rp).read()); o1=stats.get(f"{run}_1.fastq.gz"); o2=stats.get(f"{run}_2.fastq.gz")
        ok = (o1==o2==n-rem)
        if not ok: bad+=1
        f.write(f"{run}\t{alias}\t{mode}\t{n}\t{rem}\t{100*rem/n:.2f}\t{o1}\t{'OK' if ok else 'MISMATCH'}\n")
print("samples with problems:",bad)
sys.exit(1 if bad else 0)
EOF
RC=$?
echo "Job finished $(date); summary: ${DIAG}/depletion_summary.tsv; verification exit ${RC}"
cat "${DIAG}/depletion_summary.tsv"
exit ${RC}
