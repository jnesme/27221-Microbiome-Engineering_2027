#!/bin/bash
# Verify host depletion outputs and write depletion_summary.tsv.
# Input pair counts come from the raw FASTQ files themselves: the read_count in
# M2_sample_selection.tsv (ENA metadata) is off by a few pairs for some runs.
# A sample is OK if R1 and R2 raw counts agree, R1 and R2 output counts agree,
# and pairs_out = pairs_in - pairs_removed. Exit 1 if any sample is not OK.
# Usage: verify_depletion.sh [threads]
set -uo pipefail
J=${1:-8}
BASE="/work3/josne/github/27221-Microbiome-Engineering_2027/Module2_MetaG"
DIAG="${BASE}/coassembly56"
source /work3/josne/miniconda3/etc/profile.d/conda.sh
conda activate anvio-9

seqkit stats -T -j "$J" "${BASE}"/raw_reads/SRR*_[12].fastq.gz > "${DIAG}/raw_pair_counts.tsv" 2>/dev/null
seqkit stats -T -j "$J" "${DIAG}"/depleted_reads/*_1.fastq.gz "${DIAG}"/depleted_reads/*_2.fastq.gz > "${DIAG}/depletion_out_stats.tsv" 2>/dev/null

python3 - <<'PY'
import csv,os,sys
base="/work3/josne/github/27221-Microbiome-Engineering_2027/Module2_MetaG/coassembly56"
def load(p): return {r["file"].split("/")[-1]:int(r["num_seqs"].replace(",","")) for r in csv.DictReader(open(p),delimiter="\t")}
raw=load(f"{base}/raw_pair_counts.tsv"); out=load(f"{base}/depletion_out_stats.tsv")
bad=0
with open(f"{base}/depletion_summary.tsv","w") as f:
    f.write("run\talias\tmode\tpairs_in\tpairs_removed\tpct_removed\tpairs_out\tstatus\n")
    for line in open(f"{base}/depletion_samples.tsv"):
        run,alias,mode,_=line.rstrip("\n").split("\t")
        rp=f"{base}/depletion_tmp/{run}.removed"
        i1,i2=raw.get(f"{run}_1.fastq.gz"),raw.get(f"{run}_2.fastq.gz")
        if not os.path.exists(rp) or i1 is None:
            f.write(f"{run}\t{alias}\t{mode}\t{i1 or ''}\t\t\t\tMISSING\n"); bad+=1; continue
        rem=int(open(rp).read()); o1=out.get(f"{run}_1.fastq.gz"); o2=out.get(f"{run}_2.fastq.gz")
        ok = (i1==i2 and o1==o2==i1-rem)
        if not ok: bad+=1
        f.write(f"{run}\t{alias}\t{mode}\t{i1}\t{rem}\t{100*rem/i1:.2f}\t{o1}\t{'OK' if ok else 'MISMATCH'}\n")
print("samples with problems:",bad)
sys.exit(1 if bad else 0)
PY
