# Module 2 — Session 1: From raw reads to a contigs database
### 27221 Microbiome Engineering — Genome-resolved metagenomics of a plant-root synthetic community (SynCom)

**How to use this document:** this is a reference, not a live notebook — copy each
command into your terminal session (opened via ThinLinc), run it, and read the
output before moving on. Where a choice is marked `>> YOUR DECISION <<`, note
your reasoning; there's no single correct answer.

---

## 0. Context

A **synthetic community (SynCom)** is a defined, reconstituted set of cultured
microbial isolates inoculated onto a host — here, plant roots. Because the
identity of every member is known in advance, SynCom metagenomes are an
unusually good teaching system for genome-resolved metagenomics: you can check
whether what you *recover* computationally matches what was actually put in.

As a real-world illustration of why this matters for engineering: Carrión et
al. (2019, *Science*) showed that a two-strain SynCom (a *Flavobacterium* and
a *Chitinophaga* isolate) could reproduce the disease-suppressive effect of
the full sugar beet root microbiome against a fungal pathogen — and that
metagenome data was what pinned this protective function down to specific
genes in these two species. Genome-resolved metagenomics is what let them
move from "the community protects the plant" to "these two organisms, and
these specific genes, protect the plant."

Today's dataset is more complex than that: a multi-member root SynCom, shotgun
sequenced (Selten et al.), with the isolate genomes used to build the SynCom
available as ground truth.

---

## 1. Setup

```bash
# Activate the shared conda environment (ask your instructor for the exact name if different)
conda activate anvio-8

# --- SESSION CONFIG: adjust these paths if your instructor gives you different ones ---
export RAW_READS_DIR=/path/to/shared/syncom_metagenomics/subsampled_reads   # small, live-processing subset
export WORKDIR=~/module2_output
mkdir -p $WORKDIR
cd $WORKDIR
```

```bash
# Sanity check: what do we have?
ls -lh $RAW_READS_DIR
```

**Question:** How many samples/read files do you have, and roughly how large
are they? (This subset has been deliberately subsampled down from the full
dataset so that assembly finishes within the session — the full dataset is
what you'll explore in Session 2, already pre-assembled.)

---

## 2. Quality control

```bash
mkdir -p $WORKDIR/qc
for R1 in $RAW_READS_DIR/*_R1.fastq.gz; do
  SAMPLE=$(basename $R1 _R1.fastq.gz)
  R2=$RAW_READS_DIR/${SAMPLE}_R2.fastq.gz
  fastp \
    -i $R1 -I $R2 \
    -o $WORKDIR/qc/${SAMPLE}_R1.trimmed.fastq.gz \
    -O $WORKDIR/qc/${SAMPLE}_R2.trimmed.fastq.gz \
    --json $WORKDIR/qc/${SAMPLE}.fastp.json \
    --html $WORKDIR/qc/${SAMPLE}.fastp.html
done
```

**Question:** Open one of the `.fastp.html` reports (via the ThinLinc file
browser or a local browser pointed at the shared filesystem). What fraction of
reads/bases were removed? Does anything look unusual (adapter content,
duplication level)?

---

## 3. Assembly

```bash
# Example for a single sample — repeat, or loop, for each sample in your subset
megahit \
  -1 $WORKDIR/qc/SAMPLE_R1.trimmed.fastq.gz \
  -2 $WORKDIR/qc/SAMPLE_R2.trimmed.fastq.gz \
  -o $WORKDIR/assembly/SAMPLE \
  --min-contig-len 1000
```

**`>> YOUR DECISION <<`** — `--min-contig-len` sets the shortest contig
megahit will report. We've suggested 1000 bp as a reasonable default for
genome-resolved metagenomics (very short contigs rarely bin well and mostly
add noise). Would you expect a much smaller value (e.g. 300) to help or hurt
downstream binning, and why?

```bash
# Basic assembly stats
seqkit stats $WORKDIR/assembly/SAMPLE/final.contigs.fa
```

**Question:** How many contigs did you get, what's the N50, and what's the
longest contig? Compare notes with a classmate who assembled a different
sample — is assembly quality consistent across samples?

---

## 4. Build a contigs database

The **contigs database** is anvi'o's central data structure: it stores contig
sequences plus every piece of downstream annotation (genes, functions,
taxonomy, HMM hits) in one place.

```bash
# anvi'o requires simple contig deflines — reformat first
anvi-script-reformat-fasta \
  $WORKDIR/assembly/SAMPLE/final.contigs.fa \
  -o $WORKDIR/assembly/SAMPLE/contigs-fixed.fa \
  -l 1000 --simplify-names --report-file $WORKDIR/assembly/SAMPLE/rename-report.txt

anvi-gen-contigs-database \
  -f $WORKDIR/assembly/SAMPLE/contigs-fixed.fa \
  -o $WORKDIR/SAMPLE-CONTIGS.db \
  -n "SynCom sample SAMPLE"
```

---

## 5. Annotate: single-copy core genes, and gene calls

```bash
# HMM search for bacterial/archaeal single-copy core genes — the backbone of
# later completion/redundancy estimates for any bins you recover
anvi-run-hmms -c $WORKDIR/SAMPLE-CONTIGS.db --num-threads 4
```

```bash
# Optional, if time allows: functional annotation
anvi-run-ncbi-cogs -c $WORKDIR/SAMPLE-CONTIGS.db --num-threads 4
```

```bash
# Quick summary of what's in the contigs database so far
anvi-display-contigs-stats $WORKDIR/SAMPLE-CONTIGS.db
```

`anvi-display-contigs-stats` starts a small local web server — open the URL it
prints in a browser tab inside your ThinLinc session.

**Question:** How many single-copy core gene hits did anvi'o find, and how
does that compare to the number you'd expect from a handful of bacterial
genomes (a single bacterial genome typically carries close to one copy of
each of ~70–140 commonly used marker genes)? What would a much lower count
suggest about your assembly?

---

## Wrap-up discussion

- You've now built one small contigs database from one (or a few) subsampled
  sample(s). In Session 2 you'll work with a full-scale, pre-built, merged
  profile spanning the whole SynCom experiment — too large to assemble live.
- What was the practical bottleneck today: read volume, assembly time, or
  something else? That's exactly why the full dataset is precomputed for next
  session.
- Keep your `SAMPLE-CONTIGS.db` — you'll compare what a small-scale, single
  sample can tell you against the full profile in Session 2.
