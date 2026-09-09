# Module 2 — Session 2: From a merged profile to recovered genomes
### 27221 Microbiome Engineering — Genome-resolved metagenomics of a plant-root synthetic community (SynCom)

**How to use this document:** reference material — copy commands into your
terminal (via ThinLinc), and open anvi'o's interactive interface in a browser
tab inside the same session. `>> YOUR DECISION <<` marks a judgment call.

---

## 0. What's different today

Last session you built one contigs database from a small, subsampled read
set — enough to see every step, but too little data to responsibly bin
genomes. Today you'll work with a **precomputed, full-scale merged profile**
across the whole SynCom experiment (built the same way you did last session,
just at real scale, ahead of time by your instructor) using the actual reads
recruited from all samples against the assembled contigs.

The point of today: recover **metagenome-assembled genomes (MAGs)** from the
community, and — because this is a *defined* SynCom — check them against the
**known isolate genomes** that went into building it.

---

## 1. Setup

```bash
conda activate anvio-8

# --- SESSION CONFIG: adjust these paths if your instructor gives you different ones ---
export CONTIGS_DB=/path/to/shared/syncom_metagenomics/precomputed/CONTIGS.db
export PROFILE_DB=/path/to/shared/syncom_metagenomics/precomputed/PROFILE.db      # merged, full-scale
export ISOLATE_GENOMES_DIR=/path/to/shared/syncom_metagenomics/isolate_genomes    # ground truth
export WORKDIR=~/module2_output
cd $WORKDIR
```

```bash
# Sanity check
anvi-db-info $CONTIGS_DB
anvi-db-info $PROFILE_DB
```

**Question:** How many contigs and how many samples does the merged profile
cover? Compare this to your single-sample database from Session 1.

---

## 2. Read recruitment: the concept

Before binning, it's worth understanding what actually produced this merged
profile: reads from *every* sample were mapped ("recruited") back onto the
assembled contigs, and per-contig, per-sample coverage was recorded. A
contig's **coverage pattern across samples** — not just its sequence
composition — is one of the strongest signals for figuring out which contigs
belong to the same genome (organisms present in similar relative amounts
across samples tend to have correlated coverage).

---

## 3. Explore the profile interactively

```bash
anvi-interactive -c $CONTIGS_DB -p $PROFILE_DB --server-only -P 8080
```

Open the printed URL in a browser tab inside your ThinLinc session. You'll see
a circular display: each "spoke" is a contig (or cluster of related contigs),
arranged by sequence composition and coverage similarity, with coverage
across your samples shown as concentric rings.

**Question:** Do you see visually distinct clusters of contigs? Roughly how
many major clusters — and does that number look consistent with the number of
isolates you'd expect to have been added to this SynCom?

---

## 4. Binning

```bash
# If you want to try automatic binning as a starting point (e.g. CONCOCT),
# your instructor may have this pre-run and stored in the profile database
# already — check with anvi-show-collections-and-bins first:
anvi-show-collections-and-bins -p $PROFILE_DB
```

**`>> YOUR DECISION <<`** — if a default collection already exists, start from
it and refine; otherwise, use the interactive interface's manual binning tools
(selecting clusters of contigs by dragging on the dendrogram) to draw your own
bins. Aim for a small number of bins that correspond to visually coherent
clusters, rather than many fragmented ones.

```bash
# Once you've defined bins (named "SynCom_bins" here, adjust to your instructor's convention):
anvi-summarize -c $CONTIGS_DB -p $PROFILE_DB -C SynCom_bins -o $WORKDIR/summary
```

`anvi-summarize` produces a table of completion/redundancy and basic stats for
every bin — open `$WORKDIR/summary/bins_summary.txt`.

**Question:** For your best bin, what's its completion (% of expected
single-copy core genes found) and redundancy (% found in more than one
copy — a sign of contamination/mixing)? A commonly used "good enough for a
MAG" threshold is roughly >70% completion and <10% redundancy — how does your
bin compare?

---

## 5. Ground truth: does your bin match a known isolate?

```bash
# Compare a recovered bin's genome sequence against the known isolate genomes
anvi-compute-genome-similarity \
  --external-genomes external-genomes.txt \
  --internal-genomes internal-genomes.txt \
  --program pyANI \
  -o $WORKDIR/ani_comparison
```

(Your instructor will provide the `external-genomes.txt` / `internal-genomes.txt`
files pointing at the isolate genomes in `$ISOLATE_GENOMES_DIR` and your
recovered bins — ask if these aren't already staged.)

**Question:** Does your recovered bin match a known isolate at high average
nucleotide identity (ANI; >95% is a common species-level threshold)? If a bin
*doesn't* clearly match any single isolate, what could explain that — multiple
closely related isolates it didn't distinguish between, or a genuinely
under-assembled/composite bin?

---

## 6. Functional/taxonomic characterization

```bash
anvi-estimate-scg-taxonomy -c $CONTIGS_DB -p $PROFILE_DB -C SynCom_bins --metagenome-mode
```

**Question:** What taxonomy does anvi'o assign your bin of interest? Does it
agree with the ANI-based ground-truth match from Step 5?

---

## 7. *(Optional bridge)* Amplicon vs. shotgun on the same samples

This SynCom was also sequenced by 16S amplicon (the same approach you used in
Module 1). Discuss with your group:

- Which isolates would 16S amplicon sequencing likely have distinguished
  clearly, and which would it likely have blurred together (think about how
  closely related different SynCom members are)?
- What did genome-resolved metagenomics give you today that amplicon
  sequencing alone could not have?

---

## 8. Produce your deliverable

**Choose ONE figure**: an annotated view of your bin (e.g. a screenshot of the
`anvi-interactive` display with your bin highlighted, or a simple bar chart of
completion/redundancy across your bins) that best supports a claim about what
was recovered from this community.

**Write your interpretation (2–5 sentences):** which community member(s) did
you recover, what evidence supports the identification (completion,
redundancy, ANI match, taxonomy), and what engineering decision about this
SynCom could this inform (e.g. confirming a strain colonized/persisted,
flagging one that didn't assemble well and may need deeper sequencing)?

---

## Wrap-up discussion

- What judgment calls did you make today (binning boundaries, thresholds),
  and how might they affect your conclusions?
- If this had been a natural (non-synthetic) root microbiome instead of a
  SynCom, what would you have lost by not having ground-truth isolate genomes
  to check against?
