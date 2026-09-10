# Hands-On Modules: Metabarcoding, Metagenomics & Pangenomics
### 27221 Microbiome Engineering — DTU, MSc Bioengineering

**Status:** Modules 1 and 2 signed off; notebooks built. Module 3 drafted and
proposed, pending time-budget approval from the course coordinator.

---

## Overview

| | Module 1 — Metabarcoding | Module 2 — Genome-resolved metagenomics | Module 3 — Pangenomics *(proposed)* |
|---|---|---|---|
| **Duration** | 2 × 4h sessions | 2 × 4h sessions | 1 × 4h session |
| **System** | Wastewater treatment plant (WWTP) activated sludge | Plant-root synthetic community (SynCom) | Same SynCom isolates as Module 2 (+ optionally its recovered MAGs) |
| **Toolchain** | R / DADA2 / phyloseq | anvi'o (+ megahit, fastp for small-scale assembly) | anvi'o (pangenomics workflow) |
| **Environment** | JupyterLab + IRkernel, via a shared conda env, accessed through a browser tab inside ThinLinc on the DTU HPC (Gbar/DCC) shared filesystem | Terminal + anvi'o interactive (browser) via ThinLinc on the same shared filesystem | Same as Module 2 — no new environment |
| **Student background assumed** | No prior coding proficiency | No prior coding proficiency | Completed Module 2 |
| **Deliverable** | One figure + short interpretation (per session pair) | One figure + short interpretation (per session pair) | One figure + short interpretation |

Student baseline: molecular microbiology + microbial genetics/physiology (per
course prerequisites), **no assumed programming experience**. All notebooks
are heavily scaffolded — pre-written code with reasoning-required blanks
(`>> YOUR DECISION <<`), not blank cells.

### Environment decision (resolved)

DTU HPC (Gbar/DCC) has no Singularity/Apptainer module and no RStudio
Server/Desktop available — confirmed via `module avail` on the actual
cluster. RStudio Desktop cannot be installed by students themselves either
(no root access, and compute nodes may lack outbound internet). The working
setup for **both** R- and anvi'o-based modules is therefore:

- A **single, shared, pre-built conda environment** (not installed per
  student) staged on the shared filesystem, containing R + DADA2 + phyloseq +
  vegan + tidyverse + IRkernel + JupyterLab for Module 1, and a separate
  anvi'o environment (anvi'o's own recommended conda-based install, plus
  megahit and fastp) for Modules 2 and 3.
- Students activate the relevant environment and, for Module 1, open
  JupyterLab in a browser tab inside their ThinLinc session; for Modules 2/3,
  work from a terminal and view anvi'o's interactive interface in a browser
  tab the same way.
- *(Still open: whether to request Apptainer support from DTU HPC as a
  longer-term improvement for a nicer RStudio experience in future runs of
  the course — not blocking for this iteration.)*

### Notebook execution model

- **Module 1**: `.ipynb` notebooks, run live in JupyterLab (R kernel) via
  ThinLinc — students execute cells directly, not copy-paste. (Originally
  scoped for RStudio/`.Rmd`; converted to Jupyter/IRkernel once RStudio was
  confirmed unavailable on the HPC — see Environment decision above.)
- **Module 2 & 3**: notebooks are structured reference documents (`.ipynb`
  with bash cells, plus matching `.md`) — students copy commands into a
  terminal session, since anvi'o is CLI-driven. The anvi'o interactive
  interface is viewed through a browser tab inside the ThinLinc session.

### Files produced so far

| File | Module | Format |
|---|---|---|
| `M1_S1_reads_to_ASVs.ipynb` (+ `.Rmd`) | 1, Session 1 | Live-executable |
| `M1_S2_ASVs_to_biology.ipynb` (+ `.Rmd`) | 1, Session 2 | Live-executable |
| `M2_S1_reads_to_contigsDB.ipynb` (+ `.md`) | 2, Session 1 | Copy-paste reference |
| `M2_S2_profile_to_MAGs.ipynb` (+ `.md`) | 2, Session 2 | Copy-paste reference |
| `M3_S1_pangenomics.ipynb` (+ `.md`) | 3 *(proposed)* | Copy-paste reference |

---

## Module 1 — Metabarcoding: Wastewater Treatment Plant Microbiomes

### Dataset
- **Source**: raw 16S rRNA V4 amplicon data from the **Global Water Microbiome Consortium** survey (NCBI SRA, BioProject **PRJNA509305**), the same raw dataset underlying the MiDAS 4 database (Dueholm et al. 2022, *Nat Commun*).
- **Curated subset**: instructor pre-selects ~12–20 samples spanning a few countries and/or process configurations (e.g. conventional activated sludge vs. other configurations), staged on the shared HPC filesystem ahead of class, together with a metadata table (country, plant type/process, and any available performance parameters).
- **Metadata**: sample metadata is in Supplementary Table 1 of Wu et al. 2019 (*Nat Microbiol*), freely downloadable. Real process/environmental variables confirmed available include Continent, Influent BOD, Sludge Age, mean annual air temperature, MLSS, DO, F/M ratio, and nutrient removal efficiencies — good candidate "engineering variables" for Session 2. *(You're finalizing the exact sample subset and column names directly from the file.)*
- **Reference taxonomy**: **MiDAS 4** (attribution: Aalborg University / Center for Microbial Communities), used at the DADA2 taxonomy-assignment step — gives genus/species-level resolution that generic universal databases don't achieve well for WWTP bacteria.
- **Why this dataset**: real, globally-sourced, well-documented, metadata-rich, and directly on-topic for a biotechnology/engineering audience — students can ask "does the community track with how the plant is run?" with real data.

### Learning objectives targeted
- Describe structure/function/dynamics of an engineered-environment microbiome.
- Describe and evaluate sequencing-based microbiome characterization methods.
- Analyze/interpret microbiome data using computational and statistical approaches.

### Session 1 (4h) — Raw reads → ASVs
| Time | Segment |
|---|---|
| 0:00–1:00 | **Theory**: amplicon sequencing logic, primer choice/regions, sequencing error models, why ASVs (DADA2) supersede OTU clustering. |
| 1:00–1:15 | Break |
| 1:15–3:30 | **Hands-on** (JupyterLab, live execution): import demultiplexed FASTQs → quality profiling (`plotQualityProfile`) → filter/trim (guided reasoning on truncation length from quality plots) → error learning (`learnErrors`) → denoise (`dada`) → merge pairs → build sequence table → chimera removal → final ASV table. Every cell pre-written; blanks require a parameter choice justified from the preceding plot/output, not free coding. |
| 3:30–4:00 | **Wrap-up**: what does a row/column of an ASV table represent; sanity checks (read-tracking table through the pipeline). |

### Session 2 (4h) — ASVs → biology
| Time | Segment |
|---|---|
| 0:00–0:45 | **Theory**: taxonomy assignment logic, alpha vs. beta diversity, ordination, what a "differential abundance" claim actually means (and its pitfalls). |
| 0:45–1:00 | Break |
| 1:00–3:00 | **Hands-on**: assign taxonomy against MiDAS 4 → build `phyloseq` object with sample metadata → alpha diversity comparison across an engineering variable (e.g. plant/process type, sludge age, F/M ratio) → ordination (e.g. NMDS/PCoA on Bray-Curtis) → one simple differential abundance comparison. |
| 3:00–3:45 | Students produce their deliverable figure + draft interpretation. |
| 3:45–4:00 | Discussion / share-out of a few figures. |

### Deliverable
One figure (diversity plot, ordination, or differential abundance plot) + a short paragraph: what shifts with the chosen engineering variable, and what that would mean for someone operating the plant.

---

## Module 2 — Genome-Resolved Metagenomics: Plant-Root Synthetic Community

### Dataset
- **Primary**: Selten et al. — shotgun metagenomics of a complex, defined plant-root SynCom (NCBI SRA BioProject **PRJNA1131994**), with known constituent isolate genomes serving as **ground truth** for bin recovery/validation — pedagogically valuable since students can check "did we get back what was put in?" *(preprint, as of Sep 2026 — bioRxiv [10.1101/2024.08.22.609090](https://www.biorxiv.org/content/10.1101/2024.08.22.609090), not yet peer-reviewed; the dataset itself is public and stable regardless.)*
- **Bonus/bridge exercise**: the same samples were also sequenced by 16S amplicon (BioProject **PRJNA1191388**) — optional short comparison at the end of Session 2, contrasting what Module 1's amplicon approach would show vs. genome-resolved metagenomics on the *same* system.
- **Framing story**: open Session 1's theory block with Carrión et al. 2019 (*Science*) — the 2-strain sugarbeet SynCom where metagenomics pinned a biocontrol mechanism to specific gene clusters — as a compact, well-known illustration of *why* genome-resolved metagenomics matters for engineered communities, before moving to the more complex dataset used for hands-on work.

### Learning objectives targeted
- Describe/evaluate microbiome characterization methods (sequencing-based, model systems).
- Explain synthetic biology / systems biology advances relevant to microbiome engineering.
- Design strategies to manipulate microbiomes (ties to "designing microbial communities").
- Analyze/interpret microbiome data computationally.

### Session 1 (4h) — Small-scale, live: reads → contigs database
| Time | Segment |
|---|---|
| 0:00–1:00 | **Theory**: shotgun vs. amplicon, assembly basics, what a contigs database is, the Carrión et al. framing story, why this matters for engineered/designed communities. |
| 1:00–1:15 | Break |
| 1:15–3:30 | **Hands-on** (terminal via ThinLinc, on a small pre-subsampled read set — small enough to finish live): QC (`fastp`) → quick assembly (`megahit`) → `anvi-gen-contigs-database` → `anvi-run-hmms` → gene calling / basic functional annotation. |
| 3:30–4:00 | **Wrap-up**: first look at a contigs database and the `anvi-interactive` interface (opened in browser via ThinLinc). |

### Session 2 (4h) — Precomputed profile → genome-resolved exploration
| Time | Segment |
|---|---|
| 0:00–0:45 | **Theory**: read recruitment/mapping concept, binning logic, MAG quality metrics (completion/redundancy). |
| 0:45–1:00 | Break |
| 1:00–3:00 | **Hands-on**: explore a pre-built, full-scale merged profile (too large to generate live) in `anvi-interactive` → binning → pull out a bin of interest → compare against known isolate genomes (ground truth, via ANI) → basic functional/taxonomic characterization of the recovered bin. |
| 3:00–3:30 | *(Optional bridge)*: brief comparison against the matched 16S amplicon data (PRJNA1191388) for the same samples — what would Module 1's approach have told us, vs. what genome-resolved metagenomics reveals? |
| 3:30–4:00 | Students produce their deliverable figure + draft interpretation; share-out. |

### Deliverable
One figure (e.g. annotated bin overview, or an `anvi-interactive`/binning screenshot) + a short paragraph: which community members were recovered, one functional or taxonomic claim about them, and what engineering decision it could inform (e.g. relevant to the SynCom's design goal).

---

## Module 3 — Pangenomics *(proposed)*

**Status:** proposed, pending time-budget approval from the course coordinator.
Designed as a direct continuation of Module 2 — same dataset, same tool, no
new environment to build.

### Dataset
Same SynCom isolate genomes already staged for Module 2's ground-truth
comparison step, reused directly as the pangenome input — no new data-staging
effort. Optionally, the MAGs students recovered in Module 2 can be added as
"internal genomes" alongside the isolates, to ask whether a recovered MAG
sits where its matched isolate does in gene-content space.

### Learning objectives targeted
- Design strategies to manipulate/design microbial communities (redundancy vs.
  complementarity is a direct design question for a SynCom).
- Explain synthetic biology / systems biology advances relevant to
  microbiome engineering.
- Analyze/interpret microbiome data computationally.

### Session (4h) — Core, accessory, and the redundancy question
| Time | Segment |
|---|---|
| 0:00–0:45 | **Theory**: what a pangenome is (core vs. accessory vs. singleton gene clusters), how anvi'o builds one (protein clustering via MCL on translated ORFs), and why this matters for a *designed* community — functional redundancy vs. complementarity across strains is a direct SynCom design question ("do we need all of these members?"). |
| 0:45–1:00 | Break |
| 1:00–3:00 | **Hands-on**: build a genomes storage from the isolate genomes → run the pangenome analysis → explore interactively (`anvi-display-pan`) → extract core/accessory/singleton counts per genome → functional enrichment of the accessory genome (reusing the COG annotation from Module 2). |
| 3:00–3:45 | Students produce their deliverable figure + draft interpretation. |
| 3:45–4:00 | Discussion / share-out. |

### Deliverable
One figure (e.g. an annotated `anvi-display-pan` screenshot, or a bar plot of
core/accessory/singleton gene cluster counts per isolate) + a short paragraph:
which isolate(s) look most functionally redundant with others in the SynCom,
which look most distinct, and what that would mean for a decision to simplify
the community (drop a strain) or not.

---

## Logistics checklist (pre-class, instructor-side)

- [ ] Confirm ThinLinc access + shared filesystem paths for all enrolled students.
- [ ] Build and stage a **shared, read-only conda environment** for Module 1 (R + DADA2 + phyloseq + vegan + tidyverse + IRkernel + JupyterLab) — all students activate the same env rather than building their own.
- [ ] Build and stage a **shared anvi'o conda environment** for Modules 2/3 (anvi'o's own recommended install, plus `fastp` and `megahit`).
- [ ] Module 1: finalize the ~12–20 sample WWTP subset from PRJNA509305 + the confirmed metadata columns from Supplementary Table 1; stage MiDAS 4 reference taxonomy files.
- [ ] Module 2: select and stage a small subsampled read set from PRJNA1131994 for live Session 1 use; pre-build the full merged profile database for Session 2 (too heavy to build live); stage known isolate reference genomes and the `external-genomes.txt` manifest for ground-truth comparison; stage matched amplicon data (PRJNA1191388) if using the bonus bridge exercise.
- [ ] Module 3 *(if approved)*: no new staging needed beyond what Module 2 already requires — reuses the same isolate genomes and `external-genomes.txt`.
- [ ] Confirm anvi'o's local web server for `anvi-interactive`/`anvi-display-pan` and JupyterLab are both reachable from within a ThinLinc browser session (port/proxy considerations on the HPC).
- [ ] Sanity-check exact anvi'o CLI flags used in the Module 2/3 notebooks against whatever anvi'o version is actually installed (command syntax has changed across versions).
- [ ] Decide on final deliverable submission mechanism (shared filesystem folder, LMS upload, etc.).

## Next steps

1. Course coordinator decision on Module 3's time allocation.
2. You finalize and stage real data for Module 1 (WWTP subset + metadata) and Modules 2/3 (SynCom subsamples, precomputed profile, isolate genomes) on the shared filesystem.
3. Build the two shared conda environments (R/Jupyter for Module 1; anvi'o for Modules 2/3) — happy to draft the `environment.yml`/install steps for either whenever useful.
4. Once real paths exist, update the placeholder config values in each notebook (`raw_reads_dir`, `CONTIGS_DB`, `PROFILE_DB`, etc.) to match.
5. Optional, longer-term: ask DTU HPC whether Apptainer/Singularity support could be added, to enable a smoother RStudio Server setup in future runs.
