# Binner comparison: CONCOCT vs MetaBAT2 vs SemiBin2 vs DAS_Tool consensus

Computed from `anvi-export-collection` (split level) joined with split
lengths from the contigs DB. Regenerate with the same two commands plus the
join script described in `module2_session2_binning.md` (project memory).

| Collection | Bins | Total Mbp | % of 515.2 Mbp assembly | Median bin | Min bin | Bins <50 kb | Bins ≥1 Mbp |
|---|---|---|---|---|---|---|---|
| CONCOCT | 139 | 515.0 | 99.96% | 3.25 Mbp | 1.5 kb | 19 | 93 |
| METABAT2 | 136 | 293.4 | 56.9% | 990 kb | 200.6 kb | 0 | 67 |
| SEMIBIN2 | 517 | 349.6 | 67.9% | 8.7 kb | 2.5 kb | 339 (66%) | 89 |
| **DASTOOL_CONSENSUS** | **35** | 183.2 | 35.6% | **4.95 Mbp** | **1.74 Mbp** | 0 | 35 (100%) |

Per-bin size percentiles (kb):

| Collection | p10 | p25 | p50 | p75 | p90 |
|---|---|---|---|---|---|
| CONCOCT | 13.7 | 233.3 | 3253.5 | 4930.7 | 8792.3 |
| METABAT2 | 258.4 | 370.1 | 995.2 | 3484.6 | 5387.7 |
| SEMIBIN2 | 2.8 | 3.4 | 8.7 | 245.0 | 3195.4 |
| DASTOOL_CONSENSUS | 3253.5 | 4122.7 | 4951.6 | 6101.0 | 8043.4 |

## Interpretation (teaching point)

- **CONCOCT** force-bins almost the entire assembly (515.0 of 515.2 Mbp) —
  it never refuses to call a bin, so some bins likely merge more than one
  organism.
- **MetaBAT2**'s smallest bin is exactly 200.6 kb — not a coincidence, it
  matches its built-in `minClsSize` default (200,000 bp): it discards ~43%
  of the assembly rather than output a low-confidence bin.
- **SemiBin2** is bimodal: 339 of 517 bins (66%) are under 50 kb (fragments)
  alongside 89 genuinely genome-scale (≥1 Mbp) bins.
- **DAS_Tool's consensus** is the cleanest outcome: all 35 final bins are
  ≥1.74 Mbp — 100% genome-scale, zero fragments. It traded total recovered
  sequence (35.6% of the assembly, the least of the four) for precision,
  matching its design (Sieber et al. 2018): keep only what the
  single-copy-marker scoring across all three raw binners actually
  supports.

## Known gap

Completion/redundancy per bin (`anvi-estimate-genome-completeness -C
<collection>`) is not included here — it depended on a scikit-learn
classifier that was broken in the `anvio-9` env at the time (see memory,
fixed 2026-10-02). Worth re-running now that the env is fixed.
