# Reproducibility bundle — *toc1-1* circadian transcriptome / \[Ca²⁺]cyt ν-gap analysis

This bundle reproduces the analytical figures and tables in the manuscript

> \*A dynamical systems analysis of the toc1 1 circadian transcriptome identifies transcripts coupled to the TOC1 independent cytosolic calcium rhythm in Arabidopsis\*

from the raw input files, with a single script.

## Contents

```
toc1\_nugap\_analysis.py     the complete analysis (load → figures + tables)
data/
  100725E1c24.xlsx         ATH1 microarray, replicate experiment 1 (C24 WT + toc1-1)
  100725E2c24.xlsx         ATH1 microarray, replicate experiment 2
  input\_calcium\_data.xlsx  concurrent \[Ca2+]cyt timecourse (aequorin), WT + toc1-1
nugap/                     the ν-gap Python package (Vinnicombe metric + fitting)
results/                   reference output (the figures/tables the script produces)
```

## Input data

Affymetrix ATH1 arrays, two biological replicates, wild type (C24) and *toc1-1*,
sampled every 4 h from 49–93 h in constant light (12 time points, circadian cycles
3–4). ArrayExpress accession **E-GEOD-19271**. Each microarray file (sheet `All`)
holds the probe id, locus, annotation, then a raw log₂ block of `WT{rep}\_49..93`
and `T{rep}\_49..93` columns. The two replicate files are **not** in a common probe
order and are merged on probe id by the script. The calcium file (sheet `Sheet1`)
has columns `time`, `wt`, `toc1-1` on the same timebase.

## Requirements

* Python ≥ 3.9
* `numpy`, `scipy`, `pandas`, `matplotlib`
* our `nugap` package 

```bash
pip install numpy scipy pandas matplotlib
pip install nugap
```

The nugap package is also included in this bundle, as it has not been formerly published in an academic journal, although is available on PyPi.

## Run

```bash
python toc1\_nugap\_analysis.py data/100725E1c24.xlsx data/100725E2c24.xlsx \\
       data/input\_calcium\_data.xlsx results
```

Runtime is a few minutes (the pairwise ν-gap network is the slow step).

## Outputs

Figures (`results/`):

* `fig1\_clock\_calcium.png` — clock transcripts short-period in *toc1-1*; \[Ca²⁺]cyt unchanged
* `fig2\_genomewide.png` — period shift, rhythmicity loss, amplitude of lost vs retained genes
* `fig3\_nugap.png` — ν-gap network (global change, no significant edge) and gene→\[Ca²⁺]cyt coupling
* `fig4\_invariant\_genes.png` — all period-invariant transcripts (WT vs *toc1-1*)

Figures also include `figS1\_bidirectional.png` (forward vs reverse gene↔\[Ca²⁺]cyt ν-gap).

Tables (`results/`):

* `table\_period\_map.csv` — per-gene period, amplitude, rhythmicity calls, period\_class
* `table\_nonshortening\_genes.csv` — the 21 non-shortening transcripts (19 invariant + 2 lengthened)
* `table\_gene\_calcium.csv` — forward and reverse gene↔\[Ca²⁺]cyt ν-gap, fit quality, floors, robust\_both flag
* `table\_network\_edges.csv` — pairwise ν-gap network comparison
* `suppl\_table\_S1\_nonshortening.csv` — Table S1 (21 non-shortening transcripts, both directions)
* `suppl\_table\_S2\_robust\_both.csv` — Table S2 (14 transcripts preserved in both directions)
* `suppl\_table\_S3\_gene\_to\_calcium.csv` — Table S3 (all 274 reliable forward transcript→\[Ca²⁺]cyt)
* `suppl\_table\_S4\_calcium\_to\_gene.csv` — Table S4 (all 199 reliable reverse \[Ca²⁺]cyt→transcript)

Formatted supplementary files (in the paper outputs): `supplementary\_information.docx`
(Tables S1–S2 and Figure S1), `suppl\_table\_S3\_gene\_to\_calcium.xlsx`, `suppl\_table\_S4\_calcium\_to\_gene.xlsx`.

## Notes on method

* Rhythm parameters (period, amplitude, phase) use a variable-period cosine fit; they
are the robust readout for sustained oscillations sampled at 12 points.
* The ν-gap (Vinnicombe metric, \[0,1]) is used for *dynamic-distance* questions:
(i) how transcript–transcript relationships change between genotypes (network), and
(ii) whether a transcript→\[Ca²⁺]cyt input–output relationship is preserved in *toc1-1*.
Signals are standardised to unit variance before fitting so the metric is not
saturated by the large scale difference between luminescence and log₂ expression.
* A low-order linear model cannot represent a sustained limit cycle, so autonomous
(output-only) ν-gap fits to single transcripts are *not* used for inference.

