# A dynamical‑systems analysis of the *toc1‑1* circadian transcriptome

Code and data accompanying:

> **A dynamical‑systems analysis of the *toc1‑1* circadian transcriptome identifies
> transcripts coupled to the TOC1‑independent cytosolic calcium rhythm in *Arabidopsis*.**
> Hearn, Robertson, Briggs and Webb; [Journal]; [Year]. [DOI]

This repository reproduces every figure, table and supplementary item in the paper
from the raw input data.

## Background

The *toc1‑1* mutation shortens the period of the *Arabidopsis* circadian clock, yet
the cytosolic free‑calcium ([Ca²⁺]cyt) rhythm keeps its period. We re‑analyse a
whole‑transcriptome ATH1 microarray time‑course recorded alongside a concurrent
[Ca²⁺]cyt recording (wild type C24 and *toc1‑1*, sampled every 4 h from 49–93 h in
constant light), combining conventional rhythm analysis with the **Vinnicombe
ν‑gap**, a bounded distance between dynamical models. Two ν‑gap analyses are used:
a transcript–transcript interaction‑network comparison, and a transcript↔[Ca²⁺]cyt
input/output comparison (in both directions) that identifies transcripts whose
dynamic coupling to [Ca²⁺]cyt is preserved in *toc1‑1*.

The microarray data are deposited at ArrayExpress under accession **E‑GEOD‑19271**.

## Repository layout

The repository contains three self‑contained subdirectories, each with its own
`README.md`, input `data/`, and outputs.

| Subdirectory | Purpose | Produces |
|---|---|---|
| [`toc1_calcium_workflow_reproducibility`](toc1_calcium_workflow_reproducibility) | The main analysis pipeline: data loading/alignment, cosine rhythm analysis, the ν‑gap interaction network, and the bidirectional transcript↔[Ca²⁺]cyt coupling analysis. | Main‑text **Figures 1–4**, **Supplementary Figure S1**, and **Supplementary Tables S1–S5** |
| [`calcium_graph_visualisation`](calcium_graph_visualisation) | Standalone hub‑and‑spoke network rendering of the direction‑robust calcium‑coupled transcripts. | **Supplementary Figure S2** (`.png` and editable `.svg`) |
| [`cosopt_reproducibility`](cosopt_reproducibility) | Cross‑check substituting the original **COSOPT** period analysis into the pipeline, and the concordance with the cosine method. | **Supplementary Figure S3** and **Supplementary Table S6** |

### Where each paper item is generated

| Paper item | Subdirectory |
|---|---|
| Figures 1–4 (clock/calcium, genome‑wide effects, ν‑gap analyses, non‑shortening genes) | `toc1_calcium_workflow_reproducibility` |
| Figure S1 (forward vs reverse gene↔[Ca²⁺]cyt ν‑gap) | `toc1_calcium_workflow_reproducibility` |
| Figure S2 (calcium‑coupling network) | `calcium_graph_visualisation` |
| Figure S3 (cosine vs COSOPT concordance) | `cosopt_reproducibility` |
| Tables S1–S2 (non‑shortening set; direction‑robust candidates) | `toc1_calcium_workflow_reproducibility` |
| Tables S3–S4 (full forward / reverse gene↔[Ca²⁺]cyt results) | `toc1_calcium_workflow_reproducibility` |
| Table S5 (full per‑transcript cosine period analysis) | `toc1_calcium_workflow_reproducibility` |
| Table S6 (COSOPT period analysis and concordance) | `cosopt_reproducibility` |

## Requirements

Python ≥ 3.9 with:

```bash
pip install numpy scipy pandas matplotlib networkx openpyxl
```

The main workflow additionally uses the **`nugap`** package (the ν‑gap metric and
model fitting), which is bundled inside
`toc1_calcium_workflow_reproducibility/nugap_package/` and can be installed with
`pip install ./toc1_calcium_workflow_reproducibility/nugap_package` (the analysis
script also finds it automatically if it sits alongside).

## Reproducing the results

Each subdirectory is independent; see its `README.md` for full detail. In brief:

```bash
# 1) Main workflow — Figures 1–4, Fig S1, Tables S1–S5
cd toc1_calcium_workflow_reproducibility
pip install ./nugap_package           # or: export PYTHONPATH=$PWD/nugap_package/src
python toc1_nugap_analysis.py data/100725E1c24.xlsx data/100725E2c24.xlsx \
       data/input_calcium_data.xlsx results

# 2) Calcium‑coupling network — Figure S2
cd ../calcium_graph_visualisation
python network_figure.py

# 3) COSOPT cross‑check — Figure S3, Table S6
cd ../cosopt_reproducibility
python cosopt_period_analysis.py data/100725E1c24_cosopt.xlsx \
       data/100725E2c24_cosopt.xlsx data/cosine_period_map.csv output
```

The main workflow runs in a few minutes (the pairwise ν‑gap network is the slow
step). The other two are near‑instant.

## Notes on method

- **Rhythm parameters** (period, amplitude, phase) are estimated with a
  variable‑period cosine fit and are the robust readout for sustained oscillations
  sampled at 12 time‑points.
- **The ν‑gap** (Vinnicombe metric, range [0, 1]) is used for *dynamic‑distance*
  questions. Input and output signals are standardised to unit variance before
  fitting so the metric is not saturated by the large scale difference between
  luminescence and log₂ expression. A low‑order linear model cannot represent a
  sustained limit cycle, so autonomous (output‑only) ν‑gap fits to single
  transcripts are deliberately not used for inference.
- **Directionality.** The transcript↔[Ca²⁺]cyt analysis is run in both directions
  (transcript→[Ca²⁺]cyt and [Ca²⁺]cyt→transcript); the two are only weakly
  correlated, so the intersection (preserved in both directions) is treated as the
  higher‑confidence candidate set. These models are observational and do not
  establish causality.
- **Period method.** The COSOPT cross‑check (subdirectory 3) shows the period and
  rhythmicity results, and the period‑invariant gene assignments, are concordant
  between the cosine analysis and COSOPT (per‑gene period Pearson r ≥ 0.93).

## Data availability

- Microarray time‑course: ArrayExpress **E‑GEOD‑19271**.
- The exact input files used here (microarray, COSOPT output, and the concurrent
  [Ca²⁺]cyt recording) are included in the relevant `data/` subfolders so that the
  analyses run without external downloads.

## License

MIT

## Contact

Tim Hearn tjh70@cam.ac.uk
