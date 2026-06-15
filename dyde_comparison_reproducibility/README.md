# DyDE method comparison (standalone)

A **self-contained** bundle that reproduces the methodological comparison between
our transcript↔[Ca²⁺]cyt ν‑gap analysis and the original DyDE framework
(Mombaerts et al. 2019, *PLoS Comput Biol* 15(1):e1006674). It quantifies how
much the candidate set depends on two implementation choices, and shows that the
conclusions do not. It is kept separate from the main analysis bundle.

The main analysis uses the **more stringent** choices: the full‑spectrum
Vinnicombe ν‑gap **with** the winding‑number condition, and a reliability gate on
the simulation R² (> 0.5 in both genotypes). DyDE instead (i) measures the ν‑gap
only over the oscillation frequency band and skips the winding‑number test, and
(ii) additionally rejects models with absolute DC gain ≤ 0.1. This bundle
recomputes the gene↔[Ca²⁺]cyt results under the DyDE‑style choices and compares
them to the main analysis.

## Contents

```
nugap_method_sensitivity.py    recompute & compare the two metrics; writes output/
nugap_package/                 the nugap package (v0.2.0) with the band-limited
                               nu-gap and DC-gain-floor options
data/
  100725E1c24.xlsx             ATH1 microarray, replicate 1 (WT C24 + toc1-1)
  100725E2c24.xlsx             ATH1 microarray, replicate 2
  input_calcium_data.xlsx      concurrent [Ca2+]cyt time-course (WT, toc1-1)
  table_gene_calcium.csv       transcripts carried into the calcium analysis,
                               with the main (broadband) forward/reverse nu-gaps
output/
  figS4_nugap_method_sensitivity.png   Supplementary Figure S4
  suppl_table_S7_candidate_comparison.xlsx  Supplementary Table S7 (14 vs 29)
  supplementary_note_1.docx            Supplementary Note 1 (write-up)
  nugap_sensitivity_summary.txt        headline numbers
  nugap_sensitivity_table.csv          per-transcript values, both metrics
```

## Requirements & run

```bash
pip install numpy scipy pandas matplotlib openpyxl
pip install ./nugap_package            # or: export PYTHONPATH=$PWD/nugap_package/src
python nugap_method_sensitivity.py \
       data/100725E1c24.xlsx data/100725E2c24.xlsx \
       data/input_calcium_data.xlsx data/table_gene_calcium.csv output
```

(The script also finds `nugap_package/` automatically if the package is not
installed.)

## What it shows

- Forward (transcript→[Ca²⁺]cyt) and reverse ([Ca²⁺]cyt→transcript) ν‑gaps agree
  closely between the full‑spectrum (main) and band‑limited (DyDE‑style) metrics:
  Pearson r = 0.96 in both directions.
- The between‑genotype > within‑genotype signal is unchanged under the
  band‑limited metric (Wilcoxon p = 2×10⁻¹⁹).
- The 14 direction‑robust candidates from the main analysis are a **strict
  subset** of the 29 found with the band‑limited metric — every main‑text
  candidate is recovered, and the main analysis is the more conservative.
- All 14 main candidates have a forward‑model DC gain ≥ 0.77 (DyDE floor 0.1), so
  the DC‑gain filter removes none of them.

The band‑limited ν‑gap and the DC‑gain floor are exposed as options in the
`nugap` package: `nu_gap(..., band=(2*pi/Pmax, 2*pi/Pmin), check_winding=False)`
and `fit_first_order(..., min_dc_gain=0.1)`.
