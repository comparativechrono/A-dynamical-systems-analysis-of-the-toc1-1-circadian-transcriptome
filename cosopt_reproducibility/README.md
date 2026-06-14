# COSOPT cross-check bundle

A self-contained bundle that re-runs the period/rhythmicity layer of
the study with the original **COSOPT** analysis instead of the variable-period
cosine fit, and quantifies the concordance between the two methods. Its purpose
is to show that the conclusions do not depend on the choice of period-analysis
tool.

COSOPT affects only the period and rhythmicity calls. The ν-gap analyses
(transcript↔[Ca²⁺]cyt coupling, the interaction network) use the raw expression
time-courses and are therefore unchanged by this substitution.

## Inputs

Two replicate COSOPT result files (sheet `All`, merged on probe `Name`; columns
per WT and T = *toc1-1*: `Amp`, `P/A`, `MEL` = mean expression, `Per` = period in
hours, `Phase`, `pMMC` = pMMC-β statistic):

```
data/100725E1c24_cosopt.xlsx
data/100725E2c24_cosopt.xlsx
data/cosine_period_map.csv      cosine-method period map (= results/table_period_map.csv
                                from the main reproducibility bundle), used only for the
                                method-to-method concordance comparison
```

**Rhythmicity** is called as COSOPT pMMC-β < 0.1 in **both** replicates (the
criterion of the original study). **Period** is the mean of the two replicate
COSOPT period estimates. The downstream classification (rhythmic in each genotype
→ period change → period-invariant / shortened / lengthened → "does not shorten")
is identical to the main pipeline.

## Requirements & run

```bash
pip install pandas numpy scipy matplotlib openpyxl
python cosopt_period_analysis.py data/100725E1c24_cosopt.xlsx \
       data/100725E2c24_cosopt.xlsx data/cosine_period_map.csv output
```

## Outputs (`output/`)

- `cosopt_period_map.csv` — per-probe COSOPT period, pMMC-β, rhythmicity calls,
  period change and period class for WT and *toc1-1*
- `cosopt_nonshortening_genes.csv` — the COSOPT "does not shorten" set
- `concordance_summary.txt` — rhythmicity-call agreement, period correlation,
  period-shift agreement and gene-list overlap versus the cosine method
- `cosopt_concordance.png` — period concordance scatters (WT, *toc1-1*) and
  rhythmic-fraction comparison

## Headline concordance

- Genome-wide, COSOPT reproduces the study's effects and the original report:
  WT 16% / *toc1-1* 7% rhythmic, ~66% of WT-rhythmic transcripts lost, median
  period 24.0 → 22.0 h (cosine: 19% / 11%, 58% lost, 23.8 → 21.8 h).
- Per-gene periods are highly concordant between methods (Pearson r = 0.93 in WT,
  0.96 in *toc1-1*; median absolute difference 0.1 h), and the median period shift
  is identical (−2.1 h COSOPT, −2.0 h cosine).
- Rhythmicity calls agree for ~91% (WT) / 94% (*toc1-1*) of probes (Cohen's
  κ ≈ 0.6–0.7), with COSOPT the more conservative caller.
- The headline period-invariant calcium genes *ACA8* and *PLP2* are classified as
  period-invariant by both methods.
