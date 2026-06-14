# Calcium-coupling network figure (standalone)

This is a **self-contained** bundle for the supplementary network figure only.
It is deliberately separate from the main reproducibility bundle.

## What it shows

A hub-and-spoke network with cytosolic free calcium ([Ca²⁺]cyt) at the centre
and, as spokes, the transcripts whose dynamic coupling to [Ca²⁺]cyt is preserved
in *toc1-1* in **both** modelling directions (forward transcript→[Ca²⁺]cyt and
reverse [Ca²⁺]cyt→transcript ν-gap < 0.15, with R² > 0.5 in both genotypes in
both directions — the set in Supplementary Table S2).

- **Edge thickness** = strength of the preserved coupling (thicker = lower
  mean ν-gap = more strongly preserved).
- **Node colour** = what happened to each transcript's *own* circadian period in
  *toc1-1*: period-invariant (teal), shortened (orange), or not confidently
  rhythmic in *toc1-1* so unclassified (grey).
- **Red ring** = transcript with a directly calcium/metal-ion-related annotation
  (PLP2, the EF-hand protein AT5G54130, ATFP6).

The figure's message is that coupling to [Ca²⁺]cyt can be preserved even when a
transcript's own period shortens with the accelerated clock — most spokes are
orange (e.g. SIG5), while only PLP2, CYP81F2 and AT1G52200 are period-invariant.

## Contents

```
network_figure.py                       the plotting script
figure_legend.txt                       the figure legend (caption) text
data/gene_calcium_bidirectional.csv     per-transcript forward/reverse gene↔[Ca²⁺]cyt
                                         ν-gap, fit reliability, period class, robust_both flag
output/calcium_coupling_network.png     rendered figure (raster)
output/calcium_coupling_network.svg     rendered figure (vector, for editing/print)
```

## Requirements & run

```bash
pip install pandas numpy matplotlib networkx
python network_figure.py                       # uses data/ and writes to output/
# or:  python network_figure.py data/gene_calcium_bidirectional.csv output
```

## Data provenance

`data/gene_calcium_bidirectional.csv` is the per-transcript forward and reverse
gene↔[Ca²⁺]cyt ν-gap table produced by the main analysis (it corresponds to the
combined content of Supplementary Tables S2–S4). It is included here so this
figure regenerates without the full microarray/calcium pipeline. Columns:
`locus, name, annot, period_WT, period_toc1, dperiod, period_class,
fwd_nugap, fwd_reliable, rev_nugap, rev_reliable, robust_both`.
