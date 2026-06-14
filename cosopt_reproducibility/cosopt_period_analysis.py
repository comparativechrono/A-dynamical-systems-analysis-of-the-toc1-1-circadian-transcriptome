#!/usr/bin/env python3
"""
COSOPT cross-check for the toc1-1 / [Ca2+]cyt study.

The main analysis calls circadian period and rhythmicity with a variable-period
cosine fit. This script substitutes the ORIGINAL COSOPT period analysis into the
same downstream classification (rhythmic in each genotype -> period change ->
period-invariant / shortened / lengthened -> "does not shorten" gene list) and
quantifies the concordance with the cosine method, to show the conclusions do
not depend on the choice of period-analysis tool.

COSOPT only affects the period/rhythmicity layer; the nu-gap analyses use the raw
expression time-courses and are therefore unchanged by this substitution.

Inputs (two replicate COSOPT files, sheet "All", merged on probe Name; per WT and
T=toc1: Amp, P/A, MEL, Per (period, h), Phase, pMMC (pMMC-beta)):
  100725E1c24_cosopt.xlsx, 100725E2c24_cosopt.xlsx
plus the cosine period map for concordance:
  cosine_period_map.csv   (= results/table_period_map.csv from the main pipeline)

Usage: python cosopt_period_analysis.py E1_cosopt.xlsx E2_cosopt.xlsx cosine_period_map.csv [outdir]

Rhythmicity: COSOPT pMMC-beta < 0.1 in BOTH replicates (the criterion used in the
original study). Period = mean of the two replicate COSOPT period estimates.
"""
import sys, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr


def main(f1, f2, cosine_csv, outdir="."):
    os.makedirs(outdir, exist_ok=True)
    op = lambda n: os.path.join(outdir, n)
    e1 = pd.read_excel(f1, sheet_name="All").set_index("Name")
    e2 = pd.read_excel(f2, sheet_name="All").set_index("Name")
    m = e1.join(e2[[c for c in e2.columns if c not in ('Locus Identifier', 'Annotation')]], how='inner')

    # --- COSOPT period + rhythmicity (pMMC-beta < 0.1 in both replicates) ---
    m['period_WT'] = m[['WT1Per', 'WT2Per']].mean(axis=1)
    m['period_toc1'] = m[['T1Per', 'T2Per']].mean(axis=1)
    m['dperiod'] = m['period_toc1'] - m['period_WT']
    m['amp_WT'] = m[['WT1Amp', 'WT2Amp']].mean(axis=1)
    m['rhythmic_WT'] = (m['WT1pMMC'] < 0.1) & (m['WT2pMMC'] < 0.1)
    m['rhythmic_toc1'] = (m['T1pMMC'] < 0.1) & (m['T2pMMC'] < 0.1)

    still = m['rhythmic_WT'] & m['rhythmic_toc1']
    inv = still & (m['dperiod'].abs() < 0.6)
    leng = still & (m['dperiod'] >= 0.6)
    short = still & (m['dperiod'] <= -0.6)
    cls = np.where(inv, 'invariant', np.where(leng, 'lengthened', np.where(short, 'shortened', '')))
    m['period_class'] = cls
    m['does_not_shorten'] = inv | leng

    nW, nT = int(m['rhythmic_WT'].sum()), int(m['rhythmic_toc1'].sum())
    print(f"[COSOPT] rhythmic WT={nW} ({100*nW/len(m):.0f}%)  toc1={nT} ({100*nT/len(m):.0f}%)  "
          f"lost={100*(1-(still.sum()/nW)):.0f}%")
    print(f"[COSOPT] period median WT={m.loc[m['rhythmic_WT'],'period_WT'].median():.1f}h  "
          f"toc1={m.loc[m['rhythmic_toc1'],'period_toc1'].median():.1f}h")
    print(f"[COSOPT] rhythmic in both={int(still.sum())}  invariant={int(inv.sum())}  "
          f"lengthened={int(leng.sum())}  shortened={int(short.sum())}  do-not-shorten={int(m['does_not_shorten'].sum())}")

    keep = ['Locus Identifier', 'Annotation', 'period_WT', 'period_toc1', 'dperiod',
            'amp_WT', 'WT1pMMC', 'WT2pMMC', 'T1pMMC', 'T2pMMC',
            'rhythmic_WT', 'rhythmic_toc1', 'period_class', 'does_not_shorten']
    out = m[keep].rename(columns={'Locus Identifier': 'locus', 'Annotation': 'annot'})
    out.to_csv(op("cosopt_period_map.csv"))
    out[out['does_not_shorten']].sort_values('dperiod').to_csv(op("cosopt_nonshortening_genes.csv"))

    # --- concordance with the cosine method ---
    cos = pd.read_csv(cosine_csv).set_index('name')
    j = out.join(cos[['period_WT', 'period_toc1', 'dperiod', 'rhythmic_WT', 'rhythmic_toc1',
                      'does_not_shorten', 'period_class']], how='inner', lsuffix='_cosopt', rsuffix='_cosine')
    lines = []
    def log(s): print(s); lines.append(s)
    log("\n=== CONCORDANCE: COSOPT vs cosine ===")
    # rhythmicity-call agreement
    for g in ['WT', 'toc1']:
        a = j[f'rhythmic_{g}_cosopt'].astype(bool); b = j[f'rhythmic_{g}_cosine'].astype(bool)
        agree = (a == b).mean()
        inter = (a & b).sum(); uni = (a | b).sum()
        # Cohen's kappa
        po = agree; pe = (a.mean()*b.mean() + (1-a.mean())*(1-b.mean())); kappa = (po-pe)/(1-pe)
        log(f"rhythmic {g}: COSOPT={int(a.sum())} cosine={int(b.sum())} both={int(inter)} "
            f"| agreement={100*agree:.0f}% Jaccard={inter/uni:.2f} kappa={kappa:.2f}")
    # period correlation among genes rhythmic in that genotype by both methods
    for g in ['WT', 'toc1']:
        sel = j[f'rhythmic_{g}_cosopt'].astype(bool) & j[f'rhythmic_{g}_cosine'].astype(bool)
        x = j.loc[sel, f'period_{g}_cosopt']; y = j.loc[sel, f'period_{g}_cosine']
        ok = x.notna() & y.notna()
        r = pearsonr(x[ok], y[ok])[0]; rs = spearmanr(x[ok], y[ok])[0]
        log(f"period {g} (n={int(ok.sum())} rhythmic in both methods): Pearson r={r:.2f} Spearman={rs:.2f} "
            f"median |COSOPT-cosine|={ (x[ok]-y[ok]).abs().median():.1f}h")
    # period-shift direction
    both_both = (j['rhythmic_WT_cosopt']&j['rhythmic_toc1_cosopt']&j['rhythmic_WT_cosine']&j['rhythmic_toc1_cosine']).astype(bool)
    log(f"median period shift (rhythmic-both, both methods): COSOPT={j.loc[both_both,'dperiod_cosopt'].median():+.1f}h "
        f"cosine={j.loc[both_both,'dperiod_cosine'].median():+.1f}h")
    # non-shortening gene-list overlap
    cs = set(j.index[j['does_not_shorten_cosopt']]); cn = set(j.index[j['does_not_shorten_cosine']])
    log(f"\nnon-shortening set: COSOPT={len(cs)} cosine={len(cn)} overlap={len(cs&cn)} "
        f"(Jaccard={len(cs&cn)/len(cs|cn):.2f})")
    # headline calcium genes
    for sym, loc in [('ACA8','AT5G57110'),('PLP2','AT2G26560')]:
        rr = out[out['locus'] == loc]
        if len(rr):
            r = rr.iloc[0]
            log(f"  {sym} ({loc}): COSOPT period {r['period_WT']:.1f}->{r['period_toc1']:.1f}h "
                f"(d={r['dperiod']:+.1f}), class='{r['period_class']}'")
    open(op("concordance_summary.txt"), "w").write("\n".join(lines) + "\n")

    _figure(op, j, m, out)
    print("[done] outputs in", outdir)


def _figure(op, j, m, out):
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    for a, g, lab in [(ax[0], 'WT', 'wild type'), (ax[1], 'toc1', 'toc1-1')]:
        sel = j[f'rhythmic_{g}_cosopt'].astype(bool) & j[f'rhythmic_{g}_cosine'].astype(bool)
        x = j.loc[sel, f'period_{g}_cosine']; y = j.loc[sel, f'period_{g}_cosopt']
        ok = x.notna() & y.notna()
        a.scatter(x[ok], y[ok], s=8, alpha=.25, c='steelblue', edgecolor='none')
        a.plot([18, 30], [18, 30], 'r--', lw=1)
        r = pearsonr(x[ok], y[ok])[0]
        a.set_xlabel("cosine period (h)"); a.set_ylabel("COSOPT period (h)")
        a.set_title(f"{lab}: period concordance\n(n={int(ok.sum())}, r={r:.2f})")
        a.set_xlim(18, 30); a.set_ylim(18, 30); a.grid(alpha=.3)
    # rhythmic fraction bars: cosine vs COSOPT (with original-study reference)
    a = ax[2]
    cosW = j['rhythmic_WT_cosine'].mean()*100; cosT = j['rhythmic_toc1_cosine'].mean()*100
    copW = j['rhythmic_WT_cosopt'].mean()*100; copT = j['rhythmic_toc1_cosopt'].mean()*100
    xpos = np.arange(2); w = 0.35
    a.bar(xpos-w/2, [copW, copT], w, label='COSOPT', color='#5b8c5a')
    a.bar(xpos+w/2, [cosW, cosT], w, label='cosine', color='#c46a4a')
    a.set_xticks(xpos); a.set_xticklabels(['WT', 'toc1-1'])
    a.set_ylabel("% of merged probes rhythmic"); a.set_title("Rhythmic fraction by method")
    a.legend(fontsize=9)
    fig.tight_layout(); fig.savefig(op("cosopt_concordance.png"), dpi=130, bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__":
    a = sys.argv[1:]
    if len(a) >= 3:
        main(a[0], a[1], a[2], a[3] if len(a) > 3 else ".")
    else:
        U = "/mnt/user-data/uploads/"
        main(U+"100725E1c24_cosopt.xlsx", U+"100725E2c24_cosopt.xlsx", "cosine_period_map.csv", ".")
