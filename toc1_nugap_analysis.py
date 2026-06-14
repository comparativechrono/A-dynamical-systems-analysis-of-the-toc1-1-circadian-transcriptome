#!/usr/bin/env python3
"""
Reproducible analysis for:
  "TOC1 shortens the circadian transcriptome but a TOC1-independent subset
   tracks the unchanged [Ca2+]cyt rhythm in Arabidopsis"

Inputs (ATH1 microarray, two replicate experiments E1/E2, plus a concurrent
[Ca2+]cyt timecourse), all sampled every 4 h from 49-93 h in constant light
(12 timepoints over circadian cycles 3-4), wild type (C24) and toc1-1.

Each microarray file, sheet "All":
  cols 0-2 : Name (Affymetrix probe id), Locus Identifier, Annotation
  cols 3-14: WT{rep}_49..93   cols 15-26: T{rep}_49..93   (raw log2 block)
  cols 27-50: same headers again -> median-centred block (we use the raw block;
              all dynamic methods mean-centre per trajectory themselves).
The two replicate files are NOT row-aligned -> merge on probe Name.

Usage:  python toc1_nugap_analysis.py E1.xlsx E2.xlsx calcium.xlsx [outdir]

Requires: numpy, scipy, pandas, matplotlib, and the `nugap` package
(pip install nugap, or add its src/ to PYTHONPATH).
"""
import sys, os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, wilcoxon

# ---- locate nugap ----------------------------------------------------------
try:
    from nugap.fitting import fit_first_order
    from nugap.metric import nu_gap
    from nugap.network import compare_network
except ImportError:
    # fall back to a copy of the package bundled alongside this script
    sys.modules.pop("nugap", None)
    _here = os.path.dirname(os.path.abspath(__file__))
    for _p in (os.path.join(_here, "nugap_package", "src"),
               os.path.join(_here, "nugap", "src"), "/home/claude/nugap/src"):
        if os.path.isdir(_p):
            sys.path.insert(0, _p); break
    from nugap.fitting import fit_first_order
    from nugap.metric import nu_gap
    from nugap.network import compare_network

TP = [49, 53, 57, 61, 65, 69, 73, 77, 81, 85, 89, 93]
T = np.array(TP, float)
CLOCK = {'AT5G61380': 'TOC1', 'AT2G46830': 'CCA1', 'AT1G01060': 'LHY',
         'AT2G46790': 'PRR9', 'AT5G02810': 'PRR7', 'AT5G24470': 'PRR5',
         'AT1G22770': 'GI', 'AT2G25930': 'ELF3', 'AT2G40080': 'ELF4',
         'AT3G46640': 'LUX', 'AT1G68050': 'FKF1'}
zc = lambda x: (x - x.mean()) / (x.std() + 1e-9)


# ---- data loading ----------------------------------------------------------
def load(f1, f2, fca):
    d1 = pd.read_excel(f1, sheet_name="All").set_index('Name')
    d2 = pd.read_excel(f2, sheet_name="All").set_index('Name').reindex(d1.index)
    WT = np.stack([d1[[f"WT1_{x}" for x in TP]].to_numpy(),
                   d2[[f"WT2_{x}" for x in TP]].to_numpy()], axis=1)
    TC = np.stack([d1[[f"T1_{x}" for x in TP]].to_numpy(),
                   d2[[f"T2_{x}" for x in TP]].to_numpy()], axis=1)
    ca = pd.read_excel(fca)
    return (WT, TC, d1.index.to_numpy(),
            d1['Locus Identifier'].astype(str).to_numpy(),
            d1['Annotation'].astype(str).to_numpy(),
            ca['wt'].to_numpy(float), ca['toc1-1'].to_numpy(float))


# ---- cosine rhythm fit (period, amplitude, R2, peak time) ------------------
def cosfit(Y):
    """Y: genes x 12. Returns period, amplitude, R2 at best period."""
    Yc = Y - Y.mean(1, keepdims=True)
    sst = (Yc ** 2).sum(1)
    bR = np.full(len(Y), -9.0); bP = np.full(len(Y), np.nan); bA = np.zeros(len(Y))
    for P in np.arange(18, 30.01, 0.25):
        c = np.cos(2 * np.pi * T / P); s = np.sin(2 * np.pi * T / P)
        X = np.stack([c - c.mean(), s - s.mean()], 1)
        b, *_ = np.linalg.lstsq(X, Yc.T, rcond=None)
        r2 = 1 - ((Yc - (X @ b).T) ** 2).sum(1) / np.maximum(sst, 1e-9)
        amp = np.hypot(b[0], b[1]); u = r2 > bR
        bR[u] = r2[u]; bP[u] = P; bA[u] = amp[u]
    return bP, bA, bR


def cosfit1(y):
    P, A, R = cosfit(y[None, :]); return P[0], A[0], R[0]


def main(f1, f2, fca, outdir="."):
    os.makedirs(outdir, exist_ok=True)
    op = lambda n: os.path.join(outdir, n)
    WT, TC, names, locus, annot, caWT, caTC = load(f1, f2, fca)
    locU = np.char.upper(locus.astype(str))
    WTm, TCm = WT.mean(1), TC.mean(1)
    n = len(names)
    print(f"[load] {n} probes; replicate concordance "
          f"WT={np.nanmean([np.corrcoef(WT[:,0,k],WT[:,1,k])[0,1] for k in range(12)]):.3f} "
          f"toc1={np.nanmean([np.corrcoef(TC[:,0,k],TC[:,1,k])[0,1] for k in range(12)]):.3f}")

    # ---------- rhythm calls (genome-wide) ----------
    pW, aW, rW = cosfit(WTm); pT, aT, rT = cosfit(TCm)
    expressed = WT.mean((1, 2)) > 5.0
    TH = 0.6
    rhyW = expressed & (rW > TH); rhyT = expressed & (rT > TH)
    stay = rhyW & rhyT; lost = rhyW & ~rhyT
    print(f"[rhythm] WT rhythmic={rhyW.sum()} ({100*rhyW.sum()/n:.0f}%)  "
          f"toc1 rhythmic={rhyT.sum()} ({100*rhyT.sum()/n:.0f}%)  "
          f"lost in toc1={100*lost.sum()/rhyW.sum():.0f}%")
    print(f"[rhythm] period WT={np.median(pW[rhyW]):.1f}h toc1={np.median(pT[rhyT]):.1f}h; "
          f"amp(stay)={np.median(aW[stay]):.2f} amp(lost)={np.median(aW[lost]):.2f} "
          f"MWU p={mannwhitneyu(aW[stay],aW[lost],alternative='greater').pvalue:.1e}")

    # ---------- calcium rhythm ----------
    pcaW, _, rcaW = cosfit1(caWT); pcaT, _, rcaT = cosfit1(caTC)
    print(f"[calcium] period WT={pcaW:.1f}h (R2={rcaW:.2f})  toc1={pcaT:.1f}h (R2={rcaT:.2f})")

    # ---------- strict rhythmic set + period-invariant subset ----------
    repr_ = np.array([np.corrcoef(WT[i,0],WT[i,1])[0,1] for i in range(n)]) > 0.6
    rhythmic = expressed & repr_ & (rW > 0.7) & (aW > 0.7)
    dP = pT - pW
    # classify period change among genes confidently rhythmic in BOTH genotypes
    still = rhythmic & (rT > 0.6) & (aT > 0.6)        # still rhythmic in toc1
    invariant = still & (np.abs(dP) < 0.6)            # period essentially unchanged
    lengthened = still & (dP >= 0.6)                  # period increases
    shortened = still & (dP <= -0.6)                  # period decreases (the bulk)
    does_not_shorten = invariant | lengthened
    print(f"[period] rhythmic in both={still.sum()}  invariant(|dP|<0.6h)={invariant.sum()}  "
          f"lengthened(dP>=+0.6h)={lengthened.sum()}  shortened={shortened.sum()}  "
          f"-> do NOT shorten={does_not_shorten.sum()}")

    pclass = np.full(len(names), "", dtype=object)
    pclass[invariant] = "invariant"; pclass[lengthened] = "lengthened"; pclass[shortened] = "shortened"
    dfm = pd.DataFrame({'name': names, 'locus': locus, 'annot': annot,
                        'expr_mean': WT.mean((1,2)),
                        'period_WT': pW, 'period_toc1': pT, 'dperiod': dP,
                        'amp_WT': aW, 'amp_toc1': aT, 'cosR2_WT': rW, 'cosR2_toc1': rT,
                        'rhythmic_WT': rhyW, 'rhythmic_toc1': rhyT,
                        'rhythmic_strict': rhythmic, 'rhythmic_both': still,
                        'period_class': pclass,
                        'period_invariant': invariant, 'period_lengthened': lengthened,
                        'does_not_shorten': does_not_shorten})
    dfm.to_csv(op("table_period_map.csv"), index=False)
    dfm[does_not_shorten].sort_values('dperiod').to_csv(op("table_nonshortening_genes.csv"), index=False)
    # Supplementary Table S5: full per-transcript period analysis (cosine method, all probes)
    dfm['expressed'] = expressed
    s5cols = ['locus','name','annot','expr_mean','expressed',
              'period_WT','amp_WT','cosR2_WT','rhythmic_WT',
              'period_toc1','amp_toc1','cosR2_toc1','rhythmic_toc1','dperiod','period_class']
    dfm[s5cols].sort_values(['rhythmic_WT','expressed','cosR2_WT'],
        ascending=[False, False, False]).to_csv(op("suppl_table_S5_period_analysis.csv"), index=False)

    # ---------- gene <-> calcium nu-gap (both input/output directions) ----------
    czW, czT = zc(caWT), zc(caTC)
    def fwd(g, cal):   # transcript -> [Ca2+]cyt : output=Ca, input=gene
        m, r2 = fit_first_order(T, cal, zc(g)); return m, r2
    def rev(g, cal):   # [Ca2+]cyt -> transcript : output=gene, input=Ca
        m, r2 = fit_first_order(T, zc(g), cal); return m, r2
    ridx = np.where(rhythmic)[0]
    rows = []
    for i in ridx:
        # forward
        fW, frW = fwd(WTm[i], czW); fT, frT = fwd(TCm[i], czT)
        fgap = nu_gap(fW, fT, n=256)
        fwi = np.median([nu_gap(fwd(WT[i,0],czW)[0], fwd(WT[i,1],czW)[0], n=256),
                         nu_gap(fwd(TC[i,0],czT)[0], fwd(TC[i,1],czT)[0], n=256)])
        # reverse
        vW, vrW = rev(WTm[i], czW); vT, vrT = rev(TCm[i], czT)
        rgap = nu_gap(vW, vT, n=256)
        rwi = np.median([nu_gap(rev(WT[i,0],czW)[0], rev(WT[i,1],czW)[0], n=256),
                         nu_gap(rev(TC[i,0],czT)[0], rev(TC[i,1],czT)[0], n=256)])
        rows.append((names[i], locus[i], annot[i], frW, frT, fgap, fwi, vrW, vrT, rgap, rwi))
    gc = pd.DataFrame(rows, columns=['name','locus','annot',
        'fwd_r2_WT','fwd_r2_toc1','fwd_nugap','fwd_within',
        'rev_r2_WT','rev_r2_toc1','rev_nugap','rev_within'])
    gc = gc.merge(dfm[['name','period_WT','period_toc1','dperiod','period_class','does_not_shorten']], on='name')
    gc['fwd_reliable'] = (gc['fwd_r2_WT'] > 0.5) & (gc['fwd_r2_toc1'] > 0.5)
    gc['rev_reliable'] = (gc['rev_r2_WT'] > 0.5) & (gc['rev_r2_toc1'] > 0.5)
    gc['robust_both'] = (gc['fwd_reliable'] & gc['rev_reliable'] &
                         (gc['fwd_nugap'] < 0.15) & (gc['rev_nugap'] < 0.15))
    gc.to_csv(op("table_gene_calcium.csv"), index=False)
    relf = gc[gc['fwd_reliable']]; relr = gc[gc['rev_reliable']]
    rho = gc.loc[gc['fwd_reliable'] & gc['rev_reliable'], ['fwd_nugap','rev_nugap']].corr().iloc[0,1]
    print(f"[fwd transcript->Ca] reliable={len(relf)} nugap median={relf['fwd_nugap'].median():.2f} "
          f"within={relf['fwd_within'].median():.2f} "
          f"Wilcoxon b>w p={wilcoxon(relf['fwd_nugap'],relf['fwd_within'],alternative='greater').pvalue:.1e}")
    print(f"[rev Ca->transcript] reliable={len(relr)} nugap median={relr['rev_nugap'].median():.2f} "
          f"within={relr['rev_within'].median():.2f}")
    print(f"[bidirectional] corr(fwd,rev nugap)={rho:.2f}  robust-in-both={int(gc['robust_both'].sum())}")

    # ----- supplementary tables -----
    # S1: non-shortening transcripts with both-direction nu-gaps
    s1 = dfm[does_not_shorten].merge(
        gc[['name','fwd_nugap','fwd_reliable','rev_nugap','rev_reliable']], on='name', how='left')
    s1.sort_values(['period_class','fwd_nugap']).to_csv(op("suppl_table_S1_nonshortening.csv"), index=False)
    # S2: robust-in-both calcium-coupled candidates
    gc[gc['robust_both']].sort_values('fwd_nugap').to_csv(op("suppl_table_S2_robust_both.csv"), index=False)
    # S3: full forward gene->Ca reliable
    relf.sort_values('fwd_nugap')[['name','locus','annot','period_WT','period_toc1','dperiod',
        'fwd_nugap','fwd_within','fwd_r2_WT','fwd_r2_toc1']].to_csv(op("suppl_table_S3_gene_to_calcium.csv"), index=False)
    # S4: full reverse Ca->gene reliable
    relr.sort_values('rev_nugap')[['name','locus','annot','period_WT','period_toc1','dperiod',
        'rev_nugap','rev_within','rev_r2_WT','rev_r2_toc1']].to_csv(op("suppl_table_S4_calcium_to_gene.csv"), index=False)

    # ---------- nu-gap interaction network (distributed-change demonstration) ----------
    sel = list(dict.fromkeys(list(np.argsort(np.where(rhythmic, rW, -1))[::-1][:150])
               + [int(np.where(np.char.find(locU,l)>=0)[0][0]) for l in CLOCK
                  if np.any(np.char.find(locU,l)>=0)]))
    dWT = {names[i]: WT[i] for i in sel}; dTC = {names[i]: TC[i] for i in sel}
    net = compare_network(dWT, dTC, T, order=1, gate="either", min_r2=0.5, progress=False)
    net.to_csv(op("table_network_edges.csv"), index=False)
    above = (net['nu_gap'] > net['within_median']).mean()
    print(f"[network] {len(net)} edges, {int((net['q_global']<0.1).sum())} significant at q<0.1; "
          f"between>within in {above*100:.0f}% (global Wilcoxon p="
          f"{wilcoxon(net['nu_gap'],net['within_median'],alternative='greater').pvalue:.1e})")

    _figures(outdir, op, T, WTm, TCm, caWT, caTC, locU, dfm, gc, net,
             rhyW, rhyT, stay, lost, pW, pT, aW)
    print("[done] figures + tables written to", outdir)
    return dfm, gc, net


def _figures(outdir, op, T, WTm, TCm, caWT, caTC, locU, dfm, gc, net,
             rhyW, rhyT, stay, lost, pW, pT, aW):
    GR = {'AT5G61380':'TOC1','AT2G46830':'CCA1','AT1G01060':'LHY','AT1G22770':'GI'}
    # Fig 1: clock genes + calcium
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    for a,(loc,sym) in zip(ax[:2], list(GR.items())[:2]):
        i = int(np.where(np.char.find(locU, loc) >= 0)[0][0])
        a.plot(T, WTm[i]-WTm[i].mean(), '-o', color='seagreen', ms=3, label='WT (C24)')
        a.plot(T, TCm[i]-TCm[i].mean(), '-s', color='crimson', ms=3, label='toc1-1')
        a.set_title(f"{sym}"); a.set_xlabel("time in LL (h)"); a.grid(alpha=.3)
    ax[0].set_ylabel("mean-centred log2 expression"); ax[0].legend(fontsize=8)
    a = ax[2]
    a.plot(T, caWT, '-o', color='seagreen', ms=3, label='WT')
    a.plot(T, caTC, '-s', color='crimson', ms=3, label='toc1-1')
    a.set_title("[Ca2+]cyt (aequorin)"); a.set_xlabel("time in LL (h)")
    a.set_ylabel("luminescence"); a.legend(fontsize=8); a.grid(alpha=.3)
    fig.suptitle("Clock transcripts shorten period in toc1-1; [Ca2+]cyt rhythm is unchanged")
    fig.tight_layout(); fig.savefig(op("fig1_clock_calcium.png"), dpi=130, bbox_inches="tight"); plt.close(fig)

    # Fig 2: genome-wide period histogram + amplitude
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    ax[0].hist(pW[rhyW], bins=np.arange(18,30,0.5), color='seagreen', alpha=.6, label=f'WT (n={rhyW.sum()})')
    ax[0].hist(pT[rhyT], bins=np.arange(18,30,0.5), color='crimson', alpha=.6, label=f'toc1-1 (n={rhyT.sum()})')
    ax[0].axvline(np.median(pW[rhyW]),color='seagreen',ls='--'); ax[0].axvline(np.median(pT[rhyT]),color='crimson',ls='--')
    ax[0].set_xlabel("period (h)"); ax[0].set_ylabel("rhythmic genes"); ax[0].set_title("Period shifts shorter in toc1-1"); ax[0].legend(fontsize=8)
    ax[1].bar(['WT','toc1-1'],[rhyW.sum(),rhyT.sum()],color=['seagreen','crimson'])
    ax[1].set_ylabel("# rhythmic genes"); ax[1].set_title("Fewer rhythmic genes in toc1-1")
    ax[2].boxplot([aW[stay],aW[lost]],labels=['stay\nrhythmic','become\narrhythmic'],showfliers=False)
    ax[2].set_ylabel("WT amplitude"); ax[2].set_title("Genes lost in toc1-1 have lower amplitude")
    fig.tight_layout(); fig.savefig(op("fig2_genomewide.png"), dpi=130, bbox_inches="tight"); plt.close(fig)

    # Fig 3: nu-gap network distributed signal
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    ax[0].scatter(net['within_median'], net['nu_gap'], s=5, alpha=.15, c='gray', edgecolor='none')
    ax[0].plot([0,1],[0,1],'r--',lw=1.2,label='between = within')
    ax[0].set_xlabel("within-condition nu-gap (noise floor)"); ax[0].set_ylabel("between-condition nu-gap")
    ax[0].set_xlim(0,1); ax[0].set_ylim(0,1); ax[0].legend(fontsize=8)
    ax[0].set_title(f"Network change is global ({(net['nu_gap']>net['within_median']).mean()*100:.0f}% above diagonal)\nbut no single edge survives FDR")
    ax[1].scatter(gc.loc[gc['fwd_reliable'],'fwd_within'], gc.loc[gc['fwd_reliable'],'fwd_nugap'], s=10, alpha=.5, c='steelblue', edgecolor='none')
    ax[1].plot([0,1],[0,1],'r--',lw=1.2)
    ax[1].set_xlabel("within-condition gene->Ca nu-gap"); ax[1].set_ylabel("between-condition gene->Ca nu-gap")
    ax[1].set_xlim(0,1); ax[1].set_ylim(0,1)
    ax[1].set_title("transcript -> [Ca2+]cyt transfer functions:\nlow = relationship preserved in toc1-1")
    fig.tight_layout(); fig.savefig(op("fig3_nugap.png"), dpi=130, bbox_inches="tight"); plt.close(fig)

    # Fig 4: all genes that do NOT shorten in toc1-1 (invariant + lengthened)
    sub = dfm[dfm['does_not_shorten']].sort_values('dperiod').reset_index()
    k = len(sub); cols = 5; rows = int(np.ceil(k/cols))
    fig, axes = plt.subplots(rows, cols, figsize=(3*cols, 2.4*rows), sharex=True)
    axflat = np.atleast_1d(axes).ravel()
    for ax_, (_, r) in zip(axflat, sub.iterrows()):
        i = int(r['index'])
        ax_.plot(T, WTm[i]-WTm[i].mean(), '-o', color='seagreen', ms=2.5)
        ax_.plot(T, TCm[i]-TCm[i].mean(), '-s', color='crimson', ms=2.5)
        tag = "" if r['period_class'] == "invariant" else "  [lengthened]"
        tc = 'black' if r['period_class'] == "invariant" else 'darkorange'
        ax_.set_title(f"{r['locus']}{tag}\n{str(r['annot'])[:22]}", fontsize=6.5, color=tc)
        ax_.tick_params(labelsize=6)
    for ax_ in axflat[k:]: ax_.axis('off')
    ninv = int((sub['period_class'] == 'invariant').sum()); nlen = k - ninv
    fig.suptitle(f"All {k} transcripts whose period does NOT shorten in toc1-1 "
                 f"({ninv} invariant, {nlen} lengthened): WT (green) vs toc1-1 (red)", y=1.0)
    fig.tight_layout(); fig.savefig(op("fig4_invariant_genes.png"), dpi=130, bbox_inches="tight"); plt.close(fig)

    # Fig S1: forward vs reverse gene<->Ca nu-gap (directions are largely independent)
    bo = gc['fwd_reliable'] & gc['rev_reliable']
    fig, ax = plt.subplots(figsize=(5.2, 5))
    ax.scatter(gc.loc[bo,'fwd_nugap'], gc.loc[bo,'rev_nugap'], s=14, alpha=.5, c='slategray', edgecolor='none')
    rb = gc['robust_both']
    ax.scatter(gc.loc[rb,'fwd_nugap'], gc.loc[rb,'rev_nugap'], s=30, c='crimson', edgecolor='k', lw=.3, label='robust in both (<0.15)')
    ax.axvline(0.15, color='gray', ls=':', lw=.8); ax.axhline(0.15, color='gray', ls=':', lw=.8)
    rho = gc.loc[bo, ['fwd_nugap','rev_nugap']].corr().iloc[0,1]
    ax.set_xlabel("transcript -> [Ca2+]cyt nu-gap (forward)")
    ax.set_ylabel("[Ca2+]cyt -> transcript nu-gap (reverse)")
    ax.set_title(f"Direction matters: forward vs reverse nu-gap\n(reliable in both, Pearson r = {rho:.2f})")
    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.legend(fontsize=8); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(op("figS1_bidirectional.png"), dpi=130, bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__":
    a = sys.argv[1:]
    if len(a) >= 3:
        main(a[0], a[1], a[2], a[3] if len(a) > 3 else ".")
    else:
        U = "/mnt/user-data/uploads/"
        main(U+"100725E1c24.xlsx", U+"100725E2c24.xlsx", U+"input_calcium_data.xlsx", ".")
