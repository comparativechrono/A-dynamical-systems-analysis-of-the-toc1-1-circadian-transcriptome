#!/usr/bin/env python3
"""
Methodological sensitivity analysis for the transcript<->[Ca2+]cyt nu-gap.

The main analysis uses the more stringent choices: the full-spectrum Vinnicombe
nu-gap WITH the winding-number condition, and a reliability gate based on the
simulation R^2 (> 0.5 in both genotypes). The DyDE implementation (Mombaerts
et al. 2019) instead (i) measures the nu-gap only over the oscillation frequency
band and skips the winding-number test, and (ii) additionally rejects models
whose absolute DC gain is <= 0.1. This script recomputes the gene<->[Ca2+]cyt
results under the DyDE-style choices and quantifies the concordance, to show the
candidate set is robust to these implementation differences.

Outputs: figS4_nugap_method_sensitivity.png and nugap_sensitivity_summary.txt

Usage: python nugap_method_sensitivity.py [E1.xlsx E2.xlsx calcium.xlsx gene_list.csv] [outdir]
"""
import sys, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon, pearsonr

# find the nugap package (installed, or bundled alongside this script)
try:
    from nugap.fitting import fit_first_order, dc_gain
    from nugap.metric import nu_gap
except ModuleNotFoundError:
    here = os.path.dirname(os.path.abspath(__file__))
    for c in (os.path.join(here, "nugap_package", "src"), os.path.join(here, "..", "nugap", "src")):
        if os.path.isdir(c):
            sys.path.insert(0, c); break
    from nugap.fitting import fit_first_order, dc_gain
    from nugap.metric import nu_gap

TP = [49, 53, 57, 61, 65, 69, 73, 77, 81, 85, 89, 93]
T = np.array(TP, float)
DT = 4.0
# circadian oscillation band for the band-limited (DyDE-style) nu-gap: periods 16-36 h
BAND = (2 * np.pi / 36.0, 2 * np.pi / 16.0)
DCFLOOR = 0.1  # DyDE low-DC-gain rejection threshold
zc = lambda x: (x - x.mean()) / (x.std() + 1e-9)


def main(f1, f2, fca, gene_csv, outdir="output"):
    os.makedirs(outdir, exist_ok=True)
    op = lambda n: os.path.join(outdir, n)
    d1 = pd.read_excel(f1, sheet_name="All").set_index("Name")
    d2 = pd.read_excel(f2, sheet_name="All").set_index("Name").reindex(d1.index)
    WT = np.stack([d1[[f"WT1_{x}" for x in TP]].to_numpy(), d2[[f"WT2_{x}" for x in TP]].to_numpy()], axis=1)
    TC = np.stack([d1[[f"T1_{x}" for x in TP]].to_numpy(), d2[[f"T2_{x}" for x in TP]].to_numpy()], axis=1)
    ca = pd.read_excel(fca)
    caWT, caTC = ca["wt"].to_numpy(float), ca["toc1-1"].to_numpy(float)
    names = d1.index.to_numpy()
    WTm, TCm = WT.mean(1), TC.mean(1)
    czW, czT = zc(caWT), zc(caTC)
    idx = {n: i for i, n in enumerate(names)}

    # the transcripts carried into the calcium analysis (list + main broadband nu-gaps)
    gc = pd.read_csv(gene_csv)

    def bb(s1, s2):   # main: broadband + winding
        return nu_gap(s1, s2, n=256)
    def bl(s1, s2):   # DyDE-style: band-limited, no winding
        return nu_gap(s1, s2, n=400, band=BAND, check_winding=False)

    rows = []
    for _, r in gc.iterrows():
        i = idx[r["name"]]
        # forward transcript -> [Ca2+]cyt
        fW, _ = fit_first_order(T, czW, zc(WTm[i])); fT, _ = fit_first_order(T, czT, zc(TCm[i]))
        f_bb, f_bl = bb(fW, fT), bl(fW, fT)
        fwi = np.median([bl(fit_first_order(T, czW, zc(WT[i, 0]))[0], fit_first_order(T, czW, zc(WT[i, 1]))[0]),
                         bl(fit_first_order(T, czT, zc(TC[i, 0]))[0], fit_first_order(T, czT, zc(TC[i, 1]))[0])])
        # reverse [Ca2+]cyt -> transcript
        vW, _ = fit_first_order(T, zc(WTm[i]), czW); vT, _ = fit_first_order(T, zc(TCm[i]), czT)
        r_bb, r_bl = bb(vW, vT), bl(vW, vT)
        rows.append((r["name"], f_bb, f_bl, fwi, r_bb, r_bl, min(dc_gain(fW), dc_gain(fT))))
    b = pd.DataFrame(rows, columns=["name", "fwd_bb", "fwd_bl", "fwd_bl_within",
                                    "rev_bb", "rev_bl", "fwd_dcgain"]).merge(gc, on="name")

    relf = b[b["fwd_reliable"]]
    relr = b[b["rev_reliable"]]
    r_fwd = pearsonr(relf["fwd_bb"], relf["fwd_bl"])[0]
    r_rev = pearsonr(relr["rev_bb"], relr["rev_bl"])[0]
    p_bl = wilcoxon(relf["fwd_bl"], relf["fwd_bl_within"], alternative="greater").pvalue
    b["robust_bl"] = (b["fwd_reliable"] & b["rev_reliable"] & (b["fwd_bl"] < 0.15) & (b["rev_bl"] < 0.15))
    main_set = set(b.loc[b["robust_both"], "name"])
    bl_set = set(b.loc[b["robust_bl"], "name"])
    rr = b[b["robust_both"]]

    L = []
    def log(s): print(s); L.append(s)
    log("Methodological sensitivity: broadband+winding (main) vs band-limited, no-winding (DyDE-style)")
    log(f"  forward transcript->Ca  (n={len(relf)} reliable): Pearson r(broadband, band) = {r_fwd:.2f}")
    log(f"  reverse  Ca->transcript (n={len(relr)} reliable): Pearson r(broadband, band) = {r_rev:.2f}")
    log(f"  band-limited forward: median between = {relf['fwd_bl'].median():.2f}, within = {relf['fwd_bl_within'].median():.2f}, "
        f"between>within Wilcoxon p = {p_bl:.1e}")
    log(f"  robust-in-both candidates: main (broadband) = {len(main_set)}, band-limited = {len(bl_set)}, "
        f"overlap = {len(main_set & bl_set)}  (all main candidates retained: {main_set <= bl_set})")
    log(f"  DC gain of the {len(rr)} main candidates (forward models): "
        f"min = {rr['fwd_dcgain'].min():.2f}, median = {rr['fwd_dcgain'].median():.2f}; "
        f"number at/below DyDE floor ({DCFLOOR}) = {(rr['fwd_dcgain'] <= DCFLOOR).sum()}")
    open(op("nugap_sensitivity_summary.txt"), "w").write("\n".join(L) + "\n")
    b.to_csv(op("nugap_sensitivity_table.csv"), index=False)

    # ---- figure ----
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.8))
    a0 = ax[0]
    a0.scatter(relf["fwd_bb"], relf["fwd_bl"], s=12, alpha=.35, c="steelblue", edgecolor="none", label="reliable transcripts")
    a0.scatter(rr["fwd_bb"], rr["fwd_bl"], s=34, c="crimson", edgecolor="k", lw=.3, label="robust-in-both candidates")
    a0.plot([0, 1], [0, 1], "r--", lw=1)
    a0.set_xlabel("broadband \u03bd-gap, with winding test (main)")
    a0.set_ylabel("band-limited \u03bd-gap, no winding (DyDE-style)")
    a0.set_title(f"Forward transcript\u2192[Ca\u00b2\u207a]cyt \u03bd-gap\n(Pearson r = {r_fwd:.2f}, n = {len(relf)})")
    a0.set_xlim(0, 1); a0.set_ylim(0, 1); a0.legend(fontsize=8); a0.grid(alpha=.3)

    a1 = ax[1]
    a1.hist(relf["fwd_dcgain"], bins=30, color="lightsteelblue", edgecolor="white", label="all reliable")
    a1.hist(rr["fwd_dcgain"], bins=30, color="crimson", alpha=.8, label="robust-in-both candidates")
    a1.axvline(DCFLOOR, color="k", ls="--", lw=1.2, label=f"DyDE DC-gain floor ({DCFLOOR})")
    a1.set_xlabel("forward-model DC gain |H(1)|")
    a1.set_ylabel("number of transcripts")
    a1.set_title("Input\u2192output DC gain\n(all candidates lie far above the floor)")
    a1.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(op("figS4_nugap_method_sensitivity.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("[done] wrote figS4_nugap_method_sensitivity.png and nugap_sensitivity_summary.txt")


if __name__ == "__main__":
    a = sys.argv[1:]
    if len(a) >= 4:
        main(a[0], a[1], a[2], a[3], a[4] if len(a) > 4 else "output")
    else:
        D = "data"
        main(f"{D}/100725E1c24.xlsx", f"{D}/100725E2c24.xlsx",
             f"{D}/input_calcium_data.xlsx", f"{D}/table_gene_calcium.csv", "output")
