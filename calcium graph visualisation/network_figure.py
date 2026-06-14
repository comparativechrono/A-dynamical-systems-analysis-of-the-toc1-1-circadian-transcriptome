#!/usr/bin/env python3
"""
Calcium-coupling network figure for the toc1-1 / [Ca2+]cyt nu-gap study.

Draws [Ca2+]cyt as a central hub and the direction-robust calcium-coupled
transcripts (forward AND reverse gene<->[Ca2+]cyt nu-gap < 0.15, reliable in
both directions) as spokes. Edge thickness encodes the strength of the
preserved coupling (thicker = lower nu-gap); node colour encodes what happened
to each transcript's OWN circadian period in toc1-1 (period-invariant vs
shortened). The figure therefore shows that a transcript's coupling to
[Ca2+]cyt can be preserved even when its own period shortens (e.g. SIG5).

Usage:  python network_figure.py [data/gene_calcium_bidirectional.csv] [outdir]
Requires: pandas, numpy, matplotlib, networkx
"""
import sys, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx

# short, readable labels for the candidate loci
SYMBOL = {
    'AT2G41290': 'strictosidine synthase', 'AT2G26560': 'PLP2', 'AT5G67480': 'BT4',
    'AT5G57220': 'CYP81F2', 'AT5G26340': 'MSS1/STP13', 'AT1G69830': 'AMY3',
    'AT4G38580': 'ATFP6', 'AT4G32190': 'centromeric protein', 'AT5G54130': 'EF-hand (AT5G54130)',
    'AT5G24120': 'SIG5', 'AT1G52200': 'AT1G52200', 'AT5G55380': 'MBOAT',
    'AT4G37760': 'SQE3', 'AT1G66330': 'senescence-assoc.',
}
# transcripts whose annotation is directly calcium / metal-ion related (marked with a ring)
CA_RELATED = {'AT2G26560', 'AT5G54130', 'AT4G38580', 'AT5G57110'}
CLASS_COLOUR = {'invariant': '#2a9d8f', 'shortened': '#e76f51', 'lengthened': '#7b2cbf',
                'unclassified': '#b0b0b0'}
CLASS_LABEL = {'invariant': 'period-invariant (|\u0394P| < 0.6 h)',
               'shortened': 'period shortened in toc1-1',
               'lengthened': 'period lengthened in toc1-1',
               'unclassified': 'not confidently rhythmic in toc1-1'}


def main(csv="data/gene_calcium_bidirectional.csv", outdir="output"):
    os.makedirs(outdir, exist_ok=True)
    df = pd.read_csv(csv)
    s = df[df['robust_both']].copy()
    s['period_class'] = s['period_class'].fillna('unclassified')
    s['mean_nugap'] = s[['fwd_nugap', 'rev_nugap']].mean(axis=1)
    # order: group by class, then by strongest coupling (lowest nu-gap) first
    order = {'invariant': 0, 'shortened': 1, 'lengthened': 2, 'unclassified': 3}
    s['ord'] = s['period_class'].map(order)
    s = s.sort_values(['ord', 'mean_nugap']).reset_index(drop=True)
    n = len(s)

    G = nx.Graph()
    G.add_node('Ca')
    pos = {'Ca': (0.0, 0.0)}
    # radial layout, leaving a small gap so class groups read as arcs
    angles = np.linspace(90, 90 + 360, n, endpoint=False) * np.pi / 180.0
    for k, (_, r) in enumerate(s.iterrows()):
        node = r['locus']
        G.add_node(node)
        G.add_edge('Ca', node, w=r['mean_nugap'])
        pos[node] = (np.cos(angles[k]), np.sin(angles[k]))

    fig, ax = plt.subplots(figsize=(9.5, 9.5))
    # edges: width scales with preserved coupling strength (lower nu-gap -> thicker)
    for u, v, d in G.edges(data=True):
        w = 1.0 + (0.16 - d['w']) / 0.16 * 6.0
        ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                color='#9aa0a6', lw=w, zorder=1, solid_capstyle='round')
    # candidate nodes
    for _, r in s.iterrows():
        x, y = pos[r['locus']]
        col = CLASS_COLOUR.get(r['period_class'], '#888')
        ax.scatter([x], [y], s=620, c=col, edgecolors='k', linewidths=0.6, zorder=3)
        if r['locus'] in CA_RELATED:   # ring calcium/metal-related transcripts
            ax.scatter([x], [y], s=1050, facecolors='none', edgecolors='crimson',
                       linewidths=1.8, zorder=4)
        lab = SYMBOL.get(r['locus'], r['locus'])
        ha = 'left' if x > 0.05 else ('right' if x < -0.05 else 'center')
        ax.text(x * 1.18, y * 1.18, lab, fontsize=9, ha=ha, va='center', zorder=5)
    # hub
    ax.scatter([0], [0], s=6800, c='#f4d35e', edgecolors='k', linewidths=1.2, zorder=6)
    ax.text(0, 0, '[Ca$^{2+}$]$_{cyt}$', fontsize=11, ha='center', va='center',
            fontweight='bold', zorder=7)

    # legends
    node_leg = [Line2D([0], [0], marker='o', color='w', label=CLASS_LABEL[c],
                       markerfacecolor=CLASS_COLOUR[c], markeredgecolor='k', markersize=12)
                for c in ['invariant', 'shortened', 'lengthened', 'unclassified'] if (s['period_class'] == c).any()]
    node_leg.append(Line2D([0], [0], marker='o', color='w', label='calcium / metal-ion related',
                           markerfacecolor='none', markeredgecolor='crimson', markersize=14, markeredgewidth=1.8))
    edge_leg = [Line2D([0], [0], color='#9aa0a6', lw=1.0 + (0.16 - q) / 0.16 * 6.0,
                       label=f'\u03bd-gap \u2248 {q:.2f}') for q in (0.04, 0.10, 0.15)]
    l1 = ax.legend(handles=node_leg, loc='upper left', fontsize=8.5, frameon=False,
                   title='node colour = transcript period', title_fontsize=9, bbox_to_anchor=(-0.02, 1.02))
    ax.add_artist(l1)
    ax.legend(handles=edge_leg, loc='lower left', fontsize=8.5, frameon=False,
              title='edge = preserved coupling\n(thicker = stronger)', title_fontsize=9,
              bbox_to_anchor=(-0.02, -0.02))

    ax.set_xlim(-1.45, 1.45); ax.set_ylim(-1.45, 1.45); ax.set_aspect('equal'); ax.axis('off')
    fig.tight_layout()
    for ext in ('png', 'svg'):
        fig.savefig(os.path.join(outdir, f"calcium_coupling_network.{ext}"),
                    dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"drew {n} candidate transcripts around [Ca2+]cyt "
          f"({(s['period_class']=='invariant').sum()} invariant, "
          f"{(s['period_class']=='shortened').sum()} shortened)")
    print("wrote", os.path.join(outdir, "calcium_coupling_network.png/.svg"))


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0] if len(a) > 0 else "data/gene_calcium_bidirectional.csv",
         a[1] if len(a) > 1 else "output")
