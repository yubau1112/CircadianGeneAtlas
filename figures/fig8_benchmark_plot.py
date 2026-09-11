#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure 8（兩面板）繪圖：只讀 fig8_roc_points.tsv、fig8_auc.tsv，圖上數字皆由檔案讀入"""
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.spines.top": False,
                     "axes.spines.right": False, "svg.fonttype": "none"})
COLORS = {"CGDB": "#3b0764", "GO": "#f59e0b", "KEGG": "#2563eb", "Reactome": "#059669"}
ORDER = ["CGDB", "GO", "KEGG", "Reactome"]
TITLES = {"A": ("A", "Literature-derived circadian score", "anchor genes excluded"),
          "B": ("B", "Machine-learning prediction score", "anchor and training genes excluded")}
pts = pd.read_csv("fig8_roc_points.tsv", sep="\t")
auc = pd.read_csv("fig8_auc.tsv", sep="\t")

fig, axes = plt.subplots(1, 2, figsize=(12.5, 6.2), dpi=300)
for ax, panel in zip(axes, ["A", "B"]):
    for g in ORDER:
        p = pts[(pts.panel == panel) & (pts.gene_set == g)]
        a = auc[(auc.panel == panel) & (auc.gene_set == g)].iloc[0]
        ax.plot(p.fpr, p.tpr, lw=2.0, color=COLORS[g], label=f"{g} (AUC = {a.AUC:.3f}, n = {int(a.n_pos)})")
    ax.plot([0, 1], [0, 1], "--", color="#999999", lw=1)
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("False positive rate"); ax.set_ylabel("True positive rate")
    tag, title, sub = TITLES[panel]
    ax.set_title(f"{title}\n({sub})", fontsize=11.5)
    ax.text(-0.14, 1.06, tag, transform=ax.transAxes, fontsize=16, fontweight="bold")
    ax.legend(loc="lower right", frameon=True, fontsize=9)
plt.tight_layout()
plt.savefig("CGA_Figure8_benchmark_2panel.png", dpi=300, facecolor="white")
plt.savefig("CGA_Figure8_benchmark_2panel.svg", facecolor="white")
print("已輸出 CGA_Figure8_benchmark_2panel.png、.svg")
