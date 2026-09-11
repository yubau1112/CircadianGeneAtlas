#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CGA Figure 2（Machine-learning workflow）
來源：7/9–7/10「CGA 手稿投稿與進度交接」對話框第 110 則的最終程式（育葆第 111 則確認「圖二可以了」）
9/11 修改（二）：背景校正的比較對象更正為 HIGH vs VERY_LOW literature（原寫 clock vs random，與 Figure S3 及原始分析不符）
9/11 修改（Checklist 第 27 項，只改兩個字串，版面不動）：
  clock vs random p≈6e-10          → p = 5.8e-10（與正文 P = 5.8 × 10⁻¹⁰ 一致）
  ROC-AUC ≈ 0.80                   → ROC-AUC = 0.805 ± 0.018（10 個種子的平均）
  新字串以 get_window_extent 實測：80.7–319.3，在框內可用範圍 73–327 之內
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
plt.rcParams['svg.fonttype'] = 'none'

GRAYL = ("#F1EFE8", "#5F5E5A"); GREEN = ("#EAF3DE", "#3B6D11")
TEAL = ("#E1F5EE", "#0F6E56");  RED = ("#FCEBEB", "#A32D2D"); INK = "#2C2C2A"

fig = plt.figure(figsize=(8.5, 9.0), dpi=300)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 669.6); ax.set_ylim(0, 691.2)
ax.invert_yaxis(); ax.axis("off"); fig.patch.set_facecolor("white")


def box(x, y, w, h, pair, title, subs, tsize=10.5, ssize=8.6):
    fill, edge = pair
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=8",
                                linewidth=0.8, edgecolor=edge, facecolor=fill))
    cx = x + w / 2
    ax.text(cx, y + 18, title, ha="center", va="center", fontsize=tsize, color=edge, weight="medium")
    for i, s in enumerate(subs):
        ax.text(cx, y + 36 + i * 15, s, ha="center", va="center", fontsize=ssize, color=edge)


def arrow(x1, y1, x2, y2, color=INK, dashed=False, lw=1.3):
    if dashed:
        dx, dy = x2 - x1, y2 - y1; L = (dx * dx + dy * dy) ** 0.5; ux, uy = dx / L, dy / L
        xm, ym = x2 - ux * 9, y2 - uy * 9
        ax.plot([x1, xm], [y1, ym], linestyle=(0, (4, 3)), color=color, linewidth=lw, solid_capstyle='butt', zorder=1)
        ax.add_patch(FancyArrowPatch((xm, ym), (x2, y2), arrowstyle="-|>", mutation_scale=11,
                                     linewidth=lw, color=color, linestyle='-', shrinkA=0, shrinkB=0, zorder=2))
    else:
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=11,
                                     linewidth=lw, color=color, shrinkA=0, shrinkB=0))


# 主流程（x 中心 = 200，所有直箭頭頭尾同 x）
box(65, 47, 270, 51, GRAYL, "Data sources", ["TCGA \u00b7 Circadian DB \u00b7 CircaDB"])
arrow(200, 98, 200, 123)
box(65, 123, 270, 65, GRAYL, "Anchors & gene universe",
    ["11,922 genes (19 anchors + 11,903 cand.)", "12,617 lit. \u2192 11,922 with TCGA expr."])
arrow(200, 188, 200, 213)
box(65, 213, 270, 51, GREEN, "Co-expression feature construction", ["TCGA 24 cancers \u00b7 Fisher-z"])
arrow(200, 264, 200, 296)
box(65, 296, 270, 51, GREEN, "Background correction", ["HIGH vs VERY_LOW literature p = 5.8e-10"])
arrow(200, 347, 200, 379)
box(65, 379, 270, 76, GREEN, "Supervised classification",
    ["CircaDB high-amplitude positives", "rAMP top 40% + meth \u00b7 ROC-AUC = 0.805 \u00b1 0.018"])
arrow(200, 455, 200, 484)
box(65, 484, 270, 51, GREEN, "Candidate prioritization", ["final 399 (high AI \u00b7 low lit)"])
arrow(200, 535, 200, 567)
box(65, 567, 270, 76, TEAL, "Biological validation",
    ["anchor sanity \u00b7 TIPARP, IEGs", "candidate\u00d7anchor heatmap (phase)"])

# 側支（水平箭頭：頭尾同 y）
box(371, 366, 270, 49, GREEN, "Methylation branch (done, used)", ["promoter meth \u2192 meth_radj_"], ssize=8.4)
arrow(371, 390.5, 335, 390.5, color=GREEN[1])
box(371, 419, 270, 49, RED, "GTEx-derived features", ["(excluded: information leakage)"], ssize=8.4)
arrow(371, 443.5, 335, 443.5, color=RED[1], dashed=True)
box(371, 482, 270, 64, GRAYL, "Alternative features evaluated",
    ["STRING \u00b7 motif \u00b7 per-cancer \u00b7 adj-normal", "GO 0.61 \u00b7 CGDB 0.65 positives"], ssize=8.4)


def sq(x, c, t):
    ax.add_patch(Rectangle((x, 663), 12, 12, facecolor=c, edgecolor="none"))
    ax.text(x + 18, 669, t, ha="left", va="center", fontsize=8.4, color=INK)


sq(65, GREEN[1], "Done / used"); sq(255, RED[1], "Circular (not used)"); sq(465, GRAYL[1], "Input / tested-not-used")

fig.savefig("CGA_Figure2_workflow_20260911.svg", facecolor="white", bbox_inches="tight", pad_inches=0.15)
fig.savefig("CGA_Figure2_workflow_20260911.png", facecolor="white", dpi=300, bbox_inches="tight", pad_inches=0.15)
print("完成：CGA_Figure2_workflow_20260911.svg、.png")
