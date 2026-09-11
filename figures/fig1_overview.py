#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CGA Figure 1（Overview：Stage 1 資料庫建置 + Stage 2 AI 發掘）
來源：7/6–7/7「CircadianGeneAtlas 投影片規劃與策略討論」對話框
  第 285 則初版 make_figure1.py，加上其後的修改（平台名正體、Stage 1 數字分層、Anchors 框加高、
  11,922、分類器與候選框文字、箭頭接框、紅色虛線與箭頭分開、漏斗註記移到底部）
9/11 以本文定稿圖（image1.png）疊圖驗證重建正確後，再做修正。

用法：python make_figure1.py orig   → 重現 7 月定稿（驗證用）
      python make_figure1.py        → 9/11 修正版（Checklist 第 26、27 項）：
          漏斗註記刪除誤植的「11,420 +GTEx」；p≈6e-10 → p = 5.8e-10；ROC-AUC ≈ 0.80 → 0.805 ± 0.018
"""
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.lines import Line2D

ORIG = len(sys.argv) > 1 and sys.argv[1] == "orig"
OUT = "CGA_Figure1_orig_rebuild" if ORIG else "CGA_Figure1_20260911"
plt.rcParams['svg.fonttype'] = 'none'

# ---- palette ----
TEAL   = ("#E1F5EE", "#0F6E56")   # stage 1
BLUE   = ("#E6F1FB", "#185FA5")   # TCGA omics
PURPLE = ("#EEEDFE", "#534AB7")   # anchors
GREEN  = ("#EAF3DE", "#3B6D11")   # stage 2 done
RED    = ("#FCEBEB", "#A32D2D")   # circular
GRAYL  = "#F1EFE8"
GRAYB  = "#B4B2A9"
INK    = "#2C2C2A"

fig = plt.figure(figsize=(8.5, 10.6), dpi=300)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 680); ax.set_ylim(0, 838)
ax.invert_yaxis(); ax.axis("off")
fig.patch.set_facecolor("white")


def box(x, y, w, h, fill, edge, title, subs, tsize=10.5, ssize=8.6, align="center"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=8",
                                linewidth=0.8, edgecolor=edge, facecolor=fill, mutation_aspect=1))
    cx = x + w / 2
    if align == "center":
        ax.text(cx, y + 16, title, ha="center", va="center", fontsize=tsize, color=edge, weight="medium")
        for i, s in enumerate(subs):
            ax.text(cx, y + 34 + i * 15, s, ha="center", va="center", fontsize=ssize, color=edge)
    else:
        ax.text(x + 16, y + 16, title, ha="left", va="center", fontsize=tsize, color=edge, weight="medium")
        for i, s in enumerate(subs):
            ax.text(x + 16, y + 34 + i * 15, s, ha="left", va="center", fontsize=ssize, color=edge)


def panel(x, y, w, h, label, right):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=14",
                                linewidth=0.8, edgecolor=GRAYB, facecolor=GRAYL, alpha=0.55))
    ax.text(x + 22, y + 22, label, ha="left", va="center", fontsize=12.5, color=INK, weight="medium")
    ax.text(x + w - 20, y + 22, right, ha="right", va="center", fontsize=8.6, color="#5F5E5A")


def arrow(x1, y1, x2, y2, color=INK, lw=1.3):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=11,
                                 linewidth=lw, color=color, shrinkA=0, shrinkB=0))


def dashed_arrow(x1, y1, x2, y2, color=INK, lw=1.3):
    """虛線與實心箭頭分開畫，避免虛線疊在箭頭上"""
    dx, dy = x2 - x1, y2 - y1; L = (dx * dx + dy * dy) ** 0.5; ux, uy = dx / L, dy / L
    xm, ym = x2 - ux * 9, y2 - uy * 9
    ax.plot([x1, xm], [y1, ym], linestyle=(0, (4, 3)), color=color, linewidth=lw, solid_capstyle='butt', zorder=1)
    ax.add_patch(FancyArrowPatch((xm, ym), (x2, y2), arrowstyle="-|>", mutation_scale=11,
                                 linewidth=lw, color=color, shrinkA=0, shrinkB=0, zorder=2))


def elbow(pts, color=INK, lw=1.3):
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
    ax.add_line(Line2D(xs, ys, color=color, linewidth=lw, solid_capstyle="round"))
    arrow(pts[-2][0], pts[-2][1], pts[-1][0], pts[-1][1], color=color, lw=lw)


P_TXT = "clock vs random p\u22486e-10" if ORIG else "HIGH vs VERY_LOW literature p = 5.8e-10"
AUC_TXT = "coexpr+meth \u00b7 ROC-AUC \u2248 0.80" if ORIG else "coexpr+meth \u00b7 ROC-AUC = 0.805 \u00b1 0.018"
FUNNEL = ("gene funnel:  12,617 lit \u2192 11,922 TCGA \u2192 11,420 +GTEx \u2192 10,370 +methylation \u2192 399 candidates"
          if ORIG else "gene funnel:  12,617 lit \u2192 11,922 TCGA \u2192 10,370 +methylation \u2192 399 candidates")

# ===== Panel A =====
panel(30, 30, 620, 378, "A \u00b7 Database construction (stage 1)", "Circadian Gene Atlas")
box(60, 80, 240, 58, *TEAL, "PubMed literature search", ["6 circadian terms \u00b7 149,754 articles"])
arrow(180, 138, 180, 168)
box(60, 168, 240, 58, *TEAL, "PubTator3 gene extraction", ["12,629 genes \u2192 12,617 in circadian table"])
arrow(180, 226, 180, 256)
box(60, 256, 240, 70, *TEAL, "Dual literature scoring",
    ["context + NCBI validation", "HIGH 377 \u00b7 MED 922 \u00b7 LOW \u00b7 VERY_LOW"])
box(360, 168, 260, 158, *BLUE, "TCGA multi-omics (24 cancers)", [], align="left")
for yy, t in [(214, "expression \u00b7 DEG \u00b7 methylation"), (234, "miRNA \u00b7 CNV \u00b7 mutation"), (254, "survival \u00b7 COSMIC"),
              (284, "Web platform: browse, search,"), (302, "per-gene views, downloads")]:
    ax.text(384, yy, t, ha="left", va="center", fontsize=8.6, color=BLUE[1])
arrow(300, 285, 360, 285)
arrow(180, 326, 180, 362)
ax.text(196, 345, "12,617 genes + literature scores", ha="left", va="center", fontsize=8.6, color=INK)

# ===== Panel B =====
panel(30, 430, 620, 378, "B \u00b7 AI-enhanced discovery (stage 2)", "under-valued clock genes")
box(60, 480, 240, 70, *PURPLE, "Anchors & gene universe",
    ["19 anchors + 11,903 candidates", "12,617 \u2192 11,922 with TCGA expr."])
arrow(180, 550, 180, 568)
box(60, 568, 240, 70, *GREEN, "Co-expression + correction", ["TCGA 24 cancers \u00b7 Fisher-z", P_TXT])
arrow(180, 638, 180, 668)
box(60, 668, 240, 70, *GREEN, "Supervised classifier", ["CircaDB high-amplitude positives", AUC_TXT])
box(360, 602, 260, 72, *GREEN, "Under-valued candidates",
    ["high AI \u00b7 low literature \u00b7 final 399", "IEG \u00b7 CCN \u00b7 RGS families, THBS1/MCL1"])
elbow([(300, 703), (344, 703), (344, 638), (360, 638)])
box(360, 700, 260, 44, *RED, "GTEx branch = circular", ["0.871 not used (feature = label src)"], ssize=8.6)
dashed_arrow(360, 722, 300, 722, color=RED[1])

# ===== legend =====
def sq(x, c, txt):
    ax.add_patch(Rectangle((x, 762), 12, 12, facecolor=c, edgecolor="none"))
    ax.text(x + 18, 768, txt, ha="left", va="center", fontsize=8.4, color=INK)


sq(60, TEAL[1], "Stage 1"); sq(152, GREEN[1], "Stage 2 (done)"); sq(286, RED[1], "Circular (not used)"); sq(446, BLUE[1], "TCGA multi-omics")
ax.text(60, 792, FUNNEL, fontsize=7.5, color="#5F5E5A", ha="left", va="center")

fig.savefig(f"{OUT}.svg", facecolor="white", bbox_inches="tight", pad_inches=0.15)
fig.savefig(f"{OUT}.pdf", facecolor="white", bbox_inches="tight", pad_inches=0.15)
fig.savefig(f"{OUT}.png", facecolor="white", dpi=300, bbox_inches="tight", pad_inches=0.15)
print(f"完成：{OUT}.svg／.pdf／.png")
