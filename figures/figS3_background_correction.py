#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CGA MS3 — 論文圖產生器（一次產出全部 figure）
================================================================
讀現有 CSV，輸出 6 張 300dpi PNG 到 ./figures/。每張都有明確 message：
  fig1_roc                ROC：TCGA(0.756) vs GTEx(0.853) vs 合併(0.871)
  fig2_feature_importance 合併模型特徵重要性（前 12 全是 GTEx）
  fig3_ai_vs_literature   AI 分數 vs 文獻分數散點（右上=被低估候選）
  fig4_score_dist         已知節律 / 候選 / VERY_LOW 的 ai_score 分布
  fig5_background_correct 背景校正前後（HIGH vs VERY_LOW 由分不開→分得開）
  fig6_strategy_auc       各策略 AUC 比較（GTEx 換訊源勝出）

輸入（同資料夾）：coexpr_features_adj.csv、gtex_features_adj.csv、
  coexpr_features_by_gene.csv、CircaDB_HA_query.csv、
  step2_final_scores_all_genes.csv、step3_candidates.csv
輸出：figures/*.png
================================================================
需要：matplotlib（pip install matplotlib）
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_curve, roc_auc_score

OUTDIR = "figures"
RS, FDR, LVL = 42, 0.05, "annotation_confidence_level"
PURPLE, ORANGE, GREY, BLUE = "#3b0764", "#f59e0b", "#9ca3af", "#2563eb"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 11, "axes.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 300,
})
os.makedirs(OUTDIR, exist_ok=True)


def load():
    d = {}
    d["adj"]  = pd.read_csv("coexpr_features_adj.csv")
    d["gtex"] = pd.read_csv("gtex_features_adj.csv")
    d["raw"]  = pd.read_csv("coexpr_features_by_gene.csv")
    cd = pd.read_csv("CircaDB_HA_query.csv")
    d["circ"] = set(cd[cd["FDR"] <= FDR]["Gene_Symbol"].astype(str).str.strip())
    d["scores"] = pd.read_csv("step2_final_scores_all_genes.csv")
    d["cand"]   = pd.read_csv("step3_candidates.csv")
    return d


def build_train(adj, gtex, circ):
    tcga = [c for c in adj.columns if c.startswith("radj_")]
    gt   = [c for c in gtex.columns if c.startswith("gtex_radj_")]
    df = adj.merge(gtex[["gene_symbol"]+gt], on="gene_symbol", how="inner")
    in_circ = df["gene_symbol"].isin(circ)
    is_anc = df["is_anchor"] == 1
    pos = df[in_circ & ~is_anc]
    negp = df[(df[LVL]=="VERY_LOW") & ~is_anc & ~in_circ]
    neg = negp.sample(n=len(pos), random_state=RS)
    tr = pd.concat([pos, neg], ignore_index=True)
    y = np.r_[np.ones(len(pos),int), np.zeros(len(neg),int)]
    return tr, y, tcga, gt


def oof(tr, y, cols, model):
    pipe = Pipeline([("i",SimpleImputer(strategy="median")),("s",StandardScaler()),("m",model)])
    cv = StratifiedKFold(5, shuffle=True, random_state=RS)
    return cross_val_predict(pipe, tr[cols].to_numpy(float), y, cv=cv, method="predict_proba", n_jobs=-1)[:,1]


def fig1_roc(tr, y, tcga, gt):
    rf = lambda: RandomForestClassifier(n_estimators=500, random_state=RS, class_weight="balanced", n_jobs=-1)
    sets = [("TCGA tumour", tcga, GREY), ("GTEx normal", gt, ORANGE), ("Combined", tcga+gt, PURPLE)]
    plt.figure(figsize=(5,5))
    for name, cols, col in sets:
        p = oof(tr, y, cols, rf())
        fpr, tpr, _ = roc_curve(y, p); a = roc_auc_score(y, p)
        plt.plot(fpr, tpr, color=col, lw=2, label=f"{name} (AUC={a:.3f})")
    plt.plot([0,1],[0,1],"--",color="#cccccc",lw=1)
    plt.xlabel("False positive rate"); plt.ylabel("True positive rate")
    plt.title("Classifier ROC by expression source")
    plt.legend(loc="lower right", frameon=False, fontsize=9)
    plt.tight_layout(); plt.savefig(f"{OUTDIR}/fig1_roc.png"); plt.close()
    print("  fig1_roc.png")


def fig2_importance(tr, y, tcga, gt):
    cols = tcga+gt
    pipe = Pipeline([("i",SimpleImputer(strategy="median")),("s",StandardScaler()),
                     ("m",RandomForestClassifier(n_estimators=500,random_state=RS,class_weight="balanced",n_jobs=-1))])
    pipe.fit(tr[cols].to_numpy(float), y)
    # 取前 12，最重要在最上面
    imp = pd.Series(pipe.named_steps["m"].feature_importances_, index=cols).sort_values(ascending=False).head(15)
    imp = imp.iloc[::-1]   # 反轉：barh 由下往上畫，最大值落在頂端
    colors = [ORANGE if c.startswith("gtex_") else PURPLE for c in imp.index]
    labels = [c.replace("gtex_radj_","GTEx: ").replace("radj_","TCGA: ") for c in imp.index]
    plt.figure(figsize=(6.2,6.4))
    bars = plt.barh(range(len(imp)), imp.values, color=colors, height=0.72)
    plt.yticks(range(len(imp)), labels, fontsize=9)
    plt.xlabel("Feature importance (Gini)")
    plt.title("Top features are GTEx normal-tissue co-expression", fontsize=12)
    plt.xlim(0, imp.max()*1.15)
    for i,v in enumerate(imp.values):
        plt.text(v+imp.max()*0.01, i, f"{v:.3f}", va="center", fontsize=8, color="#555555")
    from matplotlib.patches import Patch
    plt.legend(handles=[Patch(color=ORANGE,label="GTEx normal"),Patch(color=PURPLE,label="TCGA tumour")],
               frameon=False, fontsize=9, loc="lower right")
    plt.tight_layout(); plt.savefig(f"{OUTDIR}/fig2_feature_importance.png"); plt.close()
    print("  fig2_feature_importance.png")


def fig3_scatter(scores, circ, cand):
    rng = np.random.default_rng(42)
    df = scores.copy()
    df["lit"] = pd.to_numeric(df["annotation_combined_score"], errors="coerce").fillna(0)
    df["x"] = np.log1p(df["lit"]) + rng.normal(0, 0.06, len(df))   # 水平 jitter 消柵欄
    is_circ = df["gene_symbol"].isin(circ)
    is_cand = df["gene_symbol"].isin(set(cand["gene_symbol"]))
    is_other = ~is_circ & ~is_cand

    fig, ax = plt.subplots(figsize=(7,5.2))
    # 被低估區塊底色（高 AI、低文獻）
    y_cut = df.loc[is_cand, "ai_score"].min()
    x_cut = np.log1p(20)   # LOW 上界 combined_score<20
    ax.axhspan(y_cut, 1.02, xmin=0, xmax=(x_cut+0.3)/ (df["x"].max()+0.3),
               color=ORANGE, alpha=0.06, zorder=0)
    ax.text(0.05, 1.005, "under-valued zone\n(high AI · low literature)",
            fontsize=8, color="#b45309", va="top")

    # 疊圖順序：其他→已知→候選(最上層、最醒目)
    ax.scatter(df.loc[is_other,"x"], df.loc[is_other,"ai_score"],
               s=5, c=GREY, alpha=0.18, linewidths=0, label="Other genes", zorder=1)
    ax.scatter(df.loc[is_circ,"x"], df.loc[is_circ,"ai_score"],
               s=12, c=BLUE, alpha=0.45, linewidths=0, label="Known rhythmic (CircaDB)", zorder=2)
    ax.scatter(df.loc[is_cand,"x"], df.loc[is_cand,"ai_score"],
               s=13, c=ORANGE, edgecolor="#7c2d12", lw=0.3, alpha=0.85,
               label="Under-valued candidates (n=%d)"%is_cand.sum(), zorder=3)

    # 標幾個代表基因（帶細箭頭，避免疊字）
    for g,(dx,dy) in {"MCL1":(0.25,0.02),"THBS1":(0.25,0.02),"STC1":(0.25,-0.02),
                      "MBD5":(0.25,0.02),"OPHN1":(-0.55,0.02)}.items():
        r = df[df["gene_symbol"]==g]
        if len(r):
            x0,y0 = r["x"].iloc[0], r["ai_score"].iloc[0]
            ax.annotate(g, (x0,y0), xytext=(x0+dx, min(y0+dy,1.0)), fontsize=8,
                        arrowprops=dict(arrowstyle="-", lw=0.5, color="#7c2d12"))

    ax.set_xlabel("Literature score   log(1 + combined score)")
    ax.set_ylabel("AI rhythmicity score")
    ax.set_title("High AI + low literature = under-valued candidates", fontsize=12)
    ax.set_ylim(-0.02, 1.06)
    ax.legend(loc="lower right", frameon=False, fontsize=8)
    plt.tight_layout(); plt.savefig(f"{OUTDIR}/fig3_ai_vs_literature.png"); plt.close()
    print("  fig3_ai_vs_literature.png")


def fig4_dist(scores, circ, cand):
    df = scores.copy()
    known = df[df["gene_symbol"].isin(circ)]["ai_score"].dropna()
    cands = df[df["gene_symbol"].isin(set(cand["gene_symbol"]))]["ai_score"].dropna()
    vlow  = df[(df[LVL]=="VERY_LOW") & (~df["gene_symbol"].isin(circ))]["ai_score"].dropna()
    data = [known, cands, vlow]
    labels = [f"Known rhythmic\n(n={len(known)})", f"Candidates\n(n={len(cands)})", f"VERY_LOW\n(n={len(vlow)})"]
    plt.figure(figsize=(5.5,5))
    parts = plt.violinplot(data, showmedians=True, widths=0.8)
    for pc, c in zip(parts["bodies"], [BLUE, ORANGE, GREY]):
        pc.set_facecolor(c); pc.set_alpha(0.6)
    plt.xticks([1,2,3], labels, fontsize=9)
    plt.ylabel("AI rhythmicity score")
    plt.title("Candidates track known rhythmic genes")
    plt.tight_layout(); plt.savefig(f"{OUTDIR}/fig4_score_dist.png"); plt.close()
    print("  fig4_score_dist.png")


def _cos_to_template(feat_df, cols, lvlcol):
    hi = feat_df[feat_df[lvlcol]=="HIGH"]
    tmpl = np.nanmean(hi[cols].to_numpy(float), axis=0)
    def cos(x):
        m = ~np.isnan(x) & ~np.isnan(tmpl)
        if m.sum()<5: return np.nan
        a,b = x[m], tmpl[m]; return a@b/(np.linalg.norm(a)*np.linalg.norm(b)+1e-9)
    return np.array([cos(r) for r in feat_df[cols].to_numpy(float)])


def fig5_bgcorrect(raw, adj):
    raw_rho = [c for c in raw.columns if c.startswith("rho_")]
    adj_rho = [c for c in adj.columns if c.startswith("radj_")]
    raw = raw[~raw["gene_symbol"].isin(
        {"CLOCK","ARNTL","ARNTL2","PER1","PER2","PER3","CRY1","CRY2","NR1D1","NR1D2",
         "RORA","RORB","RORC","NPAS2","CSNK1D","CSNK1E","TIMELESS","BHLHE40","BHLHE41"})].copy()
    raw["sim"] = _cos_to_template(raw, raw_rho, LVL)
    adj2 = adj[adj["is_anchor"]==0].copy()
    adj2["sim"] = _cos_to_template(adj2, adj_rho, LVL)

    fig, axes = plt.subplots(1,2, figsize=(9,4.2), sharey=False)
    for ax, d, title in [(axes[0],raw,"Before correction"),(axes[1],adj2,"After correction")]:
        h = d[d[LVL]=="HIGH"]["sim"].dropna(); l = d[d[LVL]=="VERY_LOW"]["sim"].dropna()
        try: _,p = mannwhitneyu(h,l,alternative="greater")
        except: p=np.nan
        ax.hist(l, bins=40, color=GREY, alpha=0.6, density=True, label="VERY_LOW")
        ax.hist(h, bins=40, color=PURPLE, alpha=0.6, density=True, label="HIGH")
        ax.set_title(f"{title}\n(HIGH vs VERY_LOW, p={p:.1e})", fontsize=10)
        ax.set_xlabel("Cosine similarity to clock template"); ax.legend(frameon=False, fontsize=8)
    axes[0].set_ylabel("Density")
    fig.suptitle("Background correction reveals HIGH vs VERY_LOW separation", fontsize=12)
    plt.tight_layout(); plt.savefig(f"{OUTDIR}/fig5_background_correct.png"); plt.close()
    print("  fig5_background_correct.png")


def fig6_strategy():
    # 本 session 實測 AUC（各 step 報告）
    strat = [
        ("TCGA tumour", 0.756, GREY),
        ("+ derived feats", 0.762, GREY),
        ("STRING network", 0.527, GREY),
        ("per-cancer max", 0.712, GREY),
        ("E-box motif", 0.539, GREY),
        ("GTEx normal", 0.853, ORANGE),
        ("Combined (final)", 0.871, PURPLE),
    ]
    names = [s[0] for s in strat]; vals=[s[1] for s in strat]; cols=[s[2] for s in strat]
    plt.figure(figsize=(7,4.5))
    plt.bar(range(len(vals)), vals, color=cols)
    plt.axhline(0.5, ls="--", color="#cccccc", lw=1)
    plt.xticks(range(len(names)), names, rotation=35, ha="right", fontsize=9)
    plt.ylabel("CV ROC-AUC"); plt.ylim(0.4,0.95)
    for i,v in enumerate(vals):
        plt.text(i, v+0.008, f"{v:.3f}", ha="center", fontsize=8)
    plt.title("Strategy comparison (positives = CircaDB, no circularity)")
    plt.tight_layout(); plt.savefig(f"{OUTDIR}/fig6_strategy_auc.png"); plt.close()
    print("  fig6_strategy_auc.png")


def _run():
    print("讀資料..."); d = load()
    tr, y, tcga, gt = build_train(d["adj"], d["gtex"], d["circ"])
    print(f"訓練集={len(tr)}；輸出到 ./{OUTDIR}/")
    steps = [("fig1", lambda: fig1_roc(tr,y,tcga,gt)),
             ("fig2", lambda: fig2_importance(tr,y,tcga,gt)),
             ("fig3", lambda: fig3_scatter(d["scores"], d["circ"], d["cand"])),
             ("fig4", lambda: fig4_dist(d["scores"], d["circ"], d["cand"])),
             ("fig5", lambda: fig5_bgcorrect(d["raw"], d["adj"])),
             ("fig6", lambda: fig6_strategy())]
    for name, fn in steps:
        try: fn()
        except Exception as e: print(f"  [{name} 跳過] {e}")
    print("完成，圖在 figures/")

if __name__ == "__main__":
    _run()