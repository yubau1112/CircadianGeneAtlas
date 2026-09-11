#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 9/11 修改（Checklist 第 29 項）：圖上文字統一為 clock-associated prediction score、移除內部稱呼 breakthrough／clean model、Known rhythmic 改 CircaDB rhythmic；輸出改新檔名，其餘不變
"""
CGA MS3 — fig4 + fig10 突破版重繪（無循環，定版 399 候選）
================================================================
取代含 GTEx 循環的舊版 fig4/fig10（那版是 217 候選）。
用定版模型分數 step_final_scores.csv（11,922 全基因 ai_score）重畫。

fig4_score_dist_v3：候選 vs CircaDB已知 vs 其他 的 AI 分數分布（violin/box）
fig10_threshold_sens_v3：AI 分數門檻 vs 候選數（95% 百分位那點 = 399）

輸入：step_final_scores.csv、step_final_candidates.csv
輸出：figures/fig4_score_dist_v3.png、fig10_threshold_sens_v3.png
================================================================
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTDIR="figures"; os.makedirs(OUTDIR,exist_ok=True)
PURPLE,ORANGE,GREY,BLUE="#3b0764","#f59e0b","#9ca3af","#2563eb"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":11,"axes.spines.top":False,
                     "axes.spines.right":False,"figure.dpi":300})
SCORES="step_final_scores.csv"; CAND="step_final_candidates.csv"
LVL="annotation_confidence_level"; LOW={"VERY_LOW","LOW"}

def fig4():
    df=pd.read_csv(SCORES)
    cand=set(pd.read_csv(CAND)["gene_symbol"])
    # 三組：定版候選 / CircaDB已知節律 / 其他
    is_cand=df["gene_symbol"].isin(cand)
    is_circ=df["is_circadb"]==1 if "is_circadb" in df.columns else pd.Series(False,index=df.index)
    grp_cand=df.loc[is_cand,"ai_score"].dropna().values
    grp_circ=df.loc[is_circ,"ai_score"].dropna().values
    grp_other=df.loc[~is_cand & ~is_circ,"ai_score"].dropna().values

    fig,ax=plt.subplots(figsize=(7,6))
    data=[grp_other,grp_circ,grp_cand]
    labels=[f"Other genes\n(n={len(grp_other)})",
            f"CircaDB rhythmic\n(n={len(grp_circ)})",
            f"Candidates\n(n={len(grp_cand)})"]
    parts=ax.violinplot(data,showmedians=True,widths=0.8)
    cols=[GREY,BLUE,ORANGE]
    for pc,c in zip(parts["bodies"],cols):
        pc.set_facecolor(c); pc.set_alpha(0.6); pc.set_edgecolor(c)
    for key in ["cmedians","cmaxes","cmins","cbars"]:
        if key in parts: parts[key].set_color("#444444")
    ax.set_xticks([1,2,3]); ax.set_xticklabels(labels)
    ax.set_ylabel("Clock-associated prediction score")
    ax.set_title("Candidates score close to CircaDB rhythmic genes")
    # 中位數標註
    for i,d in enumerate(data,1):
        med=np.median(d)
        ax.text(i+0.34,med,f"{med:.2f}",va="center",ha="left",fontsize=9,color="#333333",
                bbox=dict(boxstyle="round,pad=0.15",fc="white",ec="none",alpha=0.85))
    plt.tight_layout(); plt.savefig(f"{OUTDIR}/fig4_score_dist_v3.png"); plt.close()
    print(f"  fig4_score_dist_v3.png (candidates n={len(grp_cand)})")

def fig10():
    df=pd.read_csv(SCORES)
    # 只在「可成為候選」的池裡算(非anchor、非CircaDB、低文獻)——與定版候選定義一致
    ia=df["is_anchor"]==1 if "is_anchor" in df.columns else pd.Series(False,index=df.index)
    circ=df["is_circadb"]==1 if "is_circadb" in df.columns else pd.Series(False,index=df.index)
    pool=df[~ia & ~circ & df[LVL].isin(LOW)].copy()

    # 掃百分位門檻(對全體 ai_score)，數符合的候選數
    pcts=np.arange(90,99.6,0.5)
    counts=[]
    for p in pcts:
        thr=df["ai_score"].quantile(p/100)
        counts.append((pool["ai_score"]>=thr).sum())

    fig,ax=plt.subplots(figsize=(7.5,5.5))
    ax.plot(pcts,counts,"o-",color=PURPLE,mfc=ORANGE,mec=PURPLE,mew=1.2,ms=6,lw=2)
    # 標出 95% 那點(定版用的門檻)
    thr95=df["ai_score"].quantile(0.95)
    n95=(pool["ai_score"]>=thr95).sum()
    ax.scatter([95],[n95],s=140,facecolor=ORANGE,edgecolor=PURPLE,lw=2,zorder=5)
    # 標註放左上空白區(資料在右下方遞增),用箭頭指向 95% 那點,避免壓線
    # 文字放 399 點的左下方(近、空白)，直線箭頭指向該點
    # 文字放 399 點右下方近處,直線箭頭
    # 文字放 399 點正下方,垂直箭頭往上指(x固定=95,不與斜線平行糾纏)
    ax.annotate("top 5% = 399", xy=(95, n95),
                xytext=(95, n95-max(counts)*0.30),
                fontsize=10.5, color="#7c2d12", ha="center", va="top", fontweight="bold",
                zorder=20,
                arrowprops=dict(arrowstyle="->", color="#7c2d12", lw=1.8, shrinkA=0, shrinkB=8))
    ax.set_xlabel("Prediction-score percentile threshold")
    ax.set_ylabel("Number of candidates\n(non-anchor, non-CircaDB, low-literature)")
    ax.set_title("Candidate count vs prediction-score threshold")
    ax.invert_xaxis()  # 門檻越嚴(右→左越寬鬆)，讓 95% 在右、更嚴在左
    plt.tight_layout(); plt.savefig(f"{OUTDIR}/fig10_threshold_sens_v3.png"); plt.close()
    print(f"  fig10_threshold_sens_v3.png (top 5% = {n95})")

def main():
    print("重繪 fig4/fig10 突破版 → figures/")
    try: fig4()
    except Exception as e: print("[fig4 跳過]",e)
    try: fig10()
    except Exception as e: print("[fig10 跳過]",e)
    print("完成。")

if __name__=="__main__":
    main()