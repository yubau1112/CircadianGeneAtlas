#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 9/11 修改（Checklist 第 29 項）：圖上文字統一為 clock-associated prediction score、移除內部稱呼 breakthrough／clean model、Known rhythmic 改 CircaDB rhythmic；輸出改新檔名，其餘不變
"""
CGA MS3 — fig3 + fig7 突破版（無循環，用 TCGA 特徵 + 新候選）
================================================================
fig3：AI vs 文獻散點，標新的 399 定版候選。
fig7：Top-30 候選 × 19 錨點，改用 TCGA 共表達(radj_，非 GTEx)聚類熱圖。
輸入：step_final_scores.csv、step_final_candidates.csv、coexpr_features_adj.csv
輸出：figures/fig3_ai_vs_literature_v3.png、figures/fig7_candidate_heatmap_v3.png
================================================================
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTDIR="figures"; os.makedirs(OUTDIR,exist_ok=True)
PURPLE,ORANGE,GREY,BLUE="#3b0764","#f59e0b","#9ca3af","#2563eb"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":11,"axes.spines.top":False,
                     "axes.spines.right":False,"figure.dpi":300})

def fig3():
    df=pd.read_csv("step_final_scores.csv")
    cand=set(pd.read_csv("step_final_candidates.csv")["gene_symbol"])
    df["x"]=np.log1p(pd.to_numeric(df["annotation_combined_score"],errors="coerce").fillna(0))
    is_cand=df["gene_symbol"].isin(cand)
    is_circ=df["is_circadb"]==1
    other=~is_cand & ~is_circ
    fig,ax=plt.subplots(figsize=(9,6.5))
    ax.scatter(df.loc[other,"x"],df.loc[other,"ai_score"],s=5,c=GREY,alpha=0.25,edgecolor="none",label="Other genes")
    ax.scatter(df.loc[is_circ,"x"],df.loc[is_circ,"ai_score"],s=13,c=BLUE,alpha=0.5,edgecolor="none",label="CircaDB rhythmic genes")
    ax.scatter(df.loc[is_cand,"x"],df.loc[is_cand,"ai_score"],s=13,c=ORANGE,edgecolor="#7c2d12",lw=0.3,alpha=0.85,
               label="Under-valued candidates (n=%d)"%is_cand.sum(),zorder=3)
    # 不標個別基因名（文獻分數常相同會疊死；散點圖重點是分布 pattern，非名字）。
    # 在橘色帶內加一行說明，指出「高 AI + 低文獻 = 被低估」。
    ax.text(0.15,1.005,"under-valued zone (high prediction score · low literature)",
            fontsize=9,color="#b45309",va="bottom",ha="left")
    ax.axhspan(df[is_cand]["ai_score"].min(),1.02,xmin=0,xmax=1,color=ORANGE,alpha=0.05)
    ax.set_xlabel("Literature score   log(1 + combined score)")
    ax.set_ylabel("Clock-associated prediction score")
    ax.set_title("Prioritization of candidate circadian-associated genes")
    ax.legend(loc="lower right",frameon=False,fontsize=9)
    plt.tight_layout(); plt.savefig(f"{OUTDIR}/fig3_ai_vs_literature_v3.png"); plt.close()
    print("  fig3_ai_vs_literature_v3.png")

def fig7():
    import seaborn as sns
    co=pd.read_csv("coexpr_features_adj.csv")
    cand=pd.read_csv("step_final_candidates.csv")
    tcga=[c for c in co.columns if c.startswith("radj_")]
    top_genes=cand["gene_symbol"].head(30).tolist()
    sub=co[co["gene_symbol"].isin(top_genes)].set_index("gene_symbol").reindex(top_genes)
    M=sub[tcga].copy()
    M.columns=[c.replace("radj_","") for c in tcga]
    M=M.dropna(how="any")
    vmax=np.nanmax(np.abs(M.values))
    g=sns.clustermap(M,cmap="RdBu_r",center=0,vmin=-vmax,vmax=vmax,
                     method="average",metric="correlation",
                     figsize=(11,11),linewidths=0.3,linecolor="white")
    g.ax_heatmap.set_position(      [0.24,0.06,0.55,0.70])
    g.ax_col_dendrogram.set_position([0.24,0.77,0.55,0.14])
    g.ax_row_dendrogram.set_position([0.09,0.06,0.13,0.70])
    g.cax.set_position(             [0.92,0.20,0.02,0.42])
    g.cax.set_ylabel("TCGA co-expression",fontsize=10)
    g.ax_heatmap.set_xlabel("");g.ax_heatmap.set_ylabel("")
    g.ax_heatmap.tick_params(axis="x",labelsize=8);g.ax_heatmap.tick_params(axis="y",labelsize=8)
    g.figure.suptitle("Top candidates × core-clock anchors (TCGA co-expression)",
                      fontsize=12,y=0.955)
    g.savefig(f"{OUTDIR}/fig7_candidate_heatmap_v3.png",dpi=300)
    plt.close()
    print("  fig7_candidate_heatmap_v3.png (TCGA, clustered)")

def main():
    print("輸出到 figures/")
    try: fig3()
    except Exception as e: print("[fig3 跳過]",e)
    try: fig7()
    except Exception as e: print("[fig7 跳過]",e)
    print("完成。")

if __name__=="__main__":
    main()