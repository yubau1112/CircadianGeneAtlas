#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 9/11 修改（Checklist 第 29 項）：圖上文字統一為 clock-associated prediction score、移除內部稱呼 breakthrough／clean model、Known rhythmic 改 CircaDB rhythmic；輸出改新檔名，其餘不變
"""
CGA MS3 — fig1/2/6 突破版重繪（無循環，定版模型 AUC≈0.805）
================================================================
取代含 GTEx 循環的舊版 fig1/2/6。
定版模型：CircaDB 高幅度正例(rAMP前40%) + coexpr+deriv+meth + L1-Logistic。

fig1_roc_v3：定版模型 ROC 曲線（out-of-fold）
fig2_feature_importance_v3：L1-Logistic 係數絕對值（特徵重要性）
fig6_strategy_auc_v3：各策略 AUC 對照（誠實版，標無循環/循環）

輸入：coexpr_features_adj.csv、meth_features_adj.csv、CircaDB_HA_query.csv
輸出：figures/fig1_roc_v3.png、fig2_feature_importance_v3.png、fig6_strategy_auc_v3.png
================================================================
"""
import os, numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, roc_curve

OUTDIR="figures"; os.makedirs(OUTDIR,exist_ok=True)
PURPLE,ORANGE,GREY,RED="#3b0764","#f59e0b","#9ca3af","#dc2626"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":11,"axes.spines.top":False,
                     "axes.spines.right":False,"figure.dpi":300})
COEXPR="coexpr_features_adj.csv"; METH="meth_features_adj.csv"; CIRCADB="CircaDB_HA_query.csv"
LVL,RS,FDR="annotation_confidence_level",42,0.05; AMP_FRAC=0.40

def build():
    co=pd.read_csv(COEXPR); me=pd.read_csv(METH)
    coex=[c for c in co.columns if c.startswith("radj_")]
    meth=[c for c in me.columns if c.startswith("meth_radj_")]
    df=co.merge(me[["gene_symbol"]+meth],on="gene_symbol",how="left")
    R=df[coex].to_numpy(float)
    POS=["radj_ARNTL","radj_CLOCK","radj_NPAS2","radj_RORA","radj_RORB","radj_RORC","radj_ARNTL2"]
    NEG=["radj_PER1","radj_PER2","radj_PER3","radj_CRY1","radj_CRY2","radj_NR1D1","radj_NR1D2","radj_BHLHE40","radj_BHLHE41"]
    df["d_posneg"]=np.nanmean(df[[c for c in POS if c in coex]],axis=1)-np.nanmean(df[[c for c in NEG if c in coex]],axis=1)
    df["d_mean"]=np.nanmean(R,axis=1); df["d_max"]=np.nanmax(R,axis=1); df["d_sign"]=np.nanmean(np.sign(R),axis=1)
    cols=coex+["d_posneg","d_mean","d_max","d_sign"]+meth
    return df,cols

def get_train(df,cols):
    cd=pd.read_csv(CIRCADB); cd["Gene_Symbol"]=cd["Gene_Symbol"].astype(str).str.strip()
    circ=set(cd[cd["FDR"]<=FDR]["Gene_Symbol"])
    amp=cd[cd["FDR"]<=FDR].groupby("Gene_Symbol")["rAMP"].max()
    thr=amp.quantile(1-AMP_FRAC); hi=set(amp[amp>=thr].index)
    ia=df["is_anchor"]==1
    pos=df[df["gene_symbol"].isin(hi)&~ia]
    negp=df[(df[LVL]=="VERY_LOW")&~ia&~df["gene_symbol"].isin(circ)]
    neg=negp.sample(n=len(pos),random_state=RS)
    tr=pd.concat([pos,neg]); y=np.r_[np.ones(len(pos),int),np.zeros(len(neg),int)]
    return tr,y,len(pos)

def fig1_roc(df,cols):
    tr,y,npos=get_train(df,cols)
    pipe=Pipeline([("i",SimpleImputer(strategy="median")),("s",StandardScaler()),
                   ("m",LogisticRegression(max_iter=3000,class_weight="balanced",penalty="l1",solver="liblinear"))])
    oof=cross_val_predict(pipe,tr[cols].to_numpy(float),y,cv=StratifiedKFold(5,shuffle=True,random_state=RS),
                          method="predict_proba",n_jobs=-1)[:,1]
    auc=roc_auc_score(y,oof); fpr,tpr,_=roc_curve(y,oof)
    fig,ax=plt.subplots(figsize=(6,6))
    ax.plot(fpr,tpr,color=PURPLE,lw=2.5,label=f"Final model (AUC = {auc:.3f})")
    ax.plot([0,1],[0,1],ls="--",color=GREY,lw=1,label="Random (0.500)")
    ax.set_xlabel("False positive rate"); ax.set_ylabel("True positive rate")
    ax.set_title("Performance of the Final Circadian Gene Classifier")
    ax.legend(loc="lower right",frameon=False)
    ax.text(0.55,0.10,f"CircaDB high-amplitude positives\nTCGA co-expr + methylation\n(robust mean 0.805±0.018)",
            fontsize=8,color="#555555")
    plt.tight_layout(); plt.savefig(f"{OUTDIR}/fig1_roc_v3.png"); plt.close()
    print(f"  fig1_roc_v3.png (this run AUC={auc:.3f})")

def fig2_importance(df,cols):
    tr,y,npos=get_train(df,cols)
    pipe=Pipeline([("i",SimpleImputer(strategy="median")),("s",StandardScaler()),
                   ("m",LogisticRegression(max_iter=3000,class_weight="balanced",penalty="l1",solver="liblinear"))])
    pipe.fit(tr[cols].to_numpy(float),y)
    coef=np.abs(pipe.named_steps["m"].coef_[0])
    imp=pd.DataFrame({"feature":cols,"imp":coef}).sort_values("imp",ascending=False).head(15)
    # 特徵名美化
    def pretty(c):
        return (c.replace("radj_","co:").replace("meth_radj_","me:")
                 .replace("d_posneg","deriv:pos-neg").replace("d_mean","deriv:mean")
                 .replace("d_max","deriv:max").replace("d_sign","deriv:sign"))
    imp["label"]=imp["feature"].map(pretty)
    fig,ax=plt.subplots(figsize=(7,6))
    ax.barh(imp["label"][::-1],imp["imp"][::-1],color=PURPLE)
    ax.set_xlabel("|L1-logistic coefficient|")
    ax.set_title("Top feature importance (final model)")
    plt.tight_layout(); plt.savefig(f"{OUTDIR}/fig2_feature_importance_v3.png"); plt.close()
    print("  fig2_feature_importance_v3.png")

def fig6_strategy():
    # 誠實版策略對照（數字全部來自實際跑出的 report）
    data=[
        ("Final model",                 0.805, "clean"),
        ("Baseline",                    0.764, "clean"),
        ("+ TCGA adj-normal",           0.765, "clean"),
        ("CGDB positives",              0.646, "clean"),
        ("GO circadian positives",      0.610, "clean"),
        ("Methylation (all pos)",       0.600, "clean"),
        ("E-box motif",                 0.539, "clean"),
        ("STRING network",              0.527, "clean"),
        ("GTEx features",               0.871, "circular"),
    ]
    data=sorted(data,key=lambda x:x[1])
    labels=[d[0] for d in data]; aucs=[d[1] for d in data]
    colors=[RED if d[2]=="circular" else PURPLE for d in data]
    fig,ax=plt.subplots(figsize=(8.5,6))
    bars=ax.barh(labels,aucs,color=colors)
    ax.axvline(0.805,ls="--",color=ORANGE,lw=1.2)
    ax.text(0.805,8.5,"final 0.805",color="#b45309",fontsize=9,ha="center",va="bottom")
    ax.axvline(0.5,ls=":",color=GREY,lw=1)
    for b,a in zip(bars,aucs):
        ax.text(a+0.005,b.get_y()+b.get_height()/2,f"{a:.3f}",va="center",fontsize=8)
    ax.set_xlim(0.45,0.92); ax.set_xlabel("Cross-validated ROC-AUC")
    ax.set_title("Comparison of Machine-Learning Strategies")
    # 圖例
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=PURPLE,label="Final models"),
                       Patch(color=RED,label="Excluded (information leakage)")],
              loc="lower center",bbox_to_anchor=(0.52,0.02),
              frameon=True,framealpha=0.9,fontsize=9)
    plt.tight_layout(); plt.savefig(f"{OUTDIR}/fig6_strategy_auc_v3.png"); plt.close()
    print("  fig6_strategy_auc_v3.png")

def main():
    print("重繪 fig1/2/6 突破版 → figures/")
    df,cols=build()
    try: fig1_roc(df,cols)
    except Exception as e: print("[fig1 跳過]",e)
    try: fig2_importance(df,cols)
    except Exception as e: print("[fig2 跳過]",e)
    try: fig6_strategy()
    except Exception as e: print("[fig6 跳過]",e)
    print("完成。")

if __name__=="__main__":
    main()