#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 9/11 修改（Checklist 第 29 項）：圖上文字統一為 clock-associated prediction score、移除內部稱呼 breakthrough／clean model、Known rhythmic 改 CircaDB rhythmic；輸出改新檔名，其餘不變
"""
CGA MS3 — fig8 + fig9 突破版重繪（定版模型,高幅度正例164）
================================================================
取代舊版 fig8/fig9（那版用全部正例419,與定版不一致）。
用定版模型（CircaDB高幅度正例rAMP前40% + coexpr+deriv+meth + L1-Logistic）重跑。

fig8_pr_curve_v3：PR 曲線（out-of-fold, AP）
fig9_confusion_v3：混淆矩陣（out-of-fold, 0.5門檻）

輸入：coexpr_features_adj.csv、meth_features_adj.csv、CircaDB_HA_query.csv
輸出：figures/fig8_pr_curve_v3.png、fig9_confusion_v3.png
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
from sklearn.metrics import average_precision_score, precision_recall_curve, confusion_matrix, roc_auc_score

OUTDIR="figures"; os.makedirs(OUTDIR,exist_ok=True)
PURPLE,ORANGE,GREY="#3b0764","#f59e0b","#9ca3af"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":11,"axes.spines.top":False,
                     "axes.spines.right":False,"figure.dpi":300})
COEXPR="coexpr_features_adj.csv"; METH="meth_features_adj.csv"; CIRCADB="CircaDB_HA_query.csv"
LVL,RS,FDR="annotation_confidence_level",42,0.05; AMP_FRAC=0.40

def get_oof():
    co=pd.read_csv(COEXPR); me=pd.read_csv(METH)
    coex=[c for c in co.columns if c.startswith("radj_")]
    meth=[c for c in me.columns if c.startswith("meth_radj_")]
    d=co.merge(me[["gene_symbol"]+meth],on="gene_symbol",how="left")
    R=d[coex].to_numpy(float)
    POS=["radj_ARNTL","radj_CLOCK","radj_NPAS2","radj_RORA","radj_RORB","radj_RORC","radj_ARNTL2"]
    NEG=["radj_PER1","radj_PER2","radj_PER3","radj_CRY1","radj_CRY2","radj_NR1D1","radj_NR1D2","radj_BHLHE40","radj_BHLHE41"]
    d["d_posneg"]=np.nanmean(d[[c for c in POS if c in coex]],axis=1)-np.nanmean(d[[c for c in NEG if c in coex]],axis=1)
    d["d_mean"]=np.nanmean(R,axis=1); d["d_max"]=np.nanmax(R,axis=1); d["d_sign"]=np.nanmean(np.sign(R),axis=1)
    cols=coex+["d_posneg","d_mean","d_max","d_sign"]+meth
    cd=pd.read_csv(CIRCADB); cd["Gene_Symbol"]=cd["Gene_Symbol"].astype(str).str.strip()
    circ=set(cd[cd["FDR"]<=FDR]["Gene_Symbol"])
    amp=cd[cd["FDR"]<=FDR].groupby("Gene_Symbol")["rAMP"].max()
    thr=amp.quantile(1-AMP_FRAC); hi=set(amp[amp>=thr].index)
    ia=d["is_anchor"]==1
    pos=d[d["gene_symbol"].isin(hi)&~ia]
    negp=d[(d[LVL]=="VERY_LOW")&~ia&~d["gene_symbol"].isin(circ)]
    neg=negp.sample(n=len(pos),random_state=RS)
    tr=pd.concat([pos,neg]); y=np.r_[np.ones(len(pos),int),np.zeros(len(neg),int)]
    pipe=Pipeline([("i",SimpleImputer(strategy="median")),("s",StandardScaler()),
                   ("m",LogisticRegression(max_iter=3000,class_weight="balanced",penalty="l1",solver="liblinear"))])
    oof=cross_val_predict(pipe,tr[cols].to_numpy(float),y,cv=StratifiedKFold(5,shuffle=True,random_state=RS),
                          method="predict_proba",n_jobs=-1)[:,1]
    return y,oof,len(pos)

def main():
    print("重繪 fig8/fig9 突破版(定版,高幅度164) → figures/")
    y,oof,npos=get_oof()
    ap=average_precision_score(y,oof); auc=roc_auc_score(y,oof)
    print(f"  正例={npos} 負例={npos}  單次 AUC={auc:.3f}  AP={ap:.3f}")

    # fig8 PR 曲線
    prec,rec,_=precision_recall_curve(y,oof)
    fig,ax=plt.subplots(figsize=(6.5,6.5))
    ax.plot(rec,prec,color=PURPLE,lw=2.5,label=f"Final model (AP = {ap:.3f})")
    ax.axhline(0.5,ls="--",color=GREY,lw=1,label="Baseline (0.500)")
    ax.set_xlabel("Recall"); ax.set_ylabel("Precision")
    ax.set_title("Precision–recall curve (out-of-fold)")
    ax.set_ylim(0.45,1.02)
    ax.legend(loc="upper right",frameon=True,framealpha=0.9,fontsize=9)
    ax.text(0.04,0.515,f"CircaDB high-amplitude positives, n={npos} per class",
            fontsize=8,color="#777777",va="bottom")
    ax.set_xlim(-0.02,1.02)
    plt.tight_layout(); plt.savefig(f"{OUTDIR}/fig8_pr_curve_v3.png"); plt.close()
    print("  fig8_pr_curve_v3.png")

    # fig9 混淆矩陣
    pred=(oof>=0.5).astype(int)
    tn,fp,fn,tp=confusion_matrix(y,pred).ravel()
    acc=(tp+tn)/len(y); prec_=tp/(tp+fp); rec_=tp/(tp+fn); spec=tn/(tn+fp)
    cm=np.array([[tn,fp],[fn,tp]])
    fig,ax=plt.subplots(figsize=(7,6.5))
    im=ax.imshow(cm,cmap="Purples",vmin=0)
    ax.set_aspect("equal")
    labels=[["TN","FP"],["FN","TP"]]
    for i in range(2):
        for j in range(2):
            v=cm[i,j]; txtcol="white" if v>cm.max()*0.5 else "#3b0764"
            ax.text(j,i,f"{labels[i][j]}\n{v}",ha="center",va="center",fontsize=12,color=txtcol,fontweight="bold")
    ax.set_xticks([0,1]); ax.set_xticklabels(["Predicted negative","Predicted positive"],fontsize=9)
    ax.set_yticks([0,1]); ax.set_yticklabels(["Negative\n(VERY_LOW)","Positive\n(CircaDB)"],fontsize=9)
    ax.set_title(f"Confusion matrix (out-of-fold)\nacc {acc:.2f} · prec {prec_:.2f} · recall {rec_:.2f} · spec {spec:.2f}",fontsize=10,pad=12)
    plt.subplots_adjust(left=0.22,right=0.82,top=0.85,bottom=0.15)
    plt.savefig(f"{OUTDIR}/fig9_confusion_v3.png"); plt.close()
    print("  fig9_confusion_v3.png")
    print(f"\n=== 給圖說更新用的數字 ===")
    print(f"fig8: AP = {ap:.3f}")
    print(f"fig9: TP={tp} TN={tn} FP={fp} FN={fn}; acc={acc:.3f} prec={prec_:.3f} recall={rec_:.3f}")
    print(f"     樣本總數 = {len(y)}（{npos} 正 + {npos} 負）")
    print("完成。")

if __name__=="__main__":
    main()