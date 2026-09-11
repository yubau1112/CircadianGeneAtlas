#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CGA MS3 — 穩健性驗證：0.847 是真突破還是少樣本運氣？
================================================================
最佳組合 = coexpr+deriv+meth / CircaDB高幅度 / Logit_L1。
驗證：
  A. 固定門檻(rAMP前50%)，跑 15 個隨機種子 → AUC 平均±標準差、範圍
  B. 掃高幅度門檻(前40/50/60/70%)，各跑 10 種子 → 看門檻敏感度
  C. 對照：CircaDB全部(419)也跑 15 種子 → 確認 0.756 穩不穩
判準：若最佳組合 AUC 穩定在高位(標準差小、下限仍明顯>0.756)＝真突破。

輸入：coexpr_features_adj.csv、meth_features_adj.csv、CircaDB_HA_query.csv
輸出：step_robust_report.txt
================================================================
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score

COEXPR="coexpr_features_adj.csv"; METH="meth_features_adj.csv"; CIRCADB="CircaDB_HA_query.csv"
OUT="step_robust_report.txt"; LVL,FDR="annotation_confidence_level",0.05

report=[]
def log(m): print(m); report.append(m)

def build_df():
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

def run_once(df,cols,posset,circ,seed):
    ia=df["is_anchor"]==1
    pos=df[df["gene_symbol"].isin(posset) & ~ia]
    negp=df[(df[LVL]=="VERY_LOW") & ~ia & ~df["gene_symbol"].isin(circ)]
    if len(pos)<30 or len(negp)<len(pos): return None,0
    neg=negp.sample(n=len(pos),random_state=seed)
    tr=pd.concat([pos,neg]); X=tr[cols].to_numpy(float)
    y=np.r_[np.ones(len(pos),int),np.zeros(len(neg),int)]
    pipe=Pipeline([("i",SimpleImputer(strategy="median")),("s",StandardScaler()),
                   ("m",LogisticRegression(max_iter=3000,class_weight="balanced",penalty="l1",solver="liblinear"))])
    cv=StratifiedKFold(5,shuffle=True,random_state=seed)
    oof=cross_val_predict(pipe,X,y,cv=cv,method="predict_proba",n_jobs=-1)[:,1]
    return roc_auc_score(y,oof),len(pos)

def main():
    df,cols=build_df()
    cd=pd.read_csv(CIRCADB); circ_all=cd[cd["FDR"]<=FDR].copy()
    circ=set(circ_all["Gene_Symbol"].astype(str).str.strip())
    circ_all["Gene_Symbol"]=circ_all["Gene_Symbol"].astype(str).str.strip()
    # 每基因取最大 rAMP（同基因多組織）
    amp=circ_all.groupby("Gene_Symbol")["rAMP"].max()

    log("="*60); log("穩健性驗證：最佳組合 coexpr+deriv+meth / Logit_L1"); log("="*60)

    def hi_set(frac):
        thr=amp.quantile(1-frac)
        return set(amp[amp>=thr].index)

    # A. 固定門檻前50%，15 seeds
    log("\n[A] CircaDB高幅度(前50%)，15 個隨機種子：")
    aucs=[]
    for s in range(15):
        a,n=run_once(df,cols,hi_set(0.50),circ,s)
        if a: aucs.append(a)
    aucs=np.array(aucs)
    log(f"  正例數≈{n}  AUC: 平均={aucs.mean():.3f} 標準差={aucs.std():.3f} 範圍[{aucs.min():.3f},{aucs.max():.3f}]")

    # B. 掃門檻
    log("\n[B] 不同高幅度門檻(各 10 種子)：")
    log(f"  {'門檻':<10}{'正例數':>6}{'AUC平均':>9}{'標準差':>8}")
    for frac in [0.40,0.50,0.60,0.70,1.00]:
        aa=[]
        for s in range(10):
            a,n=run_once(df,cols,hi_set(frac) if frac<1.0 else circ,circ,s)
            if a: aa.append(a)
        aa=np.array(aa)
        label=f"前{int(frac*100)}%" if frac<1.0 else "全部"
        log(f"  {label:<10}{n:>6}{aa.mean():>9.3f}{aa.std():>8.3f}")

    # C. 對照：全部 CircaDB，15 seeds（只用共表達，對照原始 0.756）
    log("\n[C] 對照 — CircaDB全部、只共表達特徵(原始基準)，15 種子：")
    coex=[c for c in df.columns if c.startswith("radj_")]
    aucs2=[]
    for s in range(15):
        a,n=run_once(df,coex,circ,circ,s)
        if a: aucs2.append(a)
    aucs2=np.array(aucs2)
    log(f"  正例數≈{n}  AUC: 平均={aucs2.mean():.3f} 標準差={aucs2.std():.3f} 範圍[{aucs2.min():.3f},{aucs2.max():.3f}]")

    log("\n"+"="*60)
    log("判讀：")
    log(f"  · 最佳組合(前50%)平均 {aucs.mean():.3f}±{aucs.std():.3f}，下限 {aucs.min():.3f}")
    log(f"  · 原始基準(全部/共表達)平均 {aucs2.mean():.3f}±{aucs2.std():.3f}")
    if aucs.mean()-2*aucs.std() > aucs2.mean():
        log("  → 最佳組合下限明顯高於基準：真突破，穩健。")
    else:
        log("  → 兩者信賴區間重疊：提升可能部分來自少樣本變異，需保守解讀。")
    open(OUT,"w",encoding="utf-8").write("\n".join(report))
    print(f"\n[輸出] {OUT}")

if __name__=="__main__":
    main()