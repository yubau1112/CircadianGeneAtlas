#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CGA MS3 — 定版候選名單（突破版，AUC≈0.80，無循環）
================================================================
定版模型：CircaDB 高幅度正例(rAMP 前 40%) + coexpr+deriv+meth 特徵 + Logit_L1。
穩健驗證：15 seeds 平均 AUC 0.805±0.018（非單次 0.847）。無 GTEx、無循環。

對全基因打分，挑「高 ai_score + 低文獻」乾淨候選。
輸出：step_final_candidates.csv（定版候選）、step_final_scores.csv（全體分數）
================================================================
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, average_precision_score

COEXPR="coexpr_features_adj.csv"; METH="meth_features_adj.csv"; CIRCADB="CircaDB_HA_query.csv"
OUT_C="step_final_candidates.csv"; OUT_S="step_final_scores.csv"
LVL,RS,FDR="annotation_confidence_level",42,0.05
AI_PCT=0.95; LOW={"VERY_LOW","LOW"}; AMP_FRAC=0.40   # 高幅度前 40%

report=[]
def log(m): print(m); report.append(m)

def main():
    co=pd.read_csv(COEXPR); me=pd.read_csv(METH)
    coex=[c for c in co.columns if c.startswith("radj_")]
    meth=[c for c in me.columns if c.startswith("meth_radj_")]
    df=co.merge(me[["gene_symbol"]+meth],on="gene_symbol",how="left")
    # 衍生特徵
    R=df[coex].to_numpy(float)
    POS=["radj_ARNTL","radj_CLOCK","radj_NPAS2","radj_RORA","radj_RORB","radj_RORC","radj_ARNTL2"]
    NEG=["radj_PER1","radj_PER2","radj_PER3","radj_CRY1","radj_CRY2","radj_NR1D1","radj_NR1D2","radj_BHLHE40","radj_BHLHE41"]
    df["d_posneg"]=np.nanmean(df[[c for c in POS if c in coex]],axis=1)-np.nanmean(df[[c for c in NEG if c in coex]],axis=1)
    df["d_mean"]=np.nanmean(R,axis=1); df["d_max"]=np.nanmax(R,axis=1); df["d_sign"]=np.nanmean(np.sign(R),axis=1)
    cols=coex+["d_posneg","d_mean","d_max","d_sign"]+meth

    cd=pd.read_csv(CIRCADB); circ_all=cd[cd["FDR"]<=FDR].copy()
    circ_all["Gene_Symbol"]=circ_all["Gene_Symbol"].astype(str).str.strip()
    circ=set(circ_all["Gene_Symbol"])
    amp=circ_all.groupby("Gene_Symbol")["rAMP"].max()
    thr=amp.quantile(1-AMP_FRAC); hi=set(amp[amp>=thr].index)
    log(f"CircaDB 全部={len(circ)}；高幅度(前{int(AMP_FRAC*100)}%)={len(hi)}")

    ia=df["is_anchor"]==1
    in_hi=df["gene_symbol"].isin(hi)
    pos=df[in_hi & ~ia]
    negp=df[(df[LVL]=="VERY_LOW") & ~ia & ~df["gene_symbol"].isin(circ)]
    neg=negp.sample(n=len(pos),random_state=RS)
    tr=pd.concat([pos,neg]); y=np.r_[np.ones(len(pos),int),np.zeros(len(neg),int)]
    log(f"正例(高幅度,排除anchor)={len(pos)}  負例={len(neg)}")

    pipe=Pipeline([("i",SimpleImputer(strategy="median")),("s",StandardScaler()),
                   ("m",LogisticRegression(max_iter=3000,class_weight="balanced",penalty="l1",solver="liblinear"))])
    oof=cross_val_predict(pipe,tr[cols].to_numpy(float),y,
         cv=StratifiedKFold(5,shuffle=True,random_state=RS),method="predict_proba",n_jobs=-1)[:,1]
    log(f"此次 CV：ROC-AUC={roc_auc_score(y,oof):.3f} PR-AUC={average_precision_score(y,oof):.3f}（穩健平均≈0.805±0.018）\n")

    # 全基因打分
    pipe.fit(tr[cols].to_numpy(float),y)
    df["ai_score"]=pipe.predict_proba(df[cols].to_numpy(float))[:,1]
    df["is_circadb"]=df["gene_symbol"].isin(circ).astype(int)
    df["is_hi_amp"]=in_hi.astype(int)
    df[["gene_symbol","ai_score","is_anchor","is_circadb","is_hi_amp",
        "annotation_combined_score","annotation_confidence_level"]]\
      .sort_values("ai_score",ascending=False).to_csv(OUT_S,index=False,encoding="utf-8-sig")

    # 乾淨候選：高 ai + 低文獻 + 非 anchor + 不在 CircaDB
    cut=df["ai_score"].quantile(AI_PCT)
    cand=df[(df["ai_score"]>=cut)&(df[LVL].isin(LOW))&(~ia)&(~df["gene_symbol"].isin(circ))].sort_values("ai_score",ascending=False)
    cand.to_csv(OUT_C,index=False,encoding="utf-8-sig")
    log(f"定版候選數（ai 前{int((1-AI_PCT)*100)}% ∩ 低文獻 ∩ 非anchor ∩ 非CircaDB）：{len(cand)}")

    # 跟之前乾淨版(401)比重疊
    try:
        prev=set(pd.read_csv("step3_candidates_clean.csv")["gene_symbol"])
        now=set(cand["gene_symbol"]); both=prev&now
        log(f"與前乾淨版(401)重疊：{len(both)}（前{len(prev)}、今{len(now)}）")
    except Exception as e:
        log(f"[重疊對照跳過] {e}")

    log("\nTop 25 定版候選：")
    for i,(_,r) in enumerate(cand.head(25).iterrows(),1):
        log(f"  {i:>2} {r['gene_symbol']:<12} ai={r['ai_score']:.3f}  文獻={r['annotation_combined_score']}  {r[LVL]}")
    open("step_final_candidates_report.txt","w",encoding="utf-8").write("\n".join(report))
    print(f"\n[輸出] {OUT_C}、{OUT_S}、step_final_candidates_report.txt")

if __name__=="__main__":
    main()