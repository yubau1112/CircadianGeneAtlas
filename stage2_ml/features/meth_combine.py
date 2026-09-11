#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CGA MS3 — 共表達 + 甲基化 合併分類器（測表觀證據，無循環）
================================================================
把 TCGA 啟動子甲基化特徵（meth_radj_，與 GTEx/CircaDB 獨立）併進 TCGA 共表達
（radj_），看能否乾淨地把 0.756 往上推。正例 = CircaDB（無循環）。

比較：只共表達 / 只甲基化 / 兩者合併 × Logistic/RF。
若合併明顯 > 0.756 → 表觀證據有效，加入定版。

輸入：coexpr_features_adj.csv、meth_features_adj.csv（全量）、CircaDB_HA_query.csv
輸出：step_meth_report.txt、step_meth_scores_all_genes.csv（最佳模型全體分數）
================================================================
"""
import numpy as np, pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, average_precision_score

COEXPR="coexpr_features_adj.csv"; METH="meth_features_adj.csv"; CIRCADB="CircaDB_HA_query.csv"
OUT_R="step_meth_report.txt"; OUT_S="step_meth_scores_all_genes.csv"
LVL,RS,FDR="annotation_confidence_level",42,0.05
AI_PCT=0.95; LOW={"VERY_LOW","LOW"}

report=[]
def log(m): print(m); report.append(m)

def main():
    co=pd.read_csv(COEXPR); me=pd.read_csv(METH)
    coex=[c for c in co.columns if c.startswith("radj_")]
    meth=[c for c in me.columns if c.startswith("meth_radj_")]
    # 確認甲基化是全量(非煙霧測試436列)
    log(f"甲基化特徵表列數：{len(me)}（全量應約 10000+；若 ~436 為煙霧測試檔，請先跑全量）")

    df=co.merge(me[["gene_symbol"]+meth],on="gene_symbol",how="inner")
    log(f"共表達∩甲基化 共同基因：{len(df)}")

    cd=pd.read_csv(CIRCADB); circ=set(cd[cd["FDR"]<=FDR]["Gene_Symbol"].astype(str).str.strip())
    ia=df["is_anchor"]==1; ic=df["gene_symbol"].isin(circ)
    pos=df[ic&~ia].copy(); negp=df[(df[LVL]=="VERY_LOW")&~ia&~ic]
    neg=negp.sample(n=len(pos),random_state=RS)
    tr=pd.concat([pos,neg],ignore_index=True); y=np.r_[np.ones(len(pos),int),np.zeros(len(neg),int)]
    log(f"正例(CircaDB,排除anchor)={len(pos)}  負例(VERY_LOW)={len(neg)}\n")

    def cvauc(cols,mdl):
        pipe=Pipeline([("i",SimpleImputer(strategy="median")),("s",StandardScaler()),("m",mdl)])
        oof=cross_val_predict(pipe,tr[cols].to_numpy(float),y,
             cv=StratifiedKFold(5,shuffle=True,random_state=RS),method="predict_proba",n_jobs=-1)[:,1]
        return roc_auc_score(y,oof),average_precision_score(y,oof)

    log(f"{'特徵集 / 模型':<32}{'ROC-AUC':>9}{'PR-AUC':>8}")
    log("-"*52)
    sets={"只共表達(radj)":coex,"只甲基化(meth)":meth,"共表達+甲基化":coex+meth}
    models={"Logit":LogisticRegression(max_iter=2000,class_weight="balanced"),
            "RF":RandomForestClassifier(n_estimators=500,random_state=RS,class_weight="balanced",n_jobs=-1)}
    best=(0,None,None,None)
    for sn,cols in sets.items():
        for mn,md in models.items():
            auc,ap=cvauc(cols,md); log(f"{sn+' / '+mn:<32}{auc:>9.3f}{ap:>8.3f}")
            if auc>best[0]: best=(auc,sn,mn,cols)
    log("-"*52)
    log(f"最佳：AUC={best[0]:.3f} → {best[1]} / {best[2]}")
    log("對照：只共表達(乾淨基準)Logistic = 0.756")
    log("判讀：合併明顯 >0.756 → 甲基化(表觀)有效、無循環，加入定版；≈0.756 → 幫助有限。\n")

    # 最佳組合對全體打分 + 乾淨候選
    _,sn,mn,cols=best
    pipe=Pipeline([("i",SimpleImputer(strategy="median")),("s",StandardScaler()),("m",models[mn])])
    pipe.fit(tr[cols].to_numpy(float),y)
    df["ai_score"]=pipe.predict_proba(df[cols].to_numpy(float))[:,1]; df["is_circadb"]=ic.astype(int)
    df[["gene_symbol","ai_score","is_anchor","is_circadb","annotation_combined_score","annotation_confidence_level"]]\
      .sort_values("ai_score",ascending=False).to_csv(OUT_S,index=False,encoding="utf-8-sig")
    cut=df["ai_score"].quantile(AI_PCT)
    cand=df[(df["ai_score"]>=cut)&(df[LVL].isin(LOW))&(~ia)&(~ic)].sort_values("ai_score",ascending=False)
    log(f"最佳模型乾淨候選數（前 {(1-AI_PCT)*100:.0f}%）：{len(cand)}")
    log("Top 15：")
    for i,(_,r) in enumerate(cand.head(15).iterrows(),1):
        log(f"  {i:>2} {r['gene_symbol']:<12} ai={r['ai_score']:.3f}  文獻={r['annotation_combined_score']}")
    open(OUT_R,"w",encoding="utf-8").write("\n".join(report))
    print(f"\n[輸出] {OUT_R}、{OUT_S}")

if __name__=="__main__":
    main()