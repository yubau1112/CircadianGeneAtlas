#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CGA MS3 — CGDB 正例 + 重新啟用 GTEx 特徵（救回 GTEx，無循環）
================================================================
關鍵想法：GTEx 共表達訊號很強（0.867），但之前因 CircaDB(GTEx) 標籤而循環。
改用 CGDB 實驗節律基因當正例（功能註釋，與 GTEx 完全獨立），GTEx 特徵就乾淨了。

比較：TCGA / GTEx / 合併 × Logistic/RF，正例 = GO circadian（排除 anchor）。
若 GTEx 或合併明顯 > 0.756 → 成功救回 GTEx 的強訊號，且無循環。

輸入：coexpr_features_adj.csv、gtex_features_adj.csv、cgdb_human_genes.txt、CircaDB_HA_query.csv
輸出：step_cgdb_report.txt、step_cgdb_scores_all_genes.csv
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

TCGA_F="coexpr_features_adj.csv"; GTEX_F="gtex_features_adj.csv"
GO_LIST="cgdb_human_genes.txt"; CIRCADB="CircaDB_HA_query.csv"
OUT_R="step_cgdb_report.txt"; OUT_S="step_cgdb_scores_all_genes.csv"
LVL,RS,FDR="annotation_confidence_level",42,0.05
AI_PCT=0.95; LOW={"VERY_LOW","LOW"}

# GO 清單可能用舊別名，轉成官方符號（與錨點命名一致）
ALIAS={"BMAL1":"ARNTL","BMAL2":"ARNTL2","REV-ERBA":"NR1D1","REV-ERBB":"NR1D2",
       "DEC1":"BHLHE40","DEC2":"BHLHE41"}
ANCHORS={"CLOCK","ARNTL","ARNTL2","PER1","PER2","PER3","CRY1","CRY2","NR1D1","NR1D2",
         "RORA","RORB","RORC","NPAS2","CSNK1D","CSNK1E","TIMELESS","BHLHE40","BHLHE41"}

report=[]
def log(m): print(m); report.append(m)

def main():
    t=pd.read_csv(TCGA_F); g=pd.read_csv(GTEX_F)
    tcga=[c for c in t.columns if c.startswith("radj_")]
    gtex=[c for c in g.columns if c.startswith("gtex_radj_")]
    df=t.merge(g[["gene_symbol"]+gtex],on="gene_symbol",how="inner")

    # CGDB 正例（去別名、排除 anchor）
    go=set()
    with open(GO_LIST,encoding="utf-8") as f:
        for line in f:
            s=line.strip()
            if s: go.add(ALIAS.get(s,s))
    cgdb_noanchor=go-ANCHORS
    log(f"CGDB 節律基因：{len(go)}；排除 {len(go)-len(cgdb_noanchor)} 個 anchor → 正例池 {len(cgdb_noanchor)}")

    cd=pd.read_csv(CIRCADB); circ=set(cd[cd["FDR"]<=FDR]["Gene_Symbol"].astype(str).str.strip())
    ia=df["is_anchor"]==1
    in_cgdb=df["gene_symbol"].isin(cgdb_noanchor)
    pos=df[in_cgdb & ~ia].copy()
    # 負例：VERY_LOW，且不在 GO、不在 CircaDB（乾淨負例）
    negp=df[(df[LVL]=="VERY_LOW") & ~ia & ~df["gene_symbol"].isin(go) & ~df["gene_symbol"].isin(circ)]
    neg=negp.sample(n=len(pos),random_state=RS)
    tr=pd.concat([pos,neg],ignore_index=True); y=np.r_[np.ones(len(pos),int),np.zeros(len(neg),int)]
    log(f"正例(GO∩特徵表,排除anchor)={len(pos)}  負例(VERY_LOW)={len(neg)}\n")

    def cvauc(cols,mdl):
        pipe=Pipeline([("i",SimpleImputer(strategy="median")),("s",StandardScaler()),("m",mdl)])
        oof=cross_val_predict(pipe,tr[cols].to_numpy(float),y,
             cv=StratifiedKFold(5,shuffle=True,random_state=RS),method="predict_proba",n_jobs=-1)[:,1]
        return roc_auc_score(y,oof),average_precision_score(y,oof)

    log(f"{'特徵集 / 模型':<28}{'ROC-AUC':>9}{'PR-AUC':>8}")
    log("-"*48)
    sets={"TCGA":tcga,"GTEx":gtex,"TCGA+GTEx":tcga+gtex}
    models={"Logit":LogisticRegression(max_iter=2000,class_weight="balanced"),
            "RF":RandomForestClassifier(n_estimators=500,random_state=RS,class_weight="balanced",n_jobs=-1)}
    best=(0,None,None,None)
    for sn,cols in sets.items():
        for mn,md in models.items():
            auc,ap=cvauc(cols,md); log(f"{sn+' / '+mn:<28}{auc:>9.3f}{ap:>8.3f}")
            if auc>best[0]: best=(auc,sn,mn,cols)
    log("-"*48)
    log(f"最佳：AUC={best[0]:.3f} → {best[1]} / {best[2]}")
    log("對照：TCGA-only(CircaDB正例)乾淨基準 = 0.756")
    log("判讀：GTEx/合併若明顯>0.756 → 成功救回GTEx強訊號(正例=GO,與GTEx獨立,無循環)。\n")

    # 最佳模型打分 + 乾淨候選
    _,sn,mn,cols=best
    pipe=Pipeline([("i",SimpleImputer(strategy="median")),("s",StandardScaler()),("m",models[mn])])
    pipe.fit(tr[cols].to_numpy(float),y)
    df["ai_score"]=pipe.predict_proba(df[cols].to_numpy(float))[:,1]
    df["in_cgdb"]=in_cgdb.astype(int); df["is_circadb"]=df["gene_symbol"].isin(circ).astype(int)
    df[["gene_symbol","ai_score","is_anchor","in_cgdb","is_circadb","annotation_combined_score","annotation_confidence_level"]]\
      .sort_values("ai_score",ascending=False).to_csv(OUT_S,index=False,encoding="utf-8-sig")
    cut=df["ai_score"].quantile(AI_PCT)
    cand=df[(df["ai_score"]>=cut)&(df[LVL].isin(LOW))&(~ia)&(~in_cgdb)&(~df["gene_symbol"].isin(circ))].sort_values("ai_score",ascending=False)
    log(f"最佳模型乾淨候選數（前 {(1-AI_PCT)*100:.0f}%）：{len(cand)}")
    log("Top 15：")
    for i,(_,r) in enumerate(cand.head(15).iterrows(),1):
        log(f"  {i:>2} {r['gene_symbol']:<12} ai={r['ai_score']:.3f}  文獻={r['annotation_combined_score']}")
    open(OUT_R,"w",encoding="utf-8").write("\n".join(report))
    print(f"\n[輸出] {OUT_R}、{OUT_S}")

if __name__=="__main__":
    main()