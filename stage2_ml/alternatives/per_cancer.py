#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CGA MS3 — Step 2f：單癌深挖 / 跨癌最強訊號特徵（測聚合稀釋假設）
================================================================
假設：24 癌 Fisher-z 平均把訊號稀釋了。改用不被平均沖淡的特徵：
  - max signed rho：每基因×每錨點，取跨癌「絕對值最大」的 rho（保符號）
  - top-k mean：取跨癌最強 k 個癌的平均
  - 單癌：樣本最多的幾個癌別各自的 rho（不聚合）
全部背景校正（row-center + col z-score）後，比 AUC。
正例 = CircaDB（無循環）。

輸入：coexpr_per_cancer_long.csv、coexpr_features_adj.csv（拿 label/anchor）、CircaDB_HA_query.csv
輸出：step2f_report.txt
================================================================
"""
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, average_precision_score

LONG   = "coexpr_per_cancer_long.csv"
FEAT   = "coexpr_features_adj.csv"
CIRCADB= "CircaDB_HA_query.csv"
OUT    = "step2f_report.txt"
LVL, RS, FDR = "annotation_confidence_level", 42, 0.05
TOPK = 5                      # 跨癌取最強 5 癌平均
BIG_CANCERS = ["TCGA-BRCA","TCGA-UCEC","TCGA-KIRC","TCGA-LUAD","TCGA-THCA"]  # 樣本最多幾個

report=[]
def log(m):
    print(m); report.append(m)

def bg_correct(mat_df):
    """row-center + column z-score。輸入 index=gene, columns=anchor 的 rho 表。"""
    R = mat_df.to_numpy(float)
    Rc = R - np.nanmean(R, axis=1, keepdims=True)
    Rz = (Rc - np.nanmean(Rc, axis=0)) / (np.nanstd(Rc, axis=0) + 1e-9)
    return pd.DataFrame(Rz, index=mat_df.index, columns=[c+"_z" for c in mat_df.columns])

def cvauc(X, y, model):
    pipe = Pipeline([("imp",SimpleImputer(strategy="median")),("sc",StandardScaler()),("clf",model)])
    cv = StratifiedKFold(5, shuffle=True, random_state=RS)
    oof = cross_val_predict(pipe,X,y,cv=cv,method="predict_proba",n_jobs=-1)[:,1]
    return roc_auc_score(y,oof), average_precision_score(y,oof)

def main():
    log("讀長檔（523 萬列，稍等）...")
    lg = pd.read_csv(LONG)  # gene_symbol, anchor, cancer, n_samples, rho
    log(f"長檔欄位：{list(lg.columns)}  列數：{len(lg):,}")

    # --- 特徵 1：跨癌最強 signed rho（絕對值最大者，保符號）---
    idx = lg.groupby(["gene_symbol","anchor"])["rho"].apply(lambda s: s.loc[s.abs().idxmax()])
    max_rho = idx.unstack("anchor")
    # --- 特徵 2：跨癌 top-k 平均（依 |rho| 取前 k）---
    def topk_mean(s):
        return s.reindex(s.abs().sort_values(ascending=False).index).head(TOPK).mean()
    topk = lg.groupby(["gene_symbol","anchor"])["rho"].apply(topk_mean).unstack("anchor")
    # --- 特徵 3：單一大癌（BRCA）---
    brca = lg[lg["cancer"]=="TCGA-BRCA"].pivot_table(index="gene_symbol", columns="anchor", values="rho")

    # 背景校正
    Fmax = bg_correct(max_rho)
    Ftopk = bg_correct(topk); Ftopk.columns=[c+"_tk" for c in Ftopk.columns]
    Fbrca = bg_correct(brca); Fbrca.columns=[c+"_brca" for c in Fbrca.columns]

    # 標籤
    base = pd.read_csv(FEAT)[["gene_symbol","is_anchor",LVL]]
    cd = pd.read_csv(CIRCADB)
    circ = set(cd[cd["FDR"]<=FDR]["Gene_Symbol"].astype(str).str.strip())
    base["in_circ"] = base["gene_symbol"].isin(circ).astype(int)

    def build(feat_df, name):
        d = base.merge(feat_df.reset_index(), on="gene_symbol", how="inner")
        pos = d[(d["in_circ"]==1) & (d["is_anchor"]==0)]
        negp = d[(d[LVL]=="VERY_LOW") & (d["is_anchor"]==0) & (d["in_circ"]==0)]
        neg = negp.sample(n=len(pos), random_state=RS)
        tr = pd.concat([pos,neg]); 
        cols = [c for c in feat_df.columns]
        X = tr[cols].to_numpy(float); y = np.r_[np.ones(len(pos),int),np.zeros(len(neg),int)]
        return X, y, name, len(pos)

    log("\n" + "="*60)
    log("Step 2f：單癌深挖 / 跨癌最強訊號（正例=CircaDB）")
    log("="*60)
    log(f"{'特徵集 / 模型':<34}{'正例':>5}{'ROC-AUC':>9}{'PR-AUC':>8}")
    log("-"*60)
    models = {"Logit":LogisticRegression(max_iter=2000,class_weight="balanced"),
              "RF":RandomForestClassifier(n_estimators=500,random_state=RS,class_weight="balanced",n_jobs=-1)}
    best=(0,"")
    for feat_df, nm in [(Fmax,"跨癌最強(max)"),(Ftopk,f"跨癌top{TOPK}平均"),(Fbrca,"單癌BRCA"),
                        (pd.concat([Fmax,Ftopk],axis=1),"max+topk")]:
        X,y,name,npos = build(feat_df, nm)
        for mn,md in models.items():
            auc,ap = cvauc(X,y,md)
            log(f"{nm+' / '+mn:<34}{npos:>5}{auc:>9.3f}{ap:>8.3f}")
            if auc>best[0]: best=(auc,f"{nm}/{mn}")
    log("-"*60)
    log(f"最佳：AUC={best[0]:.3f} → {best[1]}")
    log("對照：24癌平均聚合 Logistic = 0.756")
    log("判讀：若明顯 >0.756 → 聚合稀釋是主因，改用此特徵；若≈0.756 → 稀釋非瓶頸。")
    with open(OUT,"w",encoding="utf-8") as f: f.write("\n".join(report))
    print(f"\n[輸出] {OUT}")

if __name__=="__main__":
    main()