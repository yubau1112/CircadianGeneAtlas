#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CGA MS3 — E-box motif 特徵（獨立於共表達的序列證據）
================================================================
時鐘 CLOCK/BMAL1 結合啟動子 E-box（CACGTG）。掃候選基因啟動子的 E-box，
當獨立特徵，看能否再提升合併模型（TCGA+GTEx，AUC 0.871）。

流程：
  1. 從 Circadian 庫 circadian_gene_info 拿座標（Chromosome/Start/End/Strand）。
  2. 依 strand 定 TSS（+:Start，-:End），取啟動子窗 = TSS 上游 1000 + 下游 200。
  3. 用 pyfaidx 讀 hg38.fa 抓序列，掃：
       ebox_canonical  CACGTG 數量
       ebox_consensus  CANNTG 數量（E-box 共識，較寬鬆）
       has_ebox        有無 canonical
       gc_content      對照
  4. 併進 TCGA+GTEx 合併特徵，比「合併 vs 合併+motif」AUC（正例=CircaDB）。

輸入：hg38.fa、coexpr_features_adj.csv、gtex_features_adj.csv、CircaDB_HA_query.csv
輸出：motif_features.csv、step2h_motif_report.txt
================================================================
需要：pyfaidx, pyodbc
"""
import os
import re
import numpy as np
import pandas as pd
from pyfaidx import Fasta
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, average_precision_score

HG38 = "hg38.fa"
TCGA_F, GTEX_F, CIRCADB = "coexpr_features_adj.csv", "gtex_features_adj.csv", "CircaDB_HA_query.csv"
MOTIF_OUT = "motif_features.csv"
OUT_R = "step2h_motif_report.txt"
UP, DOWN = 1000, 200          # 啟動子窗：TSS 上游 1000、下游 200
LVL, RS, FDR = "annotation_confidence_level", 42, 0.05

report=[]
def log(m):
    print(m); report.append(m)

def load_coords():
    import pyodbc
    cs=("DRIVER={ODBC Driver 18 for SQL Server};SERVER=" + os.environ.get("CGA_SQL_SERVER", "localhost,1433") + ";"
        "DATABASE=Circadian;UID=" + os.environ.get("CGA_SQL_USER", "") + ";PWD=" + os.environ.get("CGA_SQL_PASSWORD", "") + ";TrustServerCertificate=yes;Encrypt=no")
    c=pyodbc.connect(cs); cur=c.cursor()
    cur.execute("SELECT Gene_Symbol,Chromosome,Start_Position,End_Position,Strand FROM [dbo].[circadian_gene_info]")
    rows=cur.fetchall(); c.close()
    return pd.DataFrame([tuple(r) for r in rows],
                        columns=["gene_symbol","chrom","start","end","strand"])

def extract_features():
    fa = Fasta(HG38, sequence_always_upper=True)
    contigs = set(fa.keys())
    coords = load_coords()
    log(f"[coords] 座標基因數：{len(coords)}")

    def chrom_name(x):
        x=str(x).strip()
        if x in ("MT","M"): return "chrM"
        if not x.startswith("chr"): return "chr"+x
        return x

    canon = re.compile(r"CACGTG")
    consensus = re.compile(r"CA..TG")
    rows=[]
    miss=0
    for _,r in coords.iterrows():
        ch = chrom_name(r["chrom"])
        if ch not in contigs:
            miss+=1; continue
        try:
            start=int(r["start"]); end=int(r["end"]); strand=int(r["strand"])
        except: 
            miss+=1; continue
        clen=len(fa[ch])
        if strand==1:   tss=start
        else:           tss=end
        # 窗（1-based→0-based，夾邊界）
        if strand==1:
            a,b = tss-UP, tss+DOWN
        else:
            a,b = tss-DOWN, tss+UP
        a=max(0,a-1); b=min(clen,b)
        if b<=a: 
            miss+=1; continue
        seq=str(fa[ch][a:b])
        # 反股：canonical 是回文，consensus 需兩股
        rc=seq.translate(str.maketrans("ACGT","TGCA"))[::-1]
        n_can = len(canon.findall(seq))
        n_con = len(consensus.findall(seq)) + len(consensus.findall(rc))
        gc = (seq.count("G")+seq.count("C"))/max(len(seq),1)
        rows.append({"gene_symbol":r["gene_symbol"],
                     "ebox_canonical":n_can,"has_ebox":int(n_can>0),
                     "ebox_consensus":n_con,"promoter_gc":gc})
    mf=pd.DataFrame(rows).drop_duplicates("gene_symbol")
    mf.to_csv(MOTIF_OUT,index=False,encoding="utf-8-sig")
    log(f"[motif] 產出 {len(mf)} 基因（座標對不到 contig 跳過 {miss}）")
    log(f"[motif] 有 canonical E-box 的基因比例：{mf['has_ebox'].mean()*100:.1f}%")
    return mf

def cvauc(X,y,model):
    pipe=Pipeline([("imp",SimpleImputer(strategy="median")),("sc",StandardScaler()),("clf",model)])
    cv=StratifiedKFold(5,shuffle=True,random_state=RS)
    oof=cross_val_predict(pipe,X,y,cv=cv,method="predict_proba",n_jobs=-1)[:,1]
    return roc_auc_score(y,oof),average_precision_score(y,oof)

def main():
    log("="*60); log("CGA MS3 — E-box motif 特徵測試"); log("="*60)
    mf = extract_features()

    t=pd.read_csv(TCGA_F); g=pd.read_csv(GTEX_F)
    tcga=[c for c in t.columns if c.startswith("radj_")]
    gtex=[c for c in g.columns if c.startswith("gtex_radj_")]
    motif=["ebox_canonical","has_ebox","ebox_consensus","promoter_gc"]
    df=t.merge(g[["gene_symbol"]+gtex],on="gene_symbol",how="inner").merge(mf,on="gene_symbol",how="left")

    cd=pd.read_csv(CIRCADB); circ=set(cd[cd["FDR"]<=FDR]["Gene_Symbol"].astype(str).str.strip())
    is_anchor=df["is_anchor"]==1; in_circ=df["gene_symbol"].isin(circ)
    pos=df[in_circ & ~is_anchor].copy()
    negp=df[(df[LVL]=="VERY_LOW") & ~is_anchor & ~in_circ].copy()
    neg=negp.sample(n=len(pos),random_state=RS)
    tr=pd.concat([pos,neg],ignore_index=True)
    y=np.r_[np.ones(len(pos),int),np.zeros(len(neg),int)]

    log(f"\n正例={len(pos)} 負例={len(neg)}")
    log(f"{'特徵集 / 模型':<34}{'ROC-AUC':>9}{'PR-AUC':>8}")
    log("-"*60)
    sets={"只motif":motif,"合併(TCGA+GTEx)":tcga+gtex,"合併+motif":tcga+gtex+motif}
    models={"Logit":LogisticRegression(max_iter=2000,class_weight="balanced"),
            "RF":RandomForestClassifier(n_estimators=500,random_state=RS,class_weight="balanced",n_jobs=-1)}
    best=(0,"")
    for sn,cols in sets.items():
        for mn,md in models.items():
            auc,ap=cvauc(tr[cols].to_numpy(float),y,md)
            log(f"{sn+' / '+mn:<34}{auc:>9.3f}{ap:>8.3f}")
            if auc>best[0]: best=(auc,f"{sn}/{mn}")
    log("-"*60)
    log(f"最佳：AUC={best[0]:.3f} → {best[1]}")
    log("對照：合併(TCGA+GTEx) ≈ 0.871")
    log("判讀：合併+motif 若明顯 >0.871 → E-box 有用；若≈0.871 → 序列證據無額外幫助（合理，E-box 常見）。")
    with open(OUT_R,"w",encoding="utf-8") as f: f.write("\n".join(report))
    print(f"\n[輸出] {MOTIF_OUT}、{OUT_R}")

if __name__=="__main__":
    main()