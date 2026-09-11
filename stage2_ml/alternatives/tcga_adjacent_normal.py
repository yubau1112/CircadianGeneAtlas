#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CGA MS3 — TCGA 正常組織（Solid Tissue Normal）共表達，測是否有幫助
================================================================
用 TCGA 的 Solid Tissue Normal（癌旁正常）算共表達，背景校正後，
比 AUC：TCGA腫瘤 / TCGA正常 / GTEx正常 / 全部合併（正例=CircaDB，無循環）。

預期：TCGA 正常樣本少（癌旁、每癌數個到數十個）、且癌旁未必真「正常」，
      鑑別力可能不如 GTEx。用數據驗證。

輸入：SQL(TCGA/Circadian)、coexpr_features_adj.csv、gtex_features_adj.csv、CircaDB_HA_query.csv
輸出：tcga_normal_features_adj.csv、tcga_normal_report.txt
================================================================
"""
import os
import numpy as np, pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, average_precision_score

SQL_SERVER=os.environ.get("CGA_SQL_SERVER", "localhost,1433"); SQL_USER=os.environ.get("CGA_SQL_USER", ""); SQL_PW=os.environ.get("CGA_SQL_PASSWORD", "")
DB_TCGA="TCGA"; DB_CIRC="Circadian"; DRV="ODBC Driver 18 for SQL Server"
SAMPLE_TYPE="Solid Tissue Normal"; MIN_SAMPLES=10; MIN_FRAC=0.20; GENE_BATCH=1000
FDR, RS, LVL = 0.05, 42, "annotation_confidence_level"
ANCHORS=["CLOCK","ARNTL","ARNTL2","PER1","PER2","PER3","CRY1","CRY2","NR1D1","NR1D2",
         "RORA","RORB","RORC","NPAS2","CSNK1D","CSNK1E","TIMELESS","BHLHE40","BHLHE41"]
CANCERS=["TCGA-BLCA","TCGA-BRCA","TCGA-CESC","TCGA-CHOL","TCGA-COAD","TCGA-ESCA","TCGA-GBM",
 "TCGA-HNSC","TCGA-KICH","TCGA-KIRC","TCGA-KIRP","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC","TCGA-PAAD",
 "TCGA-PCPG","TCGA-PRAD","TCGA-READ","TCGA-SARC","TCGA-SKCM","TCGA-STAD","TCGA-THCA","TCGA-THYM","TCGA-UCEC"]

report=[]
def log(m): print(m); report.append(m)
def conn(db):
    import pyodbc
    return pyodbc.connect(f"DRIVER={{{DRV}}};SERVER={SQL_SERVER};DATABASE={db};UID={SQL_USER};PWD={SQL_PW};TrustServerCertificate=yes;Encrypt=no")
def qdf(db,sql,p=None):
    with conn(db) as c:
        cur=c.cursor(); cur.execute(sql,p or []); cols=[d[0] for d in cur.description]; rows=cur.fetchall()
    return pd.DataFrame.from_records([tuple(r) for r in rows], columns=cols)

def pull(cancer, genes):
    frames=[]
    for i in range(0,len(genes),GENE_BATCH):
        b=genes[i:i+GENE_BATCH]; ph=",".join("?" for _ in b)
        sql=f"""SELECT C.[gene_name] gene_name,B.[File_UUID] su,A.[tpm_unstranded] tpm
        FROM [{DB_TCGA}].[dbo].[Gene_Expression] A
        JOIN [{DB_TCGA}].[dbo].[File_UUID_Gene_Expression] B ON A.[File_UUID]=B.[num]
        JOIN [{DB_TCGA}].[dbo].[Gene_Info_GENCODEv36] C ON A.[gene_id]=C.[num]
        JOIN [{DB_TCGA}].[dbo].[Sample_Type_Gene_Expression] D ON B.[File_UUID]=D.[File_UUID]
        WHERE D.[Project_ID]=? AND D.[Sample_Type]=? AND A.[tpm_unstranded] IS NOT NULL
          AND A.[tpm_unstranded]>=0 AND C.[gene_name] IN ({ph})"""
        frames.append(qdf(DB_TCGA,sql,[cancer,SAMPLE_TYPE,*b]))
    df=pd.concat(frames,ignore_index=True) if frames else pd.DataFrame(columns=["gene_name","su","tpm"])
    df["tpm"]=pd.to_numeric(df["tpm"],errors="coerce"); return df.dropna(subset=["tpm"])

def spearman(M,acols,ccols):
    frac=(np.expm1(M)>0).mean(); keep=frac[frac>=MIN_FRAC].index
    acols=[a for a in acols if a in keep]; ccols=[g for g in ccols if g in keep]
    if not acols or not ccols: return None
    sub=M[list(dict.fromkeys(acols+ccols))]; R=sub.rank(axis=0).to_numpy(float); n=R.shape[0]
    if n<MIN_SAMPLES: return None
    ix={c:i for i,c in enumerate(sub.columns)}; Rc=R-R.mean(0,keepdims=True)
    den=np.sqrt((Rc**2).sum(0)); den[den==0]=np.nan
    ai=[ix[a] for a in acols]; ci=[ix[g] for g in ccols]
    rho=np.clip((Rc[:,ci].T@Rc[:,ai])/np.outer(den[ci],den[ai]),-1,1)
    return pd.DataFrame(rho,index=ccols,columns=acols), n

def main():
    base=pd.read_csv("coexpr_features_adj.csv"); uni=set(base["gene_symbol"].astype(str))
    wanted=sorted(uni|set(ANCHORS))
    log(f"[universe] {len(wanted)}")
    rho_stack,n_stack=[],[]
    for cz in CANCERS:
        d=pull(cz,wanted)
        if d.empty: log(f"  {cz}: 無正常樣本"); continue
        M=np.log1p(d.pivot_table(index="su",columns="gene_name",values="tpm",aggfunc="mean"))
        ap=[a for a in ANCHORS if a in M.columns]; cp=[g for g in M.columns if g in uni and g not in ANCHORS]
        r=spearman(M,ap,cp)
        if r is None: log(f"  {cz}: 樣本不足跳過"); continue
        rho_stack.append(r[0]); n_stack.append(r[1]); log(f"  {cz}: n={r[1]}")
    if not rho_stack: raise SystemExit("無任何癌別有足夠正常樣本")
    log(f"[聚合] 用了 {len(rho_stack)} 個癌別")
    genes=sorted(set().union(*[d.index for d in rho_stack])); anch=sorted(set().union(*[d.columns for d in rho_stack]))
    Z=np.full((len(rho_stack),len(genes),len(anch)),np.nan); W=np.zeros_like(Z)
    for k,(d,n) in enumerate(zip(rho_stack,n_stack)):
        r=np.clip(d.reindex(index=genes,columns=anch).to_numpy(float),-0.999,0.999)
        Z[k]=np.arctanh(r); W[k]=np.where(np.isfinite(Z[k]),max(n-3,0),0)
    pooled=np.tanh(np.nansum(W*Z,0)/np.where(W.sum(0)==0,np.nan,W.sum(0)))
    Rc=pooled-np.nanmean(pooled,1,keepdims=True); Rz=(Rc-np.nanmean(Rc,0))/(np.nanstd(Rc,0)+1e-9)
    out=pd.DataFrame({"gene_symbol":genes})
    for i,a in enumerate(anch): out[f"tnorm_radj_{a}"]=Rz[:,i]
    out=out.merge(base[["gene_symbol","is_anchor"]+[c for c in base.columns if c.startswith("annotation_")]],on="gene_symbol",how="left")
    out["is_anchor"]=out["is_anchor"].fillna(0).astype(int)
    out.to_csv("tcga_normal_features_adj.csv",index=False,encoding="utf-8-sig")
    log(f"[out] tcga_normal_features_adj.csv：{out.shape}")

    # ==== 比 AUC ====
    g=pd.read_csv("gtex_features_adj.csv")
    tcga=[c for c in base.columns if c.startswith("radj_")]
    gtex=[c for c in g.columns if c.startswith("gtex_radj_")]
    tnorm=[c for c in out.columns if c.startswith("tnorm_radj_")]
    df=base.merge(g[["gene_symbol"]+gtex],on="gene_symbol",how="inner").merge(out[["gene_symbol"]+tnorm],on="gene_symbol",how="left")
    cd=pd.read_csv("CircaDB_HA_query.csv"); circ=set(cd[cd["FDR"]<=FDR]["Gene_Symbol"].astype(str).str.strip())
    ia=df["is_anchor"]==1; ic=df["gene_symbol"].isin(circ)
    pos=df[ic&~ia].copy(); neg=df[(df[LVL]=="VERY_LOW")&~ia&~ic].sample(n=len(pos),random_state=RS)
    tr=pd.concat([pos,neg],ignore_index=True); y=np.r_[np.ones(len(pos),int),np.zeros(len(neg),int)]
    def auc(cols,mdl):
        pipe=Pipeline([("i",SimpleImputer(strategy="median")),("s",StandardScaler()),("m",mdl)])
        cv=StratifiedKFold(5,shuffle=True,random_state=RS)
        o=cross_val_predict(pipe,tr[cols].to_numpy(float),y,cv=cv,method="predict_proba",n_jobs=-1)[:,1]
        return roc_auc_score(y,o),average_precision_score(y,o)
    log("\n"+"="*56); log("AUC 比較（正例=CircaDB，無循環）"); log("="*56)
    log(f"{'訊號源 / 模型':<34}{'ROC-AUC':>9}{'PR-AUC':>8}"); log("-"*56)
    sets={"TCGA 腫瘤":tcga,"TCGA 正常(癌旁)":tnorm,"GTEx 正常":gtex,
          "TCGA腫瘤+GTEx(現定版)":tcga+gtex,"全部合併":tcga+tnorm+gtex}
    for sn,cols in sets.items():
        for mn,md in [("Logit",LogisticRegression(max_iter=2000,class_weight="balanced")),
                      ("RF",RandomForestClassifier(n_estimators=500,random_state=RS,class_weight="balanced",n_jobs=-1))]:
            a,p=auc(cols,md); log(f"{sn+' / '+mn:<34}{a:>9.3f}{p:>8.3f}")
    log("-"*56); log("對照：現定版 TCGA腫瘤+GTEx = 0.871")
    open("tcga_normal_report.txt","w",encoding="utf-8").write("\n".join(report))
    print("\n[輸出] tcga_normal_report.txt")

if __name__=="__main__":
    main()