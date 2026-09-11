#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CGA MS3 — TCGA 啟動子甲基化特徵（獨立證據，無循環）
================================================================
時鐘基因失調常由表觀（甲基化）驅動。用「候選基因啟動子甲基化 vs 19 錨點
啟動子甲基化」的跨樣本相關，作為與共表達不同的獨立證據。與 GTEx/CircaDB 無關。

流程（對應 Step 1，改成甲基化）：
  1. 每基因取啟動子 probe（UCSC_RefGene_Group 含 TSS200/TSS1500/5'UTR/1stExon），
     每 sample 對這些 probe 的 Beta 值取平均 → 基因層級啟動子甲基化。
  2. 每癌內算候選 vs 19 錨點的 Spearman，跨 24 癌 Fisher-z 聚合。
  3. 背景校正（row-center + col z-score）→ meth_radj_ 特徵。

輸出：meth_features_adj.csv（gene_symbol + meth_radj_* + is_anchor + annotation_*）
================================================================
"""
import os
import numpy as np, pandas as pd
from scipy import stats

SQL_SERVER=os.environ.get("CGA_SQL_SERVER", "localhost,1433"); SQL_USER=os.environ.get("CGA_SQL_USER", ""); SQL_PW=os.environ.get("CGA_SQL_PASSWORD", "")
DB_TCGA="TCGA"; DRV="ODBC Driver 18 for SQL Server"
FEAT="coexpr_features_adj.csv"        # 拿候選清單、is_anchor、annotation_
OUT="meth_features_adj.csv"
SAMPLE_TYPE="Primary Tumor"; MIN_SAMPLES=20; GENE_BATCH=300
FDR_UNUSED=None
ANCHORS=["CLOCK","ARNTL","ARNTL2","PER1","PER2","PER3","CRY1","CRY2","NR1D1","NR1D2",
         "RORA","RORB","RORC","NPAS2","CSNK1D","CSNK1E","TIMELESS","BHLHE40","BHLHE41"]
CANCERS=["TCGA-BLCA","TCGA-BRCA","TCGA-CESC","TCGA-CHOL","TCGA-COAD","TCGA-ESCA","TCGA-GBM",
 "TCGA-HNSC","TCGA-KICH","TCGA-KIRC","TCGA-KIRP","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC","TCGA-PAAD",
 "TCGA-PCPG","TCGA-PRAD","TCGA-READ","TCGA-SARC","TCGA-SKCM","TCGA-STAD","TCGA-THCA","TCGA-THYM","TCGA-UCEC"]

# ---- 煙霧測試開關：先小規模測通（3 癌 + 500 基因），OK 再改 False 全量 ----
SMOKE_TEST=False
if SMOKE_TEST:
    CANCERS=["TCGA-LIHC","TCGA-BRCA","TCGA-COAD"]

report=[]
def log(m): print(m); report.append(m)
def conn():
    import pyodbc
    return pyodbc.connect(f"DRIVER={{{DRV}}};SERVER={SQL_SERVER};DATABASE={DB_TCGA};UID={SQL_USER};PWD={SQL_PW};TrustServerCertificate=yes;Encrypt=no")
def qdf(sql,p=None):
    with conn() as c:
        cur=c.cursor(); cur.execute(sql,p or []); cols=[d[0] for d in cur.description]; rows=cur.fetchall()
    return pd.DataFrame.from_records([tuple(r) for r in rows],columns=cols)

def pull_meth(cancer, genes):
    """
    拉單癌、Primary Tumor、指定基因的啟動子 probe Beta 值（逐 sample、逐 probe）。
    回傳 gene_name / sample_uuid / beta（已限啟動子 probe）。之後在 Python 端每基因每sample平均。
    """
    frames=[]
    promoter="(pi.[UCSC_RefGene_Group] LIKE '%TSS200%' OR pi.[UCSC_RefGene_Group] LIKE '%TSS1500%' OR pi.[UCSC_RefGene_Group] LIKE '%5''UTR%' OR pi.[UCSC_RefGene_Group] LIKE '%1stExon%')"
    for i in range(0,len(genes),GENE_BATCH):
        b=genes[i:i+GENE_BATCH]; ph=",".join("?" for _ in b)
        sql=f"""
        SELECT pgm.[GeneName] gene_name, fu.[File_UUID] su, dm.[Value] beta
        FROM [dbo].[DNA_Methylation_450K_2023_08_31] dm
        INNER JOIN [dbo].[File_UUID_DNA_Methylation_450K_2023_08_31] fu ON dm.[File_UUID]=fu.[num]
        INNER JOIN [dbo].[Probe_Num_450K_2023_08_31] pn ON dm.[Probe_ID]=pn.[num]
        INNER JOIN [dbo].[Sample_Type_DNA_Methylation_450K_2023_08_31] st ON fu.[File_UUID]=st.[File_UUID]
        INNER JOIN [dbo].[Probe_Info_450K] pi ON pn.[Probe_ID]=pi.[IlmnID]
        INNER JOIN [dbo].[Probe_Gene_Mapper_450K] pgm ON pn.[Probe_ID]=pgm.[ProbeID]
        WHERE st.[Project_ID]=? AND st.[Sample_Type]=?
          AND dm.[Value] IS NOT NULL AND dm.[Value]>=0 AND dm.[Value]<=1
          AND {promoter}
          AND pgm.[GeneName] IN ({ph})
        """
        frames.append(qdf(sql,[cancer,SAMPLE_TYPE,*b]))
    df=pd.concat(frames,ignore_index=True) if frames else pd.DataFrame(columns=["gene_name","su","beta"])
    df["beta"]=pd.to_numeric(df["beta"],errors="coerce")
    return df.dropna(subset=["beta"])

def spearman(M,acols,ccols):
    acols=[a for a in acols if a in M.columns]; ccols=[g for g in ccols if g in M.columns]
    if not acols or not ccols: return None
    sub=M[list(dict.fromkeys(acols+ccols))]; R=sub.rank(axis=0).to_numpy(float); n=R.shape[0]
    if n<MIN_SAMPLES: return None
    ix={c:i for i,c in enumerate(sub.columns)}; Rc=R-R.mean(0,keepdims=True)
    den=np.sqrt((Rc**2).sum(0)); den[den==0]=np.nan
    ai=[ix[a] for a in acols]; ci=[ix[g] for g in ccols]
    rho=np.clip((Rc[:,ci].T@Rc[:,ai])/np.outer(den[ci],den[ai]),-1,1)
    return pd.DataFrame(rho,index=ccols,columns=acols), n

def main():
    base=pd.read_csv(FEAT); uni=set(base["gene_symbol"].astype(str))
    if SMOKE_TEST:
        uni=set(sorted(uni)[:500])
    wanted=sorted(uni|set(ANCHORS))
    log(f"[universe] 候選+錨點：{len(wanted)}（SMOKE_TEST={SMOKE_TEST}）")
    rho_stack,n_stack=[],[]
    for cz in CANCERS:
        d=pull_meth(cz,wanted)
        if d.empty: log(f"  {cz}: 無甲基化資料"); continue
        # 每基因每sample對啟動子probe取平均 → 基因層級啟動子甲基化
        g=d.groupby(["su","gene_name"])["beta"].mean().reset_index()
        M=g.pivot_table(index="su",columns="gene_name",values="beta")
        ap=[a for a in ANCHORS if a in M.columns]; cp=[x for x in M.columns if x in uni and x not in ANCHORS]
        r=spearman(M,ap,cp)
        if r is None: log(f"  {cz}: 樣本/基因不足跳過"); continue
        rho_stack.append(r[0]); n_stack.append(r[1]); log(f"  {cz}: n={r[1]}, 錨點{len(r[0].columns)}, 候選{r[0].shape[0]}")
    if not rho_stack: raise SystemExit("無任何癌別產出，檢查 schema/樣本")
    log(f"[聚合] 用了 {len(rho_stack)} 癌別")
    genes=sorted(set().union(*[d.index for d in rho_stack])); anch=sorted(set().union(*[d.columns for d in rho_stack]))
    Z=np.full((len(rho_stack),len(genes),len(anch)),np.nan); W=np.zeros_like(Z)
    for k,(d,n) in enumerate(zip(rho_stack,n_stack)):
        r=np.clip(d.reindex(index=genes,columns=anch).to_numpy(float),-0.999,0.999)
        Z[k]=np.arctanh(r); W[k]=np.where(np.isfinite(Z[k]),max(n-3,0),0)
    pooled=np.tanh(np.nansum(W*Z,0)/np.where(W.sum(0)==0,np.nan,W.sum(0)))
    Rc=pooled-np.nanmean(pooled,1,keepdims=True); Rz=(Rc-np.nanmean(Rc,0))/(np.nanstd(Rc,0)+1e-9)
    out=pd.DataFrame({"gene_symbol":genes})
    for i,a in enumerate(anch): out[f"meth_radj_{a}"]=Rz[:,i]
    out=out.merge(base[["gene_symbol","is_anchor"]+[c for c in base.columns if c.startswith("annotation_")]],on="gene_symbol",how="left")
    out["is_anchor"]=out["is_anchor"].fillna(0).astype(int)
    out.to_csv(OUT,index=False,encoding="utf-8-sig")
    log(f"[out] {OUT}：{out.shape}")
    log(f"錨點命中：{len([a for a in ANCHORS if a in anch])}/19")
    print("完成。下一步把 meth_radj_ 併進 TCGA 共表達比 AUC。")
    open("meth_features_report.txt","w",encoding="utf-8").write("\n".join(report))

if __name__=="__main__":
    main()