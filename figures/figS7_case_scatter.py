#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 9/11 修改（Checklist 第 29 項）：圖上文字統一為 clock-associated prediction score、移除內部稱呼 breakthrough／clean model、Known rhythmic 改 CircaDB rhythmic；輸出改新檔名，其餘不變
"""
CGA MS3 — fig11 突破版：明星候選個案散點（組織內計算，解決跨組織混淆）
================================================================
改進：舊版將 4 組織樣本合併算相關，有跨組織混淆（Simpson's paradox）。
本版每個子圖只用「單一組織」的樣本算相關，避免組織間差異造成假相關。

從新定版候選挑明星基因，各配一個代表錨點，在最合適的單一組織內畫散點。
輸入：GTEx GCT、SampleAttributes、step_final_candidates.csv
輸出：figures/fig11_case_scatter_v3.png
================================================================
"""
import os, gzip, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

GCT="GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"
ATTR="GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
OUTDIR="figures"; PURPLE,ORANGE="#3b0764","#f59e0b"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":11,"axes.spines.top":False,
                     "axes.spines.right":False,"figure.dpi":300})
os.makedirs(OUTDIR,exist_ok=True)

# 從新定版候選挑明星基因 + 代表錨點 + 指定「單一組織」(組織內計算)
# 注意用 GTEx GCT 的 gene symbol（CCN1 部分版本為 CYR61，避開改用其他候選）
PAIRS=[("IER3","NR1D1","Liver"),
       ("IER2","BHLHE40","Lung"),
       ("BTG2","PER2","Thyroid"),
       ("THBS1","BHLHE40","Breast - Mammary Tissue")]

def main():
    cand=set(pd.read_csv("step_final_candidates.csv")["gene_symbol"])
    # 確認挑的候選在名單內（否則換）
    pairs=[p for p in PAIRS if p[0] in cand]
    if len(pairs)<4:
        # 補：從候選 Top 補齊
        extra=[g for g in pd.read_csv("step_final_candidates.csv")["gene_symbol"].head(20) if g not in [p[0] for p in pairs]]
        anchors_cycle=["NR1D1","ARNTL","PER2","BHLHE40"]; tissues_cycle=["Liver","Lung","Thyroid","Breast - Mammary Tissue"]
        i=0
        while len(pairs)<4 and i<len(extra):
            pairs.append((extra[i],anchors_cycle[len(pairs)%4],tissues_cycle[len(pairs)%4])); i+=1
    print("個案基因:", [p[0] for p in pairs])

    tissues=sorted(set(p[2] for p in pairs))
    genes_needed=set([g for p in pairs for g in (p[0],p[1])])
    attr=pd.read_csv(ATTR,sep="\t",usecols=["SAMPID","SMTSD","SMAFRZE"])
    attr=attr[(attr["SMAFRZE"]=="RNASEQ") & (attr["SMTSD"].isin(tissues))]
    samp2tis=dict(zip(attr["SAMPID"],attr["SMTSD"]))
    keep=set(attr["SAMPID"])
    with gzip.open(GCT,"rt") as f:
        f.readline();f.readline();header=f.readline().rstrip("\n").split("\t")
    cols=[c for c in header[2:] if c in keep]
    print(f"讀 GTEx（{len(cols)} 樣本）...")
    chunks=[]
    for ch in pd.read_csv(GCT,sep="\t",skiprows=2,usecols=["Description"]+cols,chunksize=20000,low_memory=False):
        ch=ch[ch["Description"].isin(genes_needed)]
        if len(ch): chunks.append(ch)
    expr=pd.concat(chunks,ignore_index=True).groupby("Description").mean(numeric_only=True)
    expr=np.log1p(expr)
    # 樣本→組織對照(只留目標組織的欄)
    tis_of={c:samp2tis.get(c) for c in expr.columns}

    fig,axes=plt.subplots(2,2,figsize=(9,8))
    for ax,(cand_g,anch,tis) in zip(axes.ravel(),pairs):
        if cand_g not in expr.index or anch not in expr.index:
            ax.set_title(f"{cand_g}/{anch}(no data)");ax.axis("off");continue
        # 只取該組織的樣本(組織內)
        cols_t=[c for c in expr.columns if tis_of.get(c)==tis]
        x=expr.loc[anch,cols_t].values.astype(float); y=expr.loc[cand_g,cols_t].values.astype(float)
        m=np.isfinite(x)&np.isfinite(y); x,y=x[m],y[m]
        if len(x)<10: ax.set_title(f"{cand_g}/{anch}({tis} n<10)");ax.axis("off");continue
        rho,p=spearmanr(x,y)
        ax.scatter(x,y,s=10,c=ORANGE,alpha=0.5,edgecolor="none")
        b,a=np.polyfit(x,y,1);xs=np.array([x.min(),x.max()])
        ax.plot(xs,a+b*xs,color=PURPLE,lw=1.5)
        ax.set_xlabel(f"{anch}  log(1+TPM)");ax.set_ylabel(f"{cand_g}  log(1+TPM)")
        ax.set_title(f"{cand_g} vs {anch} in {tis.split(' - ')[0]}\nρ={rho:.2f}, p={p:.1e}, n={len(x)}",fontsize=9)
    fig.suptitle("Candidate–anchor co-expression within single tissues",fontsize=12)
    plt.tight_layout();plt.savefig(f"{OUTDIR}/fig11_case_scatter_v3.png");plt.close()
    print("  fig11_case_scatter_v3.png\n完成。")

if __name__=="__main__":
    main()