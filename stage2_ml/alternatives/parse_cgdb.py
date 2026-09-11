#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CGA MS3 — 解析 CGDB TableS1 → 人類實驗節律基因正例清單
================================================================
CGDB 是「實驗驗證會震盪」的節律基因（RT-PCR/Northern blot 等），與 CircaDB 同性質
（實驗震盪），故可能與共表達特徵對得上。且與 GTEx 完全獨立（無循環）。

篩人類 → 抽 UniProt/Ensembl ID → mygene 轉 gene symbol → 存清單。
輸出：cgdb_human_genes.txt（每行一個 gene symbol）
================================================================
需要：openpyxl, mygene
"""
import pandas as pd, mygene

XLSX="TableS1.xlsx"; OUT="cgdb_human_genes.txt"

def main():
    df=pd.read_excel(XLSX, sheet_name=0, header=2)
    df.columns=[str(c).strip() for c in df.columns]
    org=[c for c in df.columns if "Organism" in c][0]
    idc=[c for c in df.columns if "Ensembl" in c or "Uniprot" in c][0]
    hs=df[df[org].astype(str).str.contains("Homo sapiens", na=False)]
    ids=sorted(set(hs[idc].dropna().astype(str).str.strip()))
    print(f"人類唯一 ID：{len(ids)}")

    # ID 可能是 UniProt(P/Q/O 開頭) 或 Ensembl(ENSG/ENSP)。分開查。
    uniprot=[i for i in ids if not i.upper().startswith("ENS")]
    ensembl=[i for i in ids if i.upper().startswith("ENS")]
    print(f"  UniProt：{len(uniprot)}  Ensembl：{len(ensembl)}")

    mg=mygene.MyGeneInfo()
    syms=set()
    if uniprot:
        for r in mg.querymany(uniprot, scopes="uniprot", fields="symbol", species="human", verbose=False):
            if r.get("symbol"): syms.add(r["symbol"])
    if ensembl:
        for r in mg.querymany(ensembl, scopes="ensembl.gene,ensembl.protein", fields="symbol", species="human", verbose=False):
            if r.get("symbol"): syms.add(r["symbol"])
    syms=sorted(syms)
    print(f"轉出 gene symbol 數：{len(syms)}")

    with open(OUT,"w",encoding="utf-8") as f:
        for s in syms: f.write(s+"\n")
    print(f"[輸出] {OUT}")
    core=["CLOCK","ARNTL","PER1","PER2","PER3","CRY1","CRY2","NR1D1","NR1D2","RORA","BHLHE40","NPAS2","TIMELESS","DBP","TEF","HLF"]
    hit=[g for g in core if g in syms]
    print(f"核心/經典節律基因命中：{len(hit)} → {hit}")
    print("前 30 個:", syms[:30])

if __name__=="__main__":
    main()