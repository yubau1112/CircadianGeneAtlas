#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CGA MS3 — 解析 GO circadian GAF → gene symbol 正例清單
================================================================
讀 AmiGO 下載的 GO circadian 註釋（GAF 格式，人類），抽 UniProt ID，
用 mygene 轉成 gene symbol，存成「非 GTEx 來源」的節律正例清單。

可選：只留 GO:0007623（核心 circadian rhythm）或含更廣的節律 GO term。
輸出：go_circadian_genes.txt（每行一個 gene symbol）
================================================================
需要：mygene（已裝）
"""
import mygene

GAF="go_circadian.gaf"
OUT="go_circadian_genes.txt"
# 要納入的 GO term：預設納入全部節律相關（正例更充足）。
# 若只要最嚴格的核心節律，改成 {"GO:0007623"}。
KEEP_GO={"GO:0007623","GO:0032922","GO:0048512","GO:0045475","GO:0060086","GO:0097167","GO:0050802","GO:0042745","GO:0022410"}

def main():
    uni=set()
    with open(GAF,encoding="utf-8") as f:
        for line in f:
            if line.startswith("!"): continue
            p=line.rstrip("\n\r").split("\t")
            if len(p)>=5 and p[0].startswith("UniProtKB:") and p[3] in KEEP_GO:
                uni.add(p[0].replace("UniProtKB:",""))
    print(f"唯一 UniProt ID（符合 GO 條件）：{len(uni)}")

    mg=mygene.MyGeneInfo()
    res=mg.querymany(list(uni),scopes="uniprot",fields="symbol",species="human",verbose=False)
    syms=set()
    for r in res:
        if "symbol" in r and r["symbol"]:
            syms.add(r["symbol"])
    syms=sorted(syms)
    print(f"轉出 gene symbol 數：{len(syms)}")

    with open(OUT,"w",encoding="utf-8") as f:
        for s in syms: f.write(s+"\n")
    print(f"[輸出] {OUT}")
    # 抽查：核心時鐘基因在不在
    core=["CLOCK","ARNTL","PER1","PER2","PER3","CRY1","CRY2","NR1D1","NR1D2","RORA","BHLHE40","NPAS2","TIMELESS"]
    hit=[g for g in core if g in syms]
    print(f"核心時鐘基因命中：{len(hit)}/{len(core)} → {hit}")
    print("前 30 個 symbol:", syms[:30])

if __name__=="__main__":
    main()