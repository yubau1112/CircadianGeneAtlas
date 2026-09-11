#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CGA MS3 — 第二支:Reactome(via Enrichr) + GO/Pathway enrichment + Network
================================================================
老師建議後續:
1. 補 Reactome gold standard(透過 Enrichr 抓,繞過被擋的 reactome.org)
2. 399 候選 + Top20 高分基因群的 GO/KEGG/Reactome enrichment
3. Network 圖(Top20 高分基因)
全走 Enrichr(maayanlab.cloud,已確認可連)。不編造,抓不到就標註。
================================================================
"""
import os, numpy as np, pandas as pd, warnings, urllib.request, json
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
OUTDIR="figures"; os.makedirs(OUTDIR,exist_ok=True)
PURPLE,ORANGE,BLUE,GREEN="#3b0764","#f59e0b","#2563eb","#059669"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.spines.top":False,"axes.spines.right":False,"figure.dpi":300})

CAND="step_final_candidates.csv"; SCORES="step_final_scores.csv"

def enrichr_geneset(library):
    """從 Enrichr 抓某個 library 的基因集(繞過被擋的原站)"""
    url=f"https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={library}"
    hdr={"User-Agent":"Mozilla/5.0"}
    req=urllib.request.Request(url,headers=hdr)
    txt=urllib.request.urlopen(req,timeout=60).read().decode()
    terms={}
    for line in txt.strip().split("\n"):
        parts=line.split("\t")
        term=parts[0]; genes=[g.split(",")[0].strip().upper() for g in parts[2:] if g.strip()]
        terms[term]=set(genes)
    return terms

def run_enrichr(gene_list, library, top=15):
    """用 gseapy 對 gene_list 做 enrichment"""
    import gseapy as gp
    try:
        enr=gp.enrichr(gene_list=list(gene_list),gene_sets=library,outdir=None,no_plot=True)
        df=enr.results.sort_values("Adjusted P-value").head(top)
        return df
    except Exception as e:
        print(f"  enrichr {library} 失敗: {str(e)[:80]}"); return None

def main():
    print("="*60); print("CGA MS3 — Enrichment + Network (老師建議)"); print("="*60)
    cand=pd.read_csv(CAND)
    gcol=[c for c in cand.columns if 'symbol' in c.lower() or c.lower()=='gene'][0]
    cand_genes=[str(g).upper() for g in cand[gcol].dropna()]
    print(f"候選基因數: {len(cand_genes)}")
    # Top20 高分
    scores=pd.read_csv(SCORES)
    sgc=[c for c in scores.columns if 'symbol' in c.lower() or c.lower()=='gene'][0]
    ssc=[c for c in scores.columns if 'ai_score' in c.lower() or 'score' in c.lower()][0]
    cand_scores=scores[scores[sgc].astype(str).str.upper().isin(cand_genes)].nlargest(20,ssc)
    top20=[str(g).upper() for g in cand_scores[sgc]]
    print(f"Top20 高分候選: {', '.join(top20[:10])}...")

    # ── 1. Reactome gold standard (via Enrichr) ──
    print("\n[1] 抓 Reactome circadian 基因集(via Enrichr):")
    try:
        reactome_lib=enrichr_geneset("Reactome_2022")
        circ_terms={k:v for k,v in reactome_lib.items() if 'circadian' in k.lower() or 'clock' in k.lower()}
        print(f"  Reactome circadian 相關 terms: {list(circ_terms.keys())}")
        reactome_genes=set()
        for v in circ_terms.values(): reactome_genes|=v
        print(f"  Reactome circadian 基因數: {len(reactome_genes)}")
        # 存起來給 gold standard 用
        if reactome_genes:
            pd.Series(sorted(reactome_genes)).to_csv("reactome_circadian_genes.txt",index=False,header=False)
            print("  已存 reactome_circadian_genes.txt")
    except Exception as e:
        print(f"  Reactome via Enrichr 失敗: {str(e)[:80]}")

    # ── 2. Enrichment: 399候選 對 GO/KEGG/Reactome ──
    print("\n[2] 399 候選的 enrichment:")
    libs={"GO_Biological_Process_2021":"GO-BP","KEGG_2021_Human":"KEGG","Reactome_2022":"Reactome"}
    all_enr={}
    for lib,short in libs.items():
        print(f"  跑 {short}...")
        df=run_enrichr(cand_genes,lib,top=15)
        if df is not None and len(df):
            all_enr[short]=df
            df.to_csv(f"enrichment_{short}_candidates.csv",index=False)
            # 印前5
            for _,r in df.head(5).iterrows():
                print(f"    {r['Term'][:50]:50s} p.adj={r['Adjusted P-value']:.2e}")

    # ── 3. Enrichment 長條圖(GO-BP前10) ──
    if "GO-BP" in all_enr:
        import textwrap, re as _re
        df=all_enr["GO-BP"].head(10).iloc[::-1]
        # 第一行:完整 term 名稱(太長換行);最後一行:GO ID
        def fmt_term(t):
            m=_re.search(r"\((GO:\d+)\)",t)
            goid=m.group(1) if m else ""
            name=_re.sub(r"\s*\(GO:\d+\)","",t).strip()
            wrapped="\n".join(textwrap.wrap(name,53))
            return wrapped+("\n("+goid+")" if goid else "")
        labels=[fmt_term(t) for t in df["Term"]]
        fig,ax=plt.subplots(figsize=(13,6.5))
        y=range(len(df)); vals=-np.log10(df["Adjusted P-value"])
        ax.barh(list(y),vals,color=PURPLE)
        ax.set_yticks(list(y)); ax.set_yticklabels(labels,fontsize=9)
        ax.set_xlabel("-log10(adjusted P-value)")
        ax.set_title("GO Biological Process enrichment (399 candidates)")
        plt.subplots_adjust(left=0.40,right=0.97,top=0.90,bottom=0.12)  # 左邊大留白給長標籤
        plt.savefig(f"{OUTDIR}/enrichment_GO_BP.png"); plt.close()
        print(f"\n  圖: {OUTDIR}/enrichment_GO_BP.png")

    # ── 4. Network 圖(Top20 高分基因,用 STRING) ──
    print("\n[3] Top20 network(用已下載的 STRING):")
    try:
        import networkx as nx
        # 用已有的 STRING 檔
        info=pd.read_csv("9606.protein.info.v12.0.txt.gz",sep="\t")
        links=pd.read_csv("9606.protein.links.v12.0.txt.gz",sep=" ")
        # symbol -> string id
        info["pref"]=info["preferred_name"].str.upper()
        sym2id=dict(zip(info["pref"],info["#string_protein_id"]))
        top20_ids={sym2id[g] for g in top20 if g in sym2id}
        sub=links[(links["protein1"].isin(top20_ids))&(links["protein2"].isin(top20_ids))&(links["combined_score"]>=400)]
        id2sym={v:k for k,v in sym2id.items()}
        G=nx.Graph()
        for g in top20: G.add_node(g)
        for _,r in sub.iterrows():
            a=id2sym.get(r["protein1"]); b=id2sym.get(r["protein2"])
            if a and b: G.add_edge(a,b,weight=r["combined_score"])
        # 分「有連結」與「孤立」兩群分開排,避免擠成一團
        connected=[n for n in G.nodes() if G.degree(n)>0]
        isolated=[n for n in G.nodes() if G.degree(n)==0]
        Gc=G.subgraph(connected)
        fig,ax=plt.subplots(figsize=(11,8))
        import numpy as _np
        # 有連結的:用 kamada_kawai(比 spring 更展開)放在左側主區
        if len(connected)>0:
            pos_c=nx.kamada_kawai_layout(Gc)
            # 縮放到左 2/3 區域
            xs=_np.array([pos_c[n][0] for n in connected]); ys=_np.array([pos_c[n][1] for n in connected])
            xr=xs.max()-xs.min() or 1; yr=ys.max()-ys.min() or 1
            pos={n:((pos_c[n][0]-xs.min())/xr*0.6+0.03,(pos_c[n][1]-ys.min())/yr*0.9+0.05) for n in connected}
        else: pos={}
        # 孤立的:右側直欄排列
        for i,n in enumerate(isolated):
            col=i//7; row=i%7
            pos[n]=(0.72+col*0.14, 0.9-row*0.13)
        nx.draw_networkx_nodes(G,pos,nodelist=connected,node_color=ORANGE,node_size=900,edgecolors=PURPLE,linewidths=2,ax=ax)
        nx.draw_networkx_nodes(G,pos,nodelist=isolated,node_color="#e5e7eb",node_size=650,edgecolors="#9ca3af",linewidths=1.5,ax=ax)
        nx.draw_networkx_edges(G,pos,edge_color=PURPLE,width=1.6,alpha=0.5,ax=ax)
        # 標籤:連結的放節點上方避免疊在點內,孤立的放右邊
        for n in connected:
            ax.text(pos[n][0],pos[n][1]+0.035,n,fontsize=8,fontweight="bold",ha="center",va="bottom",color="#1f2937")
        for n in isolated:
            ax.text(pos[n][0]+0.035,pos[n][1],n,fontsize=7.5,ha="left",va="center",color="#6b7280")
        ax.set_title(f"Top20 candidate network (STRING, score>=400)\n{Gc.number_of_edges()} edges among {len(connected)} connected genes; {len(isolated)} singletons (grey)",fontsize=11)
        ax.set_xlim(-0.02,1.05); ax.set_ylim(-0.02,1.02); ax.axis("off")
        plt.tight_layout(); plt.savefig(f"{OUTDIR}/top20_network.png"); plt.close()
        print(f"  圖: {OUTDIR}/top20_network.png ({Gc.number_of_edges()} edges, {len(connected)} connected, {len(isolated)} isolated)")
    except Exception as e:
        print(f"  Network 失敗: {str(e)[:100]}")

    print("\n完成。enrichment 表 + 圖已輸出。")

if __name__=="__main__":
    main()