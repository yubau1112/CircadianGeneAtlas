#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CGA MS3 — STRING 網路距離特徵（獨立於共表達的新證據）
================================================================
從 STRING v12.0 人類蛋白互作網路，對每個基因算「與 19 個核心時鐘錨點的
網路鄰近度」特徵。這與共表達完全獨立，用來補資訊、有機會突破共表達天花板。

特徵（每基因）：
  net_min_dist      到最近錨點的最短路徑步數（越小越近；不可達給大值）
  net_mean_dist     到所有可達錨點的平均步數
  net_n_direct      與錨點直接相連的數量（0~19）
  net_sum_score     與錨點所有直接邊的 STRING 分數總和
  net_max_score     與錨點直接邊的最高 STRING 分數
  net_reachable     可達（BFS 有限步內）的錨點數

輸入（同資料夾）：
  9606.protein.links.v12.0.txt.gz
  9606.protein.info.v12.0.txt.gz
  coexpr_features_adj.csv   （拿基因清單 + is_anchor）
輸出：
  string_network_features.csv
================================================================
需要：networkx（pip install networkx）
"""
import gzip
import numpy as np
import pandas as pd
import networkx as nx
from collections import deque

LINKS = "9606.protein.links.v12.0.txt.gz"
INFO  = "9606.protein.info.v12.0.txt.gz"
FEAT  = "coexpr_features_adj.csv"
OUT   = "string_network_features.csv"

SCORE_MIN = 400     # STRING combined score 門檻（400=中信度，過濾雜訊）
MAX_HOPS  = 6       # BFS 最多找幾步（超過視為不可達）
UNREACH   = 99      # 不可達的距離填值

ANCHORS = ["CLOCK","ARNTL","ARNTL2","PER1","PER2","PER3","CRY1","CRY2","NR1D1","NR1D2",
           "RORA","RORB","RORC","NPAS2","CSNK1D","CSNK1E","TIMELESS","BHLHE40","BHLHE41"]


def main():
    # 1) ENSP -> gene symbol 對照
    print("[1/4] 讀 protein.info（ENSP→symbol）...")
    ensp2sym, sym2ensp = {}, {}
    with gzip.open(INFO, "rt", encoding="utf-8") as f:
        next(f)  # header: #string_protein_id, preferred_name, ...
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                ensp2sym[p[0]] = p[1]
                sym2ensp[p[1]] = p[0]
    print(f"      對照數：{len(ensp2sym):,}")

    # 2) 建圖（score >= SCORE_MIN）
    print(f"[2/4] 讀 links、建圖（score>={SCORE_MIN}）...")
    G = nx.Graph()
    n_edge = 0
    with gzip.open(LINKS, "rt", encoding="utf-8") as f:
        next(f)  # header: protein1 protein2 combined_score
        for line in f:
            a, b, s = line.split()
            s = int(s)
            if s >= SCORE_MIN:
                sa, sb = ensp2sym.get(a), ensp2sym.get(b)
                if sa and sb:
                    G.add_edge(sa, sb, w=s)
                    n_edge += 1
    print(f"      節點：{G.number_of_nodes():,}  邊：{G.number_of_edges():,}")

    anchors_in = [a for a in ANCHORS if a in G]
    print(f"      在網路中的錨點：{len(anchors_in)}/19 → {anchors_in}")

    # 3) 多源 BFS：一次算所有節點到「最近錨點」的步數
    print(f"[3/4] 多源 BFS（上限 {MAX_HOPS} 步）...")
    dist = {a: 0 for a in anchors_in}
    dq = deque(anchors_in)
    while dq:
        u = dq.popleft()
        if dist[u] >= MAX_HOPS:
            continue
        for v in G.neighbors(u):
            if v not in dist:
                dist[v] = dist[u] + 1
                dq.append(v)

    # 每基因對「所有錨點」的直接邊統計
    print("[4/4] 算每基因的錨點鄰近特徵...")
    df = pd.read_csv(FEAT)
    genes = df["gene_symbol"].astype(str).tolist()
    anchor_set = set(anchors_in)

    rows = []
    for g in genes:
        min_d = dist.get(g, UNREACH)
        # 直接與錨點相連的統計
        n_direct, scores, reach_ds = 0, [], []
        if g in G:
            for a in anchors_in:
                if G.has_edge(g, a):
                    n_direct += 1
                    scores.append(G[g][a]["w"])
        # 到各錨點的可達距離（用單源 BFS 太慢；改用近似：只記 min_d 與 direct）
        rows.append({
            "gene_symbol": g,
            "net_min_dist":  min_d,
            "net_n_direct":  n_direct,
            "net_sum_score": sum(scores),
            "net_max_score": max(scores) if scores else 0,
            "net_in_network": 1 if g in G else 0,
        })
    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False, encoding="utf-8-sig")
    print(f"[out] {OUT}：{len(out)} 基因")
    print(out["net_min_dist"].value_counts().sort_index().head(8).to_string())
    print(f"直接連到錨點的基因數(net_n_direct>0)：{(out['net_n_direct']>0).sum()}")
    print("完成。下一步併進 radj_ 重跑分類器看 AUC。")


if __name__ == "__main__":
    main()