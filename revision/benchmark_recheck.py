#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CGA 第 17 項：重算 Figure 8（gold standard ROC）v2：再加一版排除訓練負例（neg_train_seed42.txt）
原腳本 cga_ms3_goldstandard_eval_2.py 用 step_final_scores.csv 的 ai_score，且沒有排除任何基因。
這支只算數字、不畫圖；畫圖另外做。

分數 × 排除方式：
  AI 分數：不排除（重現原數字）／排除 19 個 anchor／排除 anchor 與 164 個訓練正例
  文獻綜合分數：不排除／排除 anchor（anchor 有已知基因加分）
執行位置：/home/yubau/aging/20250519_circadian/Phase2/circadian_pipeline
輸出：gold_recheck_v2_metrics.tsv、gold_recheck_v2_overlap.tsv；KEGG 讀 kegg_hsa04710_genes_20260911.txt（上午存的），確保與上午同一份清單
"""
import os, datetime, urllib.request
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

SCORES = "step_final_scores.csv"
SETS_LOCAL = {"CGDB": "cgdb_human_genes.txt", "GO": "go_circadian_genes.txt", "Reactome": "reactome_circadian_genes.txt"}
ORIG = {"CGDB": (909, 0.530), "GO": (127, 0.575), "KEGG": (33, 0.695), "Reactome": (63, 0.635)}  # 投稿版 Figure 8
TODAY = datetime.date.today().strftime("%Y%m%d")


def load_set(path):
    if not os.path.exists(path):
        raise SystemExit(f"找不到基因集檔案：{path}")
    genes = set()
    for line in open(path, encoding="utf-8"):
        g = line.strip().split("\t")[0].strip().upper()
        if g and g != "GENE_SYMBOL":
            genes.add(g)
    return genes


def fetch_kegg():
    """與原腳本相同的抓法；抓到後存檔以便追溯"""
    txt = urllib.request.urlopen("https://rest.kegg.jp/get/hsa04710", timeout=30).read().decode()
    genes, in_gene = set(), False
    for line in txt.split("\n"):
        if line.startswith("GENE"):
            in_gene = True
        elif line and not line.startswith(" ") and in_gene:
            in_gene = False
        if in_gene:
            for p in line.split():
                if ";" in p:
                    genes.add(p.replace(";", "").strip().upper())
                    break
    if not genes:
        raise SystemExit("KEGG 抓不到基因")
    out = f"kegg_hsa04710_genes_{TODAY}.txt"
    with open(out, "w", encoding="utf-8") as f:
        f.write("\n".join(sorted(genes)) + "\n")
    print(f"KEGG 清單已存：{out}（{len(genes)} 個）")
    return genes


df = pd.read_csv(SCORES)
need = ["gene_symbol", "ai_score", "is_anchor", "is_circadb", "is_hi_amp", "annotation_combined_score"]
miss = [c for c in need if c not in df.columns]
if miss:
    raise SystemExit(f"缺欄位：{miss}")
df["gene"] = df["gene_symbol"].astype(str).str.upper()
if df["gene"].duplicated().any():
    raise SystemExit(f"基因符號重複 {int(df['gene'].duplicated().sum())} 個，先停止")
anchor = df["is_anchor"] == 1
train_pos = (df["is_hi_amp"] == 1) & (~anchor)
NEG_FILE = "neg_train_seed42.txt"
if not os.path.exists(NEG_FILE):
    raise SystemExit(f"找不到 {NEG_FILE}，請先跑 cga_neg_sensitivity.py")
neg_set = set(l.strip().upper() for l in open(NEG_FILE, encoding="utf-8") if l.strip())
train_neg = df["gene"].isin(neg_set)
print(f"訓練負例：檔案 {len(neg_set)} 個，對到分數表 {int(train_neg.sum())} 個")
if int(train_neg.sum()) != 164:
    raise SystemExit("訓練負例不是 164 個，先停止")
print(f"基因數 {len(df):,}；anchor {int(anchor.sum())}；訓練正例（is_hi_amp=1 且非 anchor）{int(train_pos.sum())}")
if int(anchor.sum()) != 19:
    raise SystemExit("anchor 不是 19 個，先停止")
if int(train_pos.sum()) != 164:
    print("注意：訓練正例不是 164 個，請把這行貼給我再判斷")

sets = {k: load_set(v) for k, v in SETS_LOCAL.items()}
KEGG_FILE = "kegg_hsa04710_genes_20260911.txt"
sets["KEGG"] = load_set(KEGG_FILE) if os.path.exists(KEGG_FILE) else fetch_kegg()
print(f"KEGG：{len(sets['KEGG'])} 個（{KEGG_FILE if os.path.exists(KEGG_FILE) else '重新抓取'}）")
sets = {k: sets[k] for k in ["CGDB", "GO", "KEGG", "Reactome"]}

VARIANTS = [
    ("AI 分數", "ai_score", "不排除（重現原數字）", pd.Series(False, index=df.index)),
    ("AI 分數", "ai_score", "排除 19 個 anchor", anchor),
    ("AI 分數", "ai_score", "排除 anchor 與 164 個訓練正例", anchor | train_pos),
    ("AI 分數", "ai_score", "排除 anchor、訓練正例與訓練負例", anchor | train_pos | train_neg),
    ("文獻綜合分數", "annotation_combined_score", "不排除", pd.Series(False, index=df.index)),
    ("文獻綜合分數", "annotation_combined_score", "排除 19 個 anchor", anchor),
]

rows, overlap = [], []
for name, genes in sets.items():
    lab = df["gene"].isin(genes)
    overlap.append({"gene_set": name, "set_size": len(genes), "in_scored_genes": int(lab.sum()),
                    "overlap_anchor": int((lab & anchor).sum()), "overlap_train_pos": int((lab & train_pos).sum()),
                    "overlap_circadb": int((lab & (df["is_circadb"] == 1)).sum()), "overlap_train_neg": int((lab & train_neg).sum())})
    for score_name, col, vname, excl in VARIANTS:
        d = df[~excl]
        y = d["gene"].isin(genes).astype(int)
        n_pos = int(y.sum())
        a = roc_auc_score(y, d[col]) if n_pos >= 5 else np.nan
        rows.append({"gene_set": name, "score": score_name, "variant": vname, "n_genes": len(d), "n_pos": n_pos, "AUC": a})

res = pd.DataFrame(rows)
ov = pd.DataFrame(overlap)
res.to_csv("gold_recheck_v2_metrics.tsv", sep="\t", index=False)
ov.to_csv("gold_recheck_v2_overlap.tsv", sep="\t", index=False)

print("\n=== 各參考集與 anchor、訓練正例的重疊 ===")
print(ov.to_string(index=False))

print("\n=== 重現檢查（AI 分數、不排除）vs 投稿版 Figure 8 ===")
for name in sets:
    r = res[(res.gene_set == name) & (res.score == "AI 分數") & (res.variant.str.startswith("不排除"))].iloc[0]
    n0, a0 = ORIG[name]
    flag = "相符" if (r.n_pos == n0 and abs(r.AUC - a0) < 0.0015) else "不同，請留意"
    print(f"  {name:9s} n {r.n_pos}（原 {n0}）  AUC {r.AUC:.3f}（原 {a0:.3f}）  {flag}")

print("\n=== AUC 一覽 ===")
tab = res.pivot_table(index=["score", "variant"], columns="gene_set", values="AUC", sort=False)[list(sets)]
npos = res.pivot_table(index=["score", "variant"], columns="gene_set", values="n_pos", sort=False)[list(sets)]
for idx in tab.index:
    cells = "  ".join(f"{g} {tab.loc[idx, g]:.3f}（n={int(npos.loc[idx, g])}）" for g in sets)
    print(f"  {idx[0]}｜{idx[1]}：{cells}")
print("\n已輸出 gold_recheck_v2_metrics.tsv、gold_recheck_v2_overlap.tsv")