#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 8（兩面板）取數：A 文獻評分（排除 19 個 anchor）；B 機器學習預測分數（排除 anchor、164 個訓練正例、164 個訓練負例）
輸入：step_final_scores.csv、cgdb_human_genes.txt、go_circadian_genes.txt、kegg_hsa04710_genes_20260911.txt、
      reactome_circadian_genes.txt、neg_train_seed42.txt、gold_recheck_v2_metrics.tsv（用來核對）
輸出：fig8_roc_points.tsv（ROC 座標）、fig8_auc.tsv（AUC 與 n）
"""
import os
import pandas as pd
from sklearn.metrics import roc_curve, roc_auc_score

SETS = [("CGDB", "cgdb_human_genes.txt"), ("GO", "go_circadian_genes.txt"),
        ("KEGG", "kegg_hsa04710_genes_20260911.txt"), ("Reactome", "reactome_circadian_genes.txt")]


def load_set(path):
    if not os.path.exists(path):
        raise SystemExit(f"找不到 {path}")
    g = set()
    for line in open(path, encoding="utf-8"):
        s = line.strip().split("\t")[0].strip().upper()
        if s and s != "GENE_SYMBOL":
            g.add(s)
    return g


df = pd.read_csv("step_final_scores.csv")
df["gene"] = df["gene_symbol"].astype(str).str.upper()
anchor = df["is_anchor"] == 1
train_pos = (df["is_hi_amp"] == 1) & ~anchor
neg = set(l.strip().upper() for l in open("neg_train_seed42.txt", encoding="utf-8") if l.strip())
train_neg = df["gene"].isin(neg)
if int(anchor.sum()) != 19 or int(train_pos.sum()) != 164 or int(train_neg.sum()) != 164:
    raise SystemExit(f"anchor／訓練正例／訓練負例數量不符：{int(anchor.sum())}、{int(train_pos.sum())}、{int(train_neg.sum())}")

PANELS = [("A", "annotation_combined_score", ~anchor, "文獻綜合分數", "排除 19 個 anchor"),
          ("B", "ai_score", ~(anchor | train_pos | train_neg), "AI 分數", "排除 anchor、訓練正例與訓練負例")]
ref = pd.read_csv("gold_recheck_v2_metrics.tsv", sep="\t")
pts, aucs = [], []
for panel, col, keep, ref_score, ref_variant in PANELS:
    d = df[keep]
    for name, path in SETS:
        y = d["gene"].isin(load_set(path)).astype(int)
        fpr, tpr, _ = roc_curve(y, d[col])
        a = roc_auc_score(y, d[col])
        r = ref[(ref.gene_set == name) & (ref.score == ref_score) & (ref.variant == ref_variant)]
        if len(r) != 1 or abs(r.AUC.iloc[0] - a) > 1e-6 or int(r.n_pos.iloc[0]) != int(y.sum()):
            raise SystemExit(f"與伺服器結果不符：{panel} {name} AUC {a:.6f}／n {int(y.sum())}")
        aucs.append({"panel": panel, "gene_set": name, "n_genes": len(d), "n_pos": int(y.sum()), "AUC": a})
        pts += [{"panel": panel, "gene_set": name, "fpr": f, "tpr": t} for f, t in zip(fpr, tpr)]
        print(f"  {panel} {name:9s} n = {int(y.sum()):>4}／{len(d):,}  AUC = {a:.3f}（與伺服器相符）")
pd.DataFrame(pts).to_csv("fig8_roc_points.tsv", sep="\t", index=False)
pd.DataFrame(aucs).to_csv("fig8_auc.tsv", sep="\t", index=False)
print("已輸出 fig8_roc_points.tsv、fig8_auc.tsv")
