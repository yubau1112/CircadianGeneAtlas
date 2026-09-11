#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
查核 Figure S3（背景校正）與 2.8「cosine 0.140 vs 0.008, P = 5.8e-10」：
  1. 用與 cga_ms3_make_figures.py fig5_bgcorrect 相同的做法重現校正前、後的 p 值
  2. 印出 HIGH、VERY_LOW 兩群的平均與中位數相似度
  3. 留一法：每個 HIGH 基因改與「不含自己的 HIGH 平均模板」比較，量化自我比較造成的偏誤
只讀檔，不寫任何檔案。執行位置：Phase2/circadian_pipeline
v2：指定正式分析的特徵檔（v1 依字母序誤取 399 個候選的檔案）
"""
import glob
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu

LVL = "annotation_confidence_level"
ANCHORS = {"CLOCK", "ARNTL", "ARNTL2", "PER1", "PER2", "PER3", "CRY1", "CRY2", "NR1D1", "NR1D2",
           "RORA", "RORB", "RORC", "NPAS2", "CSNK1D", "CSNK1E", "TIMELESS", "BHLHE40", "BHLHE41"}


def find_table(prefix):
    """找出含 gene_symbol、信心等級與 prefix 開頭欄位的 csv"""
    hits = []
    for f in sorted(glob.glob("*.csv")):
        try:
            cols = pd.read_csv(f, nrows=0).columns
        except Exception:
            continue
        if "gene_symbol" in cols and LVL in cols and any(c.startswith(prefix) for c in cols):
            hits.append(f)
    return hits


def cos_rows(M, t):
    out = np.full(len(M), np.nan)
    for i, x in enumerate(M):
        m = ~np.isnan(x) & ~np.isnan(t)
        if m.sum() >= 5:
            a, b = x[m], t[m]; out[i] = a @ b / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-9)
    return out


def analyse(df, cols, label):
    M = df[cols].to_numpy(float); lv = df[LVL].to_numpy()
    hi = lv == "HIGH"; lo = lv == "VERY_LOW"
    tmpl = np.nanmean(M[hi], axis=0)
    sim = cos_rows(M, tmpl)                                  # 原做法：模板含自己
    # 留一法：HIGH 基因各自排除自己後重算模板
    n_hi = np.sum(~np.isnan(M[hi]), axis=0); s_hi = np.nansum(M[hi], axis=0)
    sim_loo = sim.copy()
    for i in np.where(hi)[0]:
        x = M[i]; ok = ~np.isnan(x)
        s = s_hi - np.where(ok, x, 0); n = n_hi - ok.astype(int)
        t = np.where(n > 0, s / np.maximum(n, 1), np.nan)
        sim_loo[i] = cos_rows(x[None, :], t)[0]
    for name, v in [("原做法（模板含自己）", sim), ("留一法（模板不含自己）", sim_loo)]:
        h = v[hi][~np.isnan(v[hi])]; l = v[lo][~np.isnan(v[lo])]
        p = mannwhitneyu(h, l, alternative="greater").pvalue
        print(f"  {label}｜{name}：HIGH n={len(h)} 平均 {h.mean():.3f}（中位數 {np.median(h):.3f}）；"
              f"VERY_LOW n={len(l)} 平均 {l.mean():.3f}（中位數 {np.median(l):.3f}）；單尾 Mann–Whitney p = {p:.2e}")


raw_f = find_table("rho_"); adj_f = find_table("radj_")
print(f"含 rho_ 欄位的檔案：{raw_f}\n含 radj_ 欄位的檔案：{adj_f}")
if not raw_f or not adj_f:
    raise SystemExit("找不到校正前或校正後的特徵檔，先停止")
RAW_PREF, ADJ_PREF = "coexpr_features_by_gene.csv", "coexpr_features_adj.csv"   # 正式分析使用的特徵檔
if RAW_PREF not in raw_f or ADJ_PREF not in adj_f:
    raise SystemExit(f"找不到正式分析的特徵檔 {RAW_PREF} 或 {ADJ_PREF}，先停止")
raw_f, adj_f = [RAW_PREF], [ADJ_PREF]
raw = pd.read_csv(raw_f[0]); adj = pd.read_csv(adj_f[0])
raw = raw[~raw["gene_symbol"].isin(ANCHORS)]
adj = adj[adj["is_anchor"] == 0] if "is_anchor" in adj.columns else adj[~adj["gene_symbol"].isin(ANCHORS)]
print(f"使用：校正前 {raw_f[0]}（{len(raw)} 個非 anchor 基因）；校正後 {adj_f[0]}（{len(adj)} 個）\n")
analyse(raw, [c for c in raw.columns if c.startswith("rho_")], "校正前")
analyse(adj, [c for c in adj.columns if c.startswith("radj_")], "校正後")
print("\n對照：Figure S3 圖上 p 值為校正前 6.8e-01、校正後 5.8e-10；2.8 寫 0.140 vs 0.008；step1b 註解寫 HIGH=0.140 vs VERY_LOW=0.008，p=5.9e-10")
