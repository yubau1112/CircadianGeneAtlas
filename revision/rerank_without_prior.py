#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CGA 第 9 項：拿掉「已知時鐘基因加分」與「基因名稱加分」後重新排序。
規則照 circadian_gene_analyzer_dual_optimized_20250721.py：
  文獻情境分數：已知時鐘基因 +20
  NCBI 摘要分數：已知時鐘基因 +25；基因名稱符合樣式再加分（只在有 NCBI 摘要時計算）
執行位置：/home/yubau/aging/20250519_circadian
輸出：cga_rerank_no_known_bonus.tsv
"""
import re
import pandas as pd

SRC = "circadian_genes_results_dual_optimized.csv"
OUT = "cga_rerank_no_known_bonus.tsv"

# ── 與分析腳本相同的清單（第 67–75 行）──
VERIFIED = ["CLOCK", "ARNTL", "BMAL1", "ARNTL2", "BMAL2",
            "PER1", "PER2", "PER3", "CRY1", "CRY2",
            "NR1D1", "REV-ERBA", "NR1D2", "REV-ERBB",
            "RORA", "RORB", "RORC", "NPAS2", "CSNK1D", "CSNK1E",
            "TIMELESS", "DEC1", "DEC2", "BHLHE40", "BHLHE41"]
# ── 與分析腳本相同的基因名稱加分（第 285–289 行），第一個符合者生效 ──
NAME_PATTERNS = [(r"\bCLOCK\b", 15), (r"\bBMAL1\b", 15), (r"\bARNTL\b", 15),
                 (r"\bPER[1-3]\b", 12), (r"\bCRY[1-2]\b", 12),
                 (r"\bREV-ERB\b", 10), (r"\bROR\b", 10),
                 (r"\bNPAS2\b", 8), (r"\bCSNK1[DE]\b", 8),
                 (r"\bTIMELESS\b", 8), (r"\bDEC[1-2]\b", 6)]
ANCHORS = ["CLOCK", "ARNTL", "ARNTL2", "PER1", "PER2", "PER3", "CRY1", "CRY2",
           "NR1D1", "NR1D2", "RORA", "RORB", "RORC", "NPAS2", "CSNK1D", "CSNK1E",
           "TIMELESS", "BHLHE40", "BHLHE41"]
NO_SUMMARY = {"", "無摘要", "無法獲取", "獲取失敗 - 達到最大重試次數"}


def name_bonus(symbol):
    for pat, pts in NAME_PATTERNS:
        if re.search(pat, symbol, re.IGNORECASE):
            return pts
    return 0


def level(score):
    if score >= 40:
        return "HIGH"
    if score >= 20:
        return "MEDIUM"
    if score >= 8:
        return "LOW"
    return "VERY_LOW"


df = pd.read_csv(SRC, encoding="utf-8-sig", escapechar="\\", dtype=str, keep_default_na=False)
need = ["Gene_ID", "Gene_Symbol", "Literature_Score", "NCBI_Score", "Combined_Score", "NCBI_Summary"]
missing = [c for c in need if c not in df.columns]
if missing:
    raise SystemExit(f"缺欄位：{missing}")
if len(df) != 12629:
    raise SystemExit(f"列數不是 12,629，而是 {len(df):,}；CSV 解析可能有問題，先停止")
for c in ["Literature_Score", "NCBI_Score", "Combined_Score"]:
    df[c] = pd.to_numeric(df[c], errors="raise")
gap = (df["Literature_Score"] + df["NCBI_Score"] - df["Combined_Score"]).abs().max()
if gap > 1e-6:
    raise SystemExit(f"文獻 + NCBI 不等於綜合分數（最大差 {gap}），先停止")

df["is_known"] = df["Gene_Symbol"].str.upper().isin(VERIFIED)
df["has_summary"] = ~df["NCBI_Summary"].str.strip().isin(NO_SUMMARY)
df["name_bonus"] = [name_bonus(s) if h else 0 for s, h in zip(df["Gene_Symbol"], df["has_summary"])]

df["lit_adj"] = (df["Literature_Score"] - 20 * df["is_known"]).clip(lower=0)
df["ncbi_adj"] = (df["NCBI_Score"] - 25 * (df["is_known"] & df["has_summary"]) - df["name_bonus"]).clip(lower=0)
df["combined_adj"] = df["lit_adj"] + df["ncbi_adj"]
df["capped"] = (df["Literature_Score"] >= 100) | (df["NCBI_Score"] >= 100)

df["rank_orig"] = df["Combined_Score"].rank(ascending=False, method="min").astype(int)
df["rank_adj"] = df["combined_adj"].rank(ascending=False, method="min").astype(int)
df["level_orig"] = df["Combined_Score"].map(level)
df["level_adj"] = df["combined_adj"].map(level)

# ── 核對：Table 2 的 CLOCK 為文獻 88、NCBI 76 ──
clock = df[df["Gene_Symbol"] == "CLOCK"]
if len(clock) != 1:
    raise SystemExit("找不到唯一的 CLOCK 列")
c = clock.iloc[0]
print(f"核對 CLOCK：原始 {c.Literature_Score} + {c.NCBI_Score} = {c.Combined_Score}；"
      f"扣除後 {c.lit_adj} + {c.ncbi_adj} = {c.combined_adj}（預期 68 + 36 = 104）")

n = len(df)
print(f"\n總基因數 {n:,}；有分數上限 100 被截的基因 {int(df['capped'].sum())} 個（這些扣分後的數值為近似）")

anc = df[df["Gene_Symbol"].isin(ANCHORS)].sort_values("rank_adj")
print(f"\n=== 19 個時鐘基因：原始排名 → 扣除加分後排名（找到 {len(anc)} 個）===")
print(f"{'基因':<10}{'原始分數':>8}{'原始名次':>8}{'扣後文獻':>8}{'扣後NCBI':>9}{'扣後分數':>8}{'扣後名次':>8}  扣後等級")
for _, r in anc.iterrows():
    print(f"{r.Gene_Symbol:<10}{r.Combined_Score:>8.1f}{r.rank_orig:>8}{r.lit_adj:>8.1f}{r.ncbi_adj:>9.1f}"
          f"{r.combined_adj:>8.1f}{r.rank_adj:>8}  {r.level_adj}")

for k in (15, 19, 50, 100, 377):
    print(f"扣除後前 {k:>3} 名中的時鐘基因數：{int((anc['rank_adj'] <= k).sum())} / {len(anc)}")
print(f"扣除後仍為 HIGH 的時鐘基因：{int((anc['level_adj'] == 'HIGH').sum())} / {len(anc)}")
print(f"HIGH 基因總數：原始 {int((df['level_orig'] == 'HIGH').sum())} → 扣除後 {int((df['level_adj'] == 'HIGH').sum())}")

print("\n=== 扣除加分後的前 20 名 ===")
for _, r in df.sort_values(["combined_adj", "Combined_Score"], ascending=False).head(20).iterrows():
    tag = "時鐘基因" if r.Gene_Symbol in ANCHORS else ""
    print(f"{r.rank_adj:>4}. {r.Gene_Symbol:<12} {r.combined_adj:>6.1f}（原始 {r.Combined_Score:.1f}，原名次 {r.rank_orig}） {tag}")

odd = df[((df["is_known"]) | (df["name_bonus"] > 0)) & (~df["Gene_Symbol"].isin(ANCHORS))]
print(f"\n=== 不是 19 個時鐘基因、卻拿到已知基因或名稱加分的基因：{len(odd)} 個 ===")
for _, r in odd.iterrows():
    print(f"  {r.Gene_ID:<10} {r.Gene_Symbol:<12} 已知清單={bool(r.is_known)} 名稱加分={r.name_bonus} 原始分數={r.Combined_Score}")

cols = ["Gene_ID", "Gene_Symbol", "Literature_Score", "NCBI_Score", "Combined_Score", "rank_orig", "level_orig",
        "is_known", "has_summary", "name_bonus", "lit_adj", "ncbi_adj", "combined_adj", "rank_adj", "level_adj", "capped"]
df.sort_values("rank_adj")[cols].to_csv(OUT, sep="\t", index=False, encoding="utf-8")
print(f"\n已輸出 {OUT}（{n:,} 列）")