#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CGA MS3 — Step 1b：共表達特徵的背景校正
================================================================
問題：原始 signed rho 有很強的「背景相關地板」（腫瘤純度、增殖等共同因素
      讓任兩基因都有 ~0.15 的相關），把時鐘基因的相位模式淹掉。
      sanity check 證實：未校正時 HIGH vs VERY_LOW 幾乎無鑑別力（p=0.72）；
      校正後相似度 HIGH=0.140 vs VERY_LOW=0.008，p=5.9e-10 → 訊號回來。

校正兩步：
  1. row-center：每個基因減掉「自己對全部 19 錨點的平均 rho」
     → 移除「這基因跟大盤一起漲跌」的成分，留下「特別偏向哪些錨點」。
  2. column z-score：每個錨點欄再標準化
     → 消掉不同錨點本身的高低基準（如 CLOCK 天生相關就高）。

輸出 radj_<ANCHOR> 共 19 欄，作為 Step 2 分類器的特徵。
annotation_* 原樣保留（只給 Step 3 挑候選，永不進特徵）。

輸入：coexpr_features_by_gene.csv（Step 1 產出）
輸出：coexpr_features_adj.csv
================================================================
"""
import sys
import numpy as np
import pandas as pd

IN_CSV  = "coexpr_features_by_gene.csv"
OUT_CSV = "coexpr_features_adj.csv"

# 19 個 canonical anchor（官方符號；Step 2 會用它排除正例、並當 sanity check）
ANCHORS = {
    "CLOCK", "ARNTL", "ARNTL2", "PER1", "PER2", "PER3", "CRY1", "CRY2",
    "NR1D1", "NR1D2", "RORA", "RORB", "RORC", "NPAS2", "CSNK1D", "CSNK1E",
    "TIMELESS", "BHLHE40", "BHLHE41",
}


def main():
    df = pd.read_csv(IN_CSV)
    rho_cols = [c for c in df.columns if c.startswith("rho_")]
    if not rho_cols:
        sys.exit("找不到 rho_* 欄，請確認輸入是 Step 1 的 coexpr_features_by_gene.csv")
    print(f"[in] {len(df)} 基因，{len(rho_cols)} 個 signed rho 錨點欄")

    R = df[rho_cols].to_numpy(dtype=float)

    # --- 校正 1：row-center（每基因減自身對全錨點的平均）---
    row_mean = np.nanmean(R, axis=1, keepdims=True)
    Rc = R - row_mean

    # --- 校正 2：column z-score（每錨點欄標準化）---
    col_mean = np.nanmean(Rc, axis=0)
    col_std  = np.nanstd(Rc, axis=0)
    col_std[col_std == 0] = np.nan
    Rz = (Rc - col_mean) / col_std

    # 組輸出：radj_* 校正特徵 + gene_symbol + annotation_*
    out = pd.DataFrame({"gene_symbol": df["gene_symbol"].values})
    for i, c in enumerate(rho_cols):
        out[c.replace("rho_", "radj_")] = Rz[:, i]

    # 標記是否為 anchor（Step 2 用來排除正例 / 挑出來做 sanity check）
    out["is_anchor"] = out["gene_symbol"].isin(ANCHORS).astype(int)

    # annotation_* 原樣帶過（只給 Step 3，訓練前 drop）
    for c in df.columns:
        if c.startswith("annotation_"):
            out[c] = df[c].values

    out.to_csv(OUT_CSV, index=False, encoding="utf-8-sig")
    print(f"[out] {OUT_CSV}：{out.shape[0]} 基因 × {out.shape[1]} 欄")
    print(f"      校正特徵欄：{[c for c in out.columns if c.startswith('radj_')]}")
    print(f"      anchor 標記數：{int(out['is_anchor'].sum())}（應為 19）")
    print("完成。radj_* 進 Step 2 當特徵；annotation_* 只給 Step 3；is_anchor 用來排除正例。")


if __name__ == "__main__":
    main()