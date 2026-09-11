#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CGA MS3 — GTEx 正常組織共表達（換訊號源，衝 0.8）
================================================================
假設：TCGA 腫瘤失去節律、共表達訊號被打亂；GTEx 正常組織節律結構較完整，
      對「區分節律基因」鑑別力可能更高。這是換訊號源，不是換角度。

流程（對應 TCGA 版 Step 1 + 1b，改成 GTEx）：
  1. 讀 GCT（56200×17382），只留「候選宇宙 + 19 錨點」的列（用 Description=symbol）。
  2. 用 SampleAttributes 把樣本分到組織（SMTSD），挑跟 TCGA 對得上的正常組織。
  3. 每組織內算候選 vs 19 錨點的 Spearman（rank 一次向量化）。
  4. 跨組織 Fisher-z 加權平均（權重 n-3）。
  5. 背景校正（row-center + col z-score）→ 產出 GTEx 版 radj 特徵。

輸出：gtex_features_adj.csv（gene_symbol + gtex_radj_* + is_anchor + annotation_*）
================================================================
需要：pandas, numpy, scipy（都已有）
"""
import gzip
import numpy as np
import pandas as pd
from scipy import stats

GCT   = "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"
ATTR  = "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
FEAT  = "coexpr_features_adj.csv"     # 拿候選清單、is_anchor、annotation_
OUT   = "gtex_features_adj.csv"

ANCHORS = ["CLOCK","ARNTL","ARNTL2","PER1","PER2","PER3","CRY1","CRY2","NR1D1","NR1D2",
           "RORA","RORB","RORC","NPAS2","CSNK1D","CSNK1E","TIMELESS","BHLHE40","BHLHE41"]

# GTEx 正常組織 → 對得上 TCGA 的（用 SMTSD 精確名；挑主要幾個，樣本多）
TISSUES = [
    "Liver", "Lung", "Colon - Transverse", "Colon - Sigmoid",
    "Esophagus - Mucosa", "Stomach", "Pancreas", "Thyroid",
    "Prostate", "Breast - Mammary Tissue", "Kidney - Cortex",
    "Bladder", "Uterus", "Skin - Sun Exposed (Lower leg)",
    "Brain - Cortex", "Adrenal Gland",
]
MIN_SAMPLES = 30
MIN_EXPR_FRACTION = 0.20
LOG1P = True


def main():
    # --- 候選宇宙 + 錨點（從既有特徵表拿 symbol 清單）---
    base = pd.read_csv(FEAT)
    universe = set(base["gene_symbol"].astype(str))
    wanted = universe | set(ANCHORS)
    print(f"[universe] 候選+錨點：{len(wanted)}")

    # --- 樣本→組織對照 ---
    attr = pd.read_csv(ATTR, sep="\t", usecols=["SAMPID","SMTSD","SMAFRZE"])
    attr = attr[attr["SMAFRZE"] == "RNASEQ"]
    samp2tis = dict(zip(attr["SAMPID"], attr["SMTSD"]))
    print(f"[attr] RNASEQ 樣本數：{len(samp2tis):,}")

    # --- 讀 GCT header 取樣本欄，決定要保留哪些欄（只留目標組織的樣本）---
    with gzip.open(GCT, "rt") as f:
        f.readline(); f.readline()               # 跳過 #1.2 與 維度行
        header = f.readline().rstrip("\n").split("\t")
    sample_cols = header[2:]
    keep_samples = [s for s in sample_cols if samp2tis.get(s) in TISSUES]
    print(f"[gct] 目標組織樣本數：{len(keep_samples):,}")

    usecols = ["Description"] + keep_samples
    print("[gct] 讀取中（只留 Description + 目標樣本欄，稍等）...")
    # 分塊讀，只保留 wanted 基因，省記憶體
    chunks = []
    reader = pd.read_csv(GCT, sep="\t", skiprows=2, usecols=usecols,
                         chunksize=20000, low_memory=False)
    for ch in reader:
        ch = ch[ch["Description"].isin(wanted)]
        if len(ch):
            chunks.append(ch)
    expr = pd.concat(chunks, ignore_index=True)
    # 同 symbol 多 Ensembl → 取平均
    expr = expr.groupby("Description", as_index=True).mean(numeric_only=True)
    print(f"[gct] 命中基因：{expr.shape[0]}  樣本欄：{expr.shape[1]}")

    anchors_present = [a for a in ANCHORS if a in expr.index]
    cand_present = [g for g in expr.index if g in universe and g not in ANCHORS]
    print(f"[gct] 錨點命中：{len(anchors_present)}/19  候選命中：{len(cand_present)}")

    # --- 每組織算 Spearman（候選 vs 錨點）---
    rho_stack, n_stack, tissue_used = [], [], []
    for tis in TISSUES:
        cols = [s for s in expr.columns if samp2tis.get(s) == tis]
        if len(cols) < MIN_SAMPLES:
            continue
        M = expr[cols].T                          # 樣本 × 基因
        if LOG1P:
            M = np.log1p(M)
        # 非零比例過濾
        frac = (M > 0).mean()
        good = frac[frac >= MIN_EXPR_FRACTION].index
        a_cols = [a for a in anchors_present if a in good]
        c_cols = [g for g in cand_present if g in good]
        if not a_cols or not c_cols:
            continue
        R = M[list(dict.fromkeys(a_cols + c_cols))].rank(axis=0).to_numpy(float)
        cols_ix = {c: i for i, c in enumerate(dict.fromkeys(a_cols + c_cols))}
        n = R.shape[0]
        Rc = R - R.mean(0, keepdims=True)
        den = np.sqrt((Rc**2).sum(0)); den[den == 0] = np.nan
        ai = [cols_ix[a] for a in a_cols]; ci = [cols_ix[g] for g in c_cols]
        rho = (Rc[:, ci].T @ Rc[:, ai]) / np.outer(den[ci], den[ai])
        rho = np.clip(rho, -1, 1)
        rho_stack.append(pd.DataFrame(rho, index=c_cols, columns=a_cols))
        n_stack.append(n); tissue_used.append(tis)
        print(f"  {tis:<32} n={n:<5} 錨點{len(a_cols)} 候選{len(c_cols)}")

    if not rho_stack:
        raise SystemExit("沒有組織產出結果，檢查組織名或樣本數門檻。")
    print(f"[聚合] 用了 {len(tissue_used)} 個組織：{tissue_used}")

    # --- 跨組織 Fisher-z 加權平均 ---
    genes = sorted(set().union(*[d.index for d in rho_stack]))
    anch = sorted(set().union(*[d.columns for d in rho_stack]))
    Z = np.full((len(rho_stack), len(genes), len(anch)), np.nan)
    W = np.zeros_like(Z)
    for k, (d, n) in enumerate(zip(rho_stack, n_stack)):
        r = np.clip(d.reindex(index=genes, columns=anch).to_numpy(float), -0.999, 0.999)
        z = np.arctanh(r); Z[k] = z
        W[k] = np.where(np.isfinite(z), max(n-3, 0), 0.0)
    zbar = np.nansum(W*Z, 0) / np.where(W.sum(0) == 0, np.nan, W.sum(0))
    pooled = pd.DataFrame(np.tanh(zbar), index=genes, columns=anch)

    # --- 背景校正 ---
    R = pooled.to_numpy(float)
    Rc = R - np.nanmean(R, 1, keepdims=True)
    Rz = (Rc - np.nanmean(Rc, 0)) / (np.nanstd(Rc, 0) + 1e-9)

    out = pd.DataFrame({"gene_symbol": genes})
    for i, a in enumerate(anch):
        out[f"gtex_radj_{a}"] = Rz[:, i]
    # 併回 is_anchor / annotation_
    meta = base[["gene_symbol","is_anchor"] + [c for c in base.columns if c.startswith("annotation_")]]
    out = out.merge(meta, on="gene_symbol", how="left")
    out["is_anchor"] = out["is_anchor"].fillna(0).astype(int)
    out.to_csv(OUT, index=False, encoding="utf-8-sig")
    print(f"\n[out] {OUT}：{out.shape[0]} 基因 × {out.shape[1]} 欄")
    print("完成。下一步把 gtex_radj_* 餵分類器（正例=CircaDB）看 AUC。")


if __name__ == "__main__":
    main()