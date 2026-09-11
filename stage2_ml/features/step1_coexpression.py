#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CGA MS3 — Step 1：TCGA 共表達分析（co-expression）
================================================================
每個候選基因，對 canonical core-clock「錨點基因」逐一算 Spearman 相關，
在每種癌症內各算一次，再跨 24 種 TCGA 癌症聚合，
輸出「每基因 × 每錨點」的相關特徵表，供 Step 2 監督式分類器使用。

【方法紀律 — 很重要】
  分類器特徵「只用共表達／組學」。文獻分數（Combined_Score）**絕不**當模型特徵，
  以避免循環論證。本程式只把文獻分數以 annotation_* 欄位「附在輸出旁邊」，
  純粹給 Step 3 挑候選用（高共表達訊號 + 低文獻分數 = 被低估的候選基因）。
  → 訓練時務必把所有 annotation_* 欄位排除在 X 之外。

【資料來源】
  - 表達：tcga 庫（get_gene_data.php 的四表 join，單位 tpm_unstranded / TPM）
  - 基因宇宙 + 文獻分數：circadian 庫 circadian_genes_results

【錨點基因（anchor）】
  用 VERIFIED_CIRCADIAN_GENES（canonical TTFL），去別名後 19 個官方符號。
  不使用舊投影片的「32」（那是手稿功能分類的人工數字，無清單支撐）。

作者交接：Yu-Pao Chou / CircadianGeneAtlas
需要套件：pyodbc, pandas, numpy, scipy
  pip install pyodbc pandas numpy scipy
================================================================
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy import stats

# ============================================================
#                         CONFIG
#      （要調的東西都在這一區，改完就能跑）
# ============================================================

# ---- SQL Server 連線（來自 db_config.php，直接用）----
SQL_SERVER   = os.environ.get("CGA_SQL_SERVER", "localhost,1433")    # host,port
SQL_USER=os.environ.get("CGA_SQL_USER", "")
SQL_PASSWORD=os.environ.get("CGA_SQL_PASSWORD", "")
DB_TCGA      = "TCGA"                   # 表達（注意大小寫：TCGA）
DB_CIRCADIAN = "Circadian"             # 分數 / 基因宇宙（注意大小寫：Circadian）
ODBC_DRIVER  = "ODBC Driver 18 for SQL Server"   # 沒裝 18 就改 17

# ---- 錨點基因：直接取自 circadian_gene_analyzer_dual_optimized_20250721.py line 67 ----
# VERIFIED_CIRCADIAN_GENES 原文（25 個符號，含別名），逐字照抄：
VERIFIED_CIRCADIAN_GENES = [
    "CLOCK", "ARNTL", "BMAL1", "ARNTL2", "BMAL2",
    "PER1", "PER2", "PER3",
    "CRY1", "CRY2",
    "NR1D1", "REV-ERBA", "NR1D2", "REV-ERBB",
    "RORA", "RORB", "RORC",
    "NPAS2", "CSNK1D", "CSNK1E",
    "TIMELESS", "DEC1", "DEC2", "BHLHE40", "BHLHE41",
]
# 別名 → 官方符號（GENCODEv36 的 gene_name 用官方符號，比對前先去別名）
ALIAS_TO_OFFICIAL = {
    "BMAL1": "ARNTL", "BMAL2": "ARNTL2",
    "REV-ERBA": "NR1D1", "REV-ERBB": "NR1D2",
    "DEC1": "BHLHE40", "DEC2": "BHLHE41",
}
# 去別名 + 保序去重 → 正好 19 個官方符號
ANCHOR_GENES = list(dict.fromkeys(
    ALIAS_TO_OFFICIAL.get(g, g) for g in VERIFIED_CIRCADIAN_GENES
))

# ============================================================
#  煙霧測試開關：True = 小規模先測（2 癌別 + 500 基因），
#                False = 全量（24 癌 + 全部基因）。
#  測通了、確認 console 正常，把這行改 False 再跑一次即可。
# ============================================================
SMOKE_TEST = False

# ---- 基因宇宙（要打分的候選）----
# 預設：circadian_genes_results 裡所有有 Combined_Score 的基因（含 VERY_LOW），
# 因為 Step 3 要找的「低文獻分數卻像時鐘基因」正落在這個集合裡。
_UNIVERSE_BASE = f"SELECT Gene_Symbol, Combined_Score, Confidence_Level FROM [{DB_CIRCADIAN}].[dbo].[circadian_genes_results]"
if SMOKE_TEST:
    # 只取前 500 個基因（錨點會另外自動併入，不受此限制）
    UNIVERSE_SQL = _UNIVERSE_BASE + " ORDER BY Gene_Symbol OFFSET 0 ROWS FETCH NEXT 500 ROWS ONLY"
else:
    UNIVERSE_SQL = _UNIVERSE_BASE

# ---- 樣本基礎 ----
# 預設只用 Primary Tumor：樣本多、多數癌別 n 足夠，Step 1 特徵問的是「腫瘤內是否與時鐘共表達」。
# Solid Tissue Normal 保留給之後的網路 rewiring 分支；要一起算就把它加進 list。
SAMPLE_TYPES = ["Primary Tumor"]        # 例：["Primary Tumor", "Solid Tissue Normal"]

# ---- 門檻 ----
MIN_SAMPLES_PER_CANCER = 20   # 某癌某錨點樣本數 < 此值就跳過（小 n 的 Spearman 是雜訊）
MIN_EXPR_FRACTION      = 0.20 # 某癌某基因非零 TPM 比例 < 此值就跳過（太多 0 → 排名不可靠）
SIG_ALPHA              = 0.05 # 「顯著癌症比例」特徵用的雙尾門檻
STRONG_RHO             = 0.30 # 摘要特徵「強相關錨點數」的 |rho| 門檻

# ---- TCGA 癌症（來自 get_gene_data.php 寫死的 VALUES；共 24 種）----
if SMOKE_TEST:
    CANCERS = ["TCGA-LIHC", "TCGA-BRCA"]     # 煙霧測試：先兩個
else:
    CANCERS = [
        "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL", "TCGA-COAD", "TCGA-ESCA",
        "TCGA-GBM",  "TCGA-HNSC", "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LIHC",
        "TCGA-LUAD", "TCGA-LUSC", "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ",
        "TCGA-SARC", "TCGA-SKCM", "TCGA-STAD", "TCGA-THCA", "TCGA-THYM", "TCGA-UCEC",
    ]

# ---- 聚合方式 ----
# 跨癌聚合：每癌 Spearman rho → Fisher z（arctanh），以 (n-3) 加權平均，再 tanh 回來。
# rho 帶正負號保留（時鐘基因有反相位關係，例：PER 與 BMAL 反相；signed rho 才抓得到）。
LOG1P_TPM = True   # 先 log1p(TPM) 再排名。對 Spearman 其實不改變 rho（單調變換），
                   # 但保留成開關，之後若改用 Pearson 會用到。

# ---- 輸出 ----
OUT_LONG     = "coexpr_per_cancer_long.csv"     # 稽核用：每 (基因, 錨點, 癌症) 一列
OUT_FEATURES = "coexpr_features_by_gene.csv"    # Step 2 用：每基因一列的特徵表

# ============================================================
#                       實作
# ============================================================

def get_conn(database):
    import pyodbc
    conn_str = (
        f"DRIVER={{{ODBC_DRIVER}}};SERVER={SQL_SERVER};DATABASE={database};"
        f"UID={SQL_USER};PWD={SQL_PASSWORD};TrustServerCertificate=yes;Encrypt=no"
    )
    return pyodbc.connect(conn_str)


def _query_df(database, sql, params=None):
    """用 pyodbc cursor 直接查，回 DataFrame（避開 pd.read_sql 的 SQLAlchemy 警告）。"""
    with get_conn(database) as c:
        cur = c.cursor()
        cur.execute(sql, params or [])
        cols = [d[0] for d in cur.description]
        rows = cur.fetchall()
    return pd.DataFrame.from_records([tuple(r) for r in rows], columns=cols)


def load_gene_universe():
    """從 circadian 庫拉候選基因 + 文獻分數（分數僅供 Step 3，不當特徵）。"""
    df = _query_df(DB_CIRCADIAN, UNIVERSE_SQL)
    df["Gene_Symbol"] = df["Gene_Symbol"].astype(str).str.strip()
    df = df.drop_duplicates(subset="Gene_Symbol")
    print(f"[universe] circadian_genes_results 基因數：{len(df)}")
    return df


# SQL Server 單句參數上限 2100，基因清單分批查（每批 1000）再合併
GENE_BATCH = 1000


def load_expression_for_cancer(cancer, wanted_genes):
    """
    拉單一癌症、指定樣本型別、限定在 wanted_genes 的逐樣本 TPM（長格式）。
    回傳 DataFrame：gene_name / sample_uuid / tpm
    join 圖同 get_gene_data.php：
        A.File_UUID = B.num；A.gene_id = C.num；B.File_UUID = D.File_UUID
    基因清單分批（每批 GENE_BATCH），避開 2100 參數上限。
    """
    placeholders_st = ",".join("?" for _ in SAMPLE_TYPES)
    frames = []
    for i in range(0, len(wanted_genes), GENE_BATCH):
        batch = wanted_genes[i:i + GENE_BATCH]
        placeholders_gn = ",".join("?" for _ in batch)
        sql = f"""
            SELECT C.[gene_name]              AS gene_name,
                   B.[File_UUID]              AS sample_uuid,
                   A.[tpm_unstranded]         AS tpm
            FROM [{DB_TCGA}].[dbo].[Gene_Expression] A
            INNER JOIN [{DB_TCGA}].[dbo].[File_UUID_Gene_Expression] B ON A.[File_UUID] = B.[num]
            INNER JOIN [{DB_TCGA}].[dbo].[Gene_Info_GENCODEv36]      C ON A.[gene_id]   = C.[num]
            INNER JOIN [{DB_TCGA}].[dbo].[Sample_Type_Gene_Expression] D ON B.[File_UUID] = D.[File_UUID]
            WHERE D.[Project_ID] = ?
              AND D.[Sample_Type] IN ({placeholders_st})
              AND A.[tpm_unstranded] IS NOT NULL
              AND A.[tpm_unstranded] >= 0
              AND C.[gene_name] IN ({placeholders_gn})
        """
        params = [cancer, *SAMPLE_TYPES, *batch]
        frames.append(_query_df(DB_TCGA, sql, params))
    if frames:
        df = pd.concat(frames, ignore_index=True)
    else:
        df = pd.DataFrame(columns=["gene_name", "sample_uuid", "tpm"])
    df["tpm"] = pd.to_numeric(df["tpm"], errors="coerce")
    return df.dropna(subset=["tpm"])


def build_matrix(long_df):
    """長格式 → sample(列) × gene(欄) 的 TPM 矩陣。重複值取平均（保險）。"""
    m = long_df.pivot_table(index="sample_uuid", columns="gene_name",
                            values="tpm", aggfunc="mean")
    if LOG1P_TPM:
        m = np.log1p(m)
    return m


def spearman_block(matrix, anchor_cols, cand_cols):
    """
    在單一癌症的 sample×gene 矩陣上，算每個候選基因 vs 每個錨點的 Spearman。
    做法：整體 rank 一次（= Spearman），再向量化算 Pearson 與 t 近似 p 值。
    回傳：rho(cand×anchor)、p(cand×anchor)、每欄有效樣本數 n_col。
    """
    # 每基因非零比例過濾（太多 0 排名不可靠）
    frac_nonzero = (np.expm1(matrix) > 0).mean() if LOG1P_TPM else (matrix > 0).mean()
    keep = frac_nonzero[frac_nonzero >= MIN_EXPR_FRACTION].index
    anchor_cols = [a for a in anchor_cols if a in keep]
    cand_cols   = [g for g in cand_cols   if g in keep]
    if not anchor_cols or not cand_cols:
        return None

    sub = matrix[list(dict.fromkeys(anchor_cols + cand_cols))]
    ranks = sub.rank(axis=0, method="average")     # Spearman = Pearson(ranks)
    R = ranks.to_numpy(dtype=float)
    cols = list(ranks.columns)
    idx = {c: i for i, c in enumerate(cols)}
    n = R.shape[0]
    if n < MIN_SAMPLES_PER_CANCER:
        return None

    Rc = R - R.mean(axis=0, keepdims=True)
    denom = np.sqrt((Rc ** 2).sum(axis=0))
    denom[denom == 0] = np.nan                     # 零變異基因 → NaN

    a_i = [idx[a] for a in anchor_cols]
    c_i = [idx[g] for g in cand_cols]
    Ra, Rc_cand = Rc[:, a_i], Rc[:, c_i]
    da, dc = denom[a_i], denom[c_i]
    rho = (Rc_cand.T @ Ra) / np.outer(dc, da)      # (cand × anchor)
    rho = np.clip(rho, -1.0, 1.0)

    with np.errstate(invalid="ignore", divide="ignore"):
        t = rho * np.sqrt((n - 2) / (1.0 - rho ** 2))
        p = 2 * stats.t.sf(np.abs(t), df=n - 2)

    rho_df = pd.DataFrame(rho, index=cand_cols, columns=anchor_cols)
    p_df   = pd.DataFrame(p,   index=cand_cols, columns=anchor_cols)
    return rho_df, p_df, n


def fisher_aggregate(rho_stack, n_stack):
    """
    跨癌聚合：Fisher z 加權平均（權重 n-3），再 tanh 回 rho。
    rho_stack / n_stack：list of (cand×anchor) DataFrame，index/columns 對齊。
    回傳：pooled_rho、mean_frac_sig(交由外層算)、n_cancers_used。
    """
    # 對齊成 3D
    genes = sorted(set().union(*[df.index for df in rho_stack]))
    anchors = sorted(set().union(*[df.columns for df in rho_stack]))
    Z = np.full((len(rho_stack), len(genes), len(anchors)), np.nan)
    W = np.zeros_like(Z)
    for k, (rdf, ndf) in enumerate(zip(rho_stack, n_stack)):
        rr = rdf.reindex(index=genes, columns=anchors)
        r = np.clip(rr.to_numpy(dtype=float), -0.999, 0.999)
        z = np.arctanh(r)
        Z[k] = z
        # ndf 是純量 n（每癌一個），廣播成權重
        w = max(ndf - 3, 0)
        W[k] = np.where(np.isfinite(z), w, 0.0)
    Wsum = W.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        zbar = np.nansum(W * Z, axis=0) / np.where(Wsum == 0, np.nan, Wsum)
    pooled = np.tanh(zbar)
    n_used = np.isfinite(Z).sum(axis=0)
    return (pd.DataFrame(pooled, index=genes, columns=anchors),
            pd.DataFrame(n_used, index=genes, columns=anchors))


def main():
    # 1) 基因宇宙 + 錨點
    uni = load_gene_universe()
    universe_syms = set(uni["Gene_Symbol"])
    anchors = ANCHOR_GENES                      # 已去別名，19 個官方符號
    wanted = sorted(universe_syms | set(anchors))
    print(f"[anchors] {len(anchors)} 個：{anchors}")
    print(f"[wanted]  拉表達的基因總數（宇宙 ∪ 錨點）：{len(wanted)}")

    # 2) 逐癌拉表達、算 Spearman
    rho_stack, p_frac_stack, n_stack, processed_cancers = [], [], [], []
    anchors_seen = set()
    for cancer in CANCERS:
        long_df = load_expression_for_cancer(cancer, wanted)
        if long_df.empty:
            print(f"  {cancer}: 無資料，跳過")
            continue
        M = build_matrix(long_df)
        present_anchors = [a for a in anchors if a in M.columns]
        cand_cols = [g for g in M.columns if g in universe_syms]
        res = spearman_block(M, present_anchors, cand_cols)
        if res is None:
            print(f"  {cancer}: 樣本或基因不足，跳過")
            continue
        rho_df, p_df, n = res
        rho_stack.append(rho_df)
        p_frac_stack.append((p_df < SIG_ALPHA).astype(float))  # 顯著 = 1
        n_stack.append(n)
        processed_cancers.append(cancer)
        anchors_seen |= set(rho_df.columns)
        print(f"  {cancer}: n={n}, 錨點 {len(rho_df.columns)}, 候選 {rho_df.shape[0]}")

    if not rho_stack:
        sys.exit("沒有任何癌症產出結果，請檢查連線 / schema / 樣本型別。")

    # 錨點解析報告（別名或 GENCODE 命名不合會在這裡現形）
    print(f"[錨點解析] 找到 {len(anchors_seen)}/{len(anchors)} 個：{sorted(anchors_seen)}")
    missing = sorted(set(anchors) - anchors_seen)
    if missing:
        print(f"[警告] 這些錨點在表達庫沒找到（可能是別名/命名不同）：{missing}")

    # 3) 跨癌聚合
    pooled_rho, n_used = fisher_aggregate(rho_stack, n_stack)
    # 顯著癌症比例（不分正負；時鐘基因有反相位，用雙尾）
    genes = pooled_rho.index
    anch = pooled_rho.columns
    sig_sum = np.zeros((len(genes), len(anch)))
    val_cnt = np.zeros_like(sig_sum)
    for sig_df, rho_df in zip(p_frac_stack, rho_stack):
        s = sig_df.reindex(index=genes, columns=anch)
        v = rho_df.reindex(index=genes, columns=anch).notna()
        sig_sum += np.nan_to_num(s.to_numpy(dtype=float))
        val_cnt += v.to_numpy(dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        frac_sig = pd.DataFrame(sig_sum / np.where(val_cnt == 0, np.nan, val_cnt),
                                index=genes, columns=anch)

    # 4) 組特徵表
    feat = pd.DataFrame(index=genes)
    for a in anch:                                   # 每錨點的 signed pooled rho（主特徵）
        feat[f"rho_{a}"] = pooled_rho[a]
    absr = pooled_rho.abs()
    feat["mean_abs_rho"]   = absr.mean(axis=1)
    feat["max_abs_rho"]    = absr.max(axis=1)
    feat["n_strong_anchor"] = (absr >= STRONG_RHO).sum(axis=1)
    feat["mean_frac_sig"]  = frac_sig.mean(axis=1)
    feat["mean_n_cancers"] = n_used.mean(axis=1)

    # annotation_*（給 Step 3，絕不當特徵；訓練前要 drop）
    ann = uni.set_index("Gene_Symbol").reindex(genes)
    feat["annotation_combined_score"]   = ann["Combined_Score"].values
    feat["annotation_confidence_level"] = ann["Confidence_Level"].values

    feat = feat.reset_index().rename(columns={"index": "gene_symbol"})
    feat.to_csv(OUT_FEATURES, index=False, encoding="utf-8-sig")
    print(f"[輸出] 特徵表：{OUT_FEATURES}（{feat.shape[0]} 基因 × {feat.shape[1]} 欄）")

    # 5) 長格式稽核檔（每 基因×錨點×癌症 一列，含真實癌別名與樣本數）
    long_rows = []
    for cancer, rho_df, n in zip(processed_cancers, rho_stack, n_stack):
        tmp = rho_df.reset_index().melt(id_vars="index", var_name="anchor", value_name="rho")
        tmp = tmp.rename(columns={"index": "gene_symbol"})
        tmp["cancer"] = cancer
        tmp["n_samples"] = n
        long_rows.append(tmp)
    long_out = pd.concat(long_rows, ignore_index=True)
    long_out = long_out[["gene_symbol", "anchor", "cancer", "n_samples", "rho"]]
    long_out.to_csv(OUT_LONG, index=False, encoding="utf-8-sig")
    print(f"[輸出] 稽核長檔：{OUT_LONG}（{len(long_out)} 列）")

    print("\n完成。features 檔可直接進 Step 2；訓練時記得 drop 掉 annotation_* 欄位。")


if __name__ == "__main__":
    main()