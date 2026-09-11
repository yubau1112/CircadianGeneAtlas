#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CGA 第 8 項：以獨立的人類血液時序轉錄體驗證候選（GSE39445，Möller-Levet et al. 2013）v2
v4：標題可無「Subject_」前綴（11 個樣本如 AF0033_SE_T2）；v3：受試者編號改由樣本標題取得（部分樣本特徵缺 subject）；v2：樣本特徵逐格讀取欄位名稱（v1 有 11 個樣本因欄位順序不同被漏掉）；三種節律門檻；加 Mann–Whitney 排名比較；加入選頻率分層
設計：
  時間 = 樣本的 circadianphase（相對於 DLMO 的生理時鐘相位，小時）
  每個條件（Sleep Extension 主分析、Sleep Restriction 重複驗證）各自分析
  每個基因：先扣掉每位受試者自己的平均（固定效應），模型含 cos、sin（24 小時週期）與清醒時數（線性趨勢）
  以 F 檢定比較「含 cos、sin」與「只含清醒時數」，得到節律 p 值，再做 BH 校正
  探針 → 基因：排除控制探針；同一基因取平均表現最高的探針
比較「有節律」的比例：候選 399、穩定核心（40 次中 ≥ 36 次）、公平對照池（文獻 LOW／VERY_LOW、非 anchor、非 CircaDB、非候選）、
  正對照（19 個 anchor、CircaDB 訓練正例、全部 CircaDB）
IEG 清單為暫定的常見立即早期基因（與 cga_ieg_precision.py 相同）。

執行位置：/home/yubau/aging/20250519_circadian/Phase2/timeseries_validation
需要：GSE39445_series_matrix.txt.gz、GPL15331_family.soft.gz
      ../circadian_pipeline/step_final_scores.csv、step_final_candidates.csv、cand_stability_genes.tsv
輸出：ts_GSE39445_gene_stats_<條件>.tsv、ts_GSE39445_summary.tsv
"""
import gzip, csv, io
import numpy as np, pandas as pd
from scipy import stats

SERIES = "GSE39445_series_matrix.txt.gz"
PLATFORM = "GPL15331_family.soft.gz"
PIPE = "../circadian_pipeline/"
Q_CUT, AMP_CUT, PERIOD = 0.05, 0.1, 24.0
IEG = ["FOS", "FOSB", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND", "EGR1", "EGR2", "EGR3", "EGR4",
       "IER2", "IER3", "IER5", "NR4A1", "NR4A2", "NR4A3", "ATF3", "DUSP1", "DUSP2", "DUSP5", "DUSP6",
       "KLF2", "KLF4", "KLF6", "KLF10", "ZFP36", "ZFP36L1", "ZFP36L2", "CCN1", "CCN2", "BTG2",
       "PPP1R15A", "GADD45B", "GADD45G", "MYC", "ARC", "NPAS4", "SGK1", "CEBPB", "CEBPD", "SOCS3",
       "HES1", "PLK2", "CSRNP1", "NFKBIA", "IRF1", "TNFAIP3", "KDM6B", "SIK1"]


def read_series(path):
    meta, table, chars = {}, [], None
    in_tab = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                in_tab = True; continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_tab:
                table.append(line)
            elif line.startswith("!Sample_"):
                row = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
                key, vals = row[0], row[1:]
                if key == "!Sample_characteristics_ch1":
                    if chars is None:
                        chars = [dict() for _ in vals]
                    for i, v in enumerate(vals):
                        if ":" in v:
                            k, val = v.split(":", 1)
                            chars[i][k.strip()] = val.strip()
                elif key in ("!Sample_geo_accession", "!Sample_title"):
                    meta[key[8:]] = vals
    expr = pd.read_csv(io.StringIO("".join(table)), sep="\t", index_col=0)
    expr.index = expr.index.astype(str)
    df = pd.DataFrame(meta)
    if chars is not None:
        df = pd.concat([df, pd.DataFrame(chars)], axis=1)
    return df, expr


def read_platform(path):
    rows, in_tab = [], False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!platform_table_begin"):
                in_tab = True; continue
            if line.startswith("!platform_table_end"):
                break
            if in_tab:
                rows.append(line)
    ann = pd.read_csv(io.StringIO("".join(rows)), sep="\t", dtype=str, usecols=["ID", "CONTROL_TYPE", "GENE_SYMBOL"])
    return ann


def within_subject(M, subj):
    """每個受試者各自扣掉平均（M：樣本 × 變數）"""
    M = M.copy()
    for s in np.unique(subj):
        idx = subj == s
        M[idx] -= M[idx].mean(axis=0)
    return M


def rhythm_test(Y, t, awake, subj):
    """Y：樣本 × 基因（log2）。回傳每個基因的 p、振幅、峰值相位"""
    w = 2 * np.pi * t / PERIOD
    Xf = within_subject(np.column_stack([np.cos(w), np.sin(w), awake]), subj)
    Xr = Xf[:, [2]]
    Yc = within_subject(Y, subj)
    n, k_sub = len(t), len(np.unique(subj))

    def rss(X):
        B, *_ = np.linalg.lstsq(X, Yc, rcond=None)
        R = Yc - X @ B
        return (R ** 2).sum(axis=0), B
    rss_f, Bf = rss(Xf)
    rss_r, _ = rss(Xr)
    df1, df2 = 2, n - k_sub - Xf.shape[1]
    F = ((rss_r - rss_f) / df1) / (rss_f / df2)
    p = stats.f.sf(F, df1, df2)
    amp = np.sqrt(Bf[0] ** 2 + Bf[1] ** 2)
    phase = (np.arctan2(Bf[1], Bf[0]) * PERIOD / (2 * np.pi)) % PERIOD
    return p, amp, phase, df2


def bh(p):
    p = np.asarray(p); o = np.argsort(p); r = np.empty_like(o); r[o] = np.arange(1, len(p) + 1)
    q = p * len(p) / r
    q_sorted = np.minimum.accumulate(q[o][::-1])[::-1]
    out = np.empty_like(q); out[o] = np.minimum(q_sorted, 1.0)
    return out


def main():
    meta, expr = read_series(SERIES)
    if "subject" not in meta.columns:
        meta["subject"] = np.nan
    need = ["geo_accession", "title", "subject", "sleepprotocol", "hoursawake", "circadianphase"]
    miss = [c for c in need if c not in meta.columns]
    if miss:
        raise SystemExit(f"樣本資訊缺欄位：{miss}；實際欄位：{list(meta.columns)}")
    import re
    t_subj = meta["title"].astype(str).str.strip().str.extract(r"^(?:Subject_)?([A-Za-z0-9]+)_S[ER]_T\d+$", expand=False)
    if t_subj.isna().any():
        raise SystemExit(f"有 {int(t_subj.isna().sum())} 個樣本的標題取不到受試者編號：{meta.loc[t_subj.isna(), 'title'].head().tolist()}")
    both = meta["subject"].notna()
    bad = both & (meta["subject"].astype(str).str.strip() != t_subj)
    if bad.any():
        raise SystemExit(f"標題與特徵欄的受試者編號不一致：{meta.loc[bad, ['title', 'subject']].head().to_string()}")
    print(f"特徵欄缺受試者編號的樣本：{int((~both).sum())} 個，已由標題補上")
    meta["subject"] = t_subj
    meta["circadianphase"] = pd.to_numeric(meta["circadianphase"], errors="coerce")
    meta["hoursawake"] = pd.to_numeric(meta["hoursawake"], errors="coerce")
    print(f"樣本 {len(meta)}；表現值 {expr.shape[0]:,} 探針 × {expr.shape[1]} 樣本")
    per = meta.groupby(["sleepprotocol", "subject"]).size()
    print(f"每位受試者每種條件的樣本數：最少 {per.min()}、最多 {per.max()}")
    print(meta.groupby("sleepprotocol").agg(樣本=("geo_accession", "size"), 受試者=("subject", "nunique"),
                                            相位最小=("circadianphase", "min"), 相位最大=("circadianphase", "max")).to_string())

    ann = read_platform(PLATFORM)
    ann = ann[(ann["CONTROL_TYPE"].str.upper() == "FALSE") & ann["GENE_SYMBOL"].notna() & (ann["GENE_SYMBOL"].str.strip() != "")]
    ann["GENE_SYMBOL"] = ann["GENE_SYMBOL"].str.strip().str.upper()
    expr = expr.loc[expr.index.intersection(ann["ID"])]
    sym = ann.set_index("ID").loc[expr.index, "GENE_SYMBOL"]
    means = expr.mean(axis=1)
    best = means.groupby(sym).idxmax()
    G = expr.loc[best.values]; G.index = best.index
    print(f"排除控制探針並合併後：{G.shape[0]:,} 個基因")

    sc = pd.read_csv(PIPE + "step_final_scores.csv")
    sc["gene"] = sc["gene_symbol"].astype(str).str.upper()
    cand = set(pd.read_csv(PIPE + "step_final_candidates.csv")["gene_symbol"].astype(str).str.upper())
    stab = pd.read_csv(PIPE + "cand_stability_genes.tsv", sep="\t")
    core = set(stab.loc[stab["sel_total"] >= 36, "gene_symbol"].astype(str).str.upper())
    freq75 = set(stab.loc[stab["sel_total"] >= 30, "gene_symbol"].astype(str).str.upper())
    anchors = set(sc.loc[sc["is_anchor"] == 1, "gene"])
    train_pos = set(sc.loc[(sc["is_hi_amp"] == 1) & (sc["is_anchor"] == 0), "gene"])
    circadb = set(sc.loc[sc["is_circadb"] == 1, "gene"])
    pool = set(sc.loc[sc["annotation_confidence_level"].isin(["LOW", "VERY_LOW"]) & (sc["is_anchor"] == 0)
                      & (sc["is_circadb"] == 0), "gene"]) - cand
    groups = {"候選 399": cand, "穩定核心（≥ 36／40）": core, "入選頻率 ≥ 30／40": freq75, "候選中的非核心": cand - core,
              "候選中的 IEG": cand & set(IEG), "候選中的非 IEG": cand - set(IEG),
              "對照池（低文獻、非 CircaDB、非候選）": pool, "全部 TCGA 分數基因": set(sc["gene"]),
              "19 個 anchor（正對照）": anchors, "CircaDB 訓練正例 164": train_pos, "CircaDB 全部": circadb}

    summary = []
    conds = [c for c in sorted(meta["sleepprotocol"].dropna().unique()) if c.lower().startswith("sleep")]
    print(f"分析條件：{conds}")
    for cond in conds:
        m = meta[(meta["sleepprotocol"] == cond) & meta["circadianphase"].notna() & meta["hoursawake"].notna()]
        Y = G[m["geo_accession"]].T
        Y = Y.loc[:, Y.notna().all(axis=0)]
        p, amp, ph, df2 = rhythm_test(Y.to_numpy(float), m["circadianphase"].to_numpy(float),
                                      m["hoursawake"].to_numpy(float), m["subject"].to_numpy())
        res = pd.DataFrame({"gene": Y.columns, "p": p, "q": bh(p), "amplitude_log2": amp, "peak_phase_h": ph})
        res["rhythmic"] = res["q"] < Q_CUT
        res["rhythmic_amp"] = res["rhythmic"] & (res["amplitude_log2"] >= AMP_CUT)
        tag = cond.replace(" ", "_")
        res.to_csv(f"ts_GSE39445_gene_stats_{tag}.tsv", sep="\t", index=False)
        print(f"\n=== {cond}：{len(m)} 個樣本、{m['subject'].nunique()} 位受試者、{len(res):,} 個基因（殘差自由度 {df2}）===")
        print(f"  全體有節律（q < {Q_CUT}）：{int(res.rhythmic.sum()):,} 個（{res.rhythmic.mean():.1%}）；"
              f"且振幅 ≥ {AMP_CUT}：{int(res.rhythmic_amp.sum()):,} 個")
        res["q001"] = res["q"] < 0.001
        measured = set(res["gene"]); pool_m = pool & measured
        nlp = pd.Series(-np.log10(np.clip(res["p"].to_numpy(), 1e-300, 1)), index=res["gene"])
        for crit, col in [("q < 0.05", "rhythmic"), ("q < 0.05 且振幅 ≥ 0.1", "rhythmic_amp"), ("q < 0.001", "q001")]:
            print(f"  ── 門檻：{crit}（全體 {res[col].mean():.1%}）")
            rset = set(res.loc[res[col], "gene"])
            for gname, gs in groups.items():
                gm = gs & measured
                k = len(gm & rset); frac = k / len(gm) if gm else np.nan
                row = {"condition": cond, "criterion": crit, "group": gname, "n_measured": len(gm), "n_rhythmic": k, "frac_rhythmic": frac}
                if gname not in ("對照池（低文獻、非 CircaDB、非候選）", "全部 TCGA 分數基因") and gm and pool_m:
                    c = len(pool_m & rset)
                    orr, fp = stats.fisher_exact([[k, len(gm) - k], [c, len(pool_m) - c]], alternative="greater")
                    u = stats.mannwhitneyu(nlp[list(gm)], nlp[list(pool_m)], alternative="greater")
                    row.update({"odds_ratio_vs_pool": orr, "fisher_p_vs_pool": fp, "mwu_p_vs_pool": u.pvalue})
                summary.append(row)
                extra = (f"；OR {row['odds_ratio_vs_pool']:.2f}，Fisher p = {row['fisher_p_vs_pool']:.1e}，排名比較 p = {row['mwu_p_vs_pool']:.1e}"
                         if "fisher_p_vs_pool" in row else "")
                print(f"    {gname:<22} {k:>4}／{len(gm):<5}（{frac:.1%}）{extra}")
        continue
        rset = set(res.loc[res.rhythmic, "gene"]); measured = set(res["gene"])
        pool_m = pool & measured
        for gname, gs in groups.items():
            gm = gs & measured
            k = len(gm & rset)
            frac = k / len(gm) if gm else np.nan
            row = {"condition": cond, "group": gname, "n_measured": len(gm), "n_rhythmic": k, "frac_rhythmic": frac,
                   "median_amp": res.loc[res.gene.isin(gm), "amplitude_log2"].median() if gm else np.nan}
            if gname not in ("對照池（低文獻、非 CircaDB、非候選）", "全部 TCGA 分數基因") and gm and pool_m:
                a, b = k, len(gm) - k
                c = len(pool_m & rset); d = len(pool_m) - c
                orr, fp = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
                row.update({"odds_ratio_vs_pool": orr, "fisher_p_vs_pool": fp})
            summary.append(row)
            extra = f"；對照池比較 OR {row['odds_ratio_vs_pool']:.2f}，單尾 p = {row['fisher_p_vs_pool']:.2e}" if "fisher_p_vs_pool" in row else ""
            print(f"  {gname:<22} 測得 {len(gm):>5}，有節律 {k:>4}（{frac:.1%}），振幅中位數 {row['median_amp']:.3f}{extra}")
    pd.DataFrame(summary).to_csv("ts_GSE39445_summary.tsv", sep="\t", index=False)
    print("\n已輸出 ts_GSE39445_gene_stats_<條件>.tsv、ts_GSE39445_summary.tsv")


if __name__ == "__main__":
    main()