#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CGA 第 8 項：以人類骨骼肌時序轉錄體重複驗證（GSE108539，Perrin et al. 2018 eLife）
資料：股外側肌切片，每 4 小時一次、共 24 小時（12、16、20、0、4、8 點）；表現值為 GEO 補充檔的外顯子定量（對數尺度）
設計：
  排除時間點不足 6 個的捐贈者（7 號只有 3 個）
  Ensembl → 基因符號：GENCODE v36 主要註解（與 TCGA 相同版本）；同符號取平均表現最高者
  過濾：平均值 > 0（原作者只分析 log2 RPKM > 0 的基因；本檔單位未明確標示，照實註明）
  模型 A：扣除捐贈者平均後，以 cos、sin（24 小時）做 F 檢定
  模型 B：再加入「距第一次切片的經過時間」線性趨勢（檢查反覆切片的壓力反應）
  比較組別與血液版相同
IEG 清單為暫定的常見立即早期基因。
執行位置：/home/yubau/aging/20250519_circadian/Phase2/timeseries_validation
輸出：ts_GSE108539_gene_stats_<模型>.tsv、ts_GSE108539_summary.tsv
"""
import gzip, re
import numpy as np, pandas as pd
from scipy import stats

COUNTS = "GSE108539_norm_reads_count.txt.gz"
GTF = "/home/yubau/CMUH_HCC_RNA-Seq/gencode.v36.annotation.gtf"
PIPE = "../circadian_pipeline/"
PERIOD, MIN_TP, EXPR_MIN = 24.0, 6, 0.0
IEG = ["FOS", "FOSB", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND", "EGR1", "EGR2", "EGR3", "EGR4",
       "IER2", "IER3", "IER5", "NR4A1", "NR4A2", "NR4A3", "ATF3", "DUSP1", "DUSP2", "DUSP5", "DUSP6",
       "KLF2", "KLF4", "KLF6", "KLF10", "ZFP36", "ZFP36L1", "ZFP36L2", "CCN1", "CCN2", "BTG2",
       "PPP1R15A", "GADD45B", "GADD45G", "MYC", "ARC", "NPAS4", "SGK1", "CEBPB", "CEBPD", "SOCS3",
       "HES1", "PLK2", "CSRNP1", "NFKBIA", "IRF1", "TNFAIP3", "KDM6B", "SIK1"]


def gtf_map(path):
    m = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split("\t", 9)
            if len(p) < 9 or p[2] != "gene":
                continue
            gid = re.search(r'gene_id "([^"]+)"', p[8]); gname = re.search(r'gene_name "([^"]+)"', p[8])
            if gid and gname:
                m[gid.group(1).split(".")[0]] = gname.group(1).upper()
    return m


def within(M, subj):
    M = M.copy()
    for s in np.unique(subj):
        i = subj == s
        M[i] -= M[i].mean(axis=0)
    return M


def ftest(Y, t, subj, trend=None):
    w = 2 * np.pi * t / PERIOD
    cols = [np.cos(w), np.sin(w)] + ([trend] if trend is not None else [])
    Xf = within(np.column_stack(cols), subj)
    Yc = within(Y, subj)
    def rss(X):
        if X.shape[1] == 0:
            return (Yc ** 2).sum(axis=0), None
        B, *_ = np.linalg.lstsq(X, Yc, rcond=None)
        return ((Yc - X @ B) ** 2).sum(axis=0), B
    rf, Bf = rss(Xf); rr, _ = rss(Xf[:, 2:])
    df2 = len(t) - len(np.unique(subj)) - Xf.shape[1]
    F = ((rr - rf) / 2) / (rf / df2)
    return stats.f.sf(F, 2, df2), np.sqrt(Bf[0] ** 2 + Bf[1] ** 2), df2


def bh(p):
    p = np.asarray(p); o = np.argsort(p); n = len(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(q, 1); return out


def main():
    df = pd.read_csv(COUNTS, sep="\t", index_col=0)
    ex = [c for c in df.columns if c.startswith("Exon_")]
    info = pd.DataFrame({"col": ex})
    parsed = info["col"].str.extract(r"^Exon_(\d+)_S(\d+)_(\d+)\.(\d+)$")
    if parsed.isna().any().any():
        raise SystemExit(f"欄名格式不符：{info.loc[parsed.isna().any(axis=1), 'col'].tolist()[:5]}")
    info["donor"] = parsed[0]; info["clock"] = parsed[2].astype(int) + parsed[3].astype(int) / 60
    info["elapsed"] = (info["clock"] - 12) % 24
    n_tp = info.groupby("donor").size()
    keep = n_tp[n_tp >= MIN_TP].index
    print(f"外顯子欄 {len(ex)}；捐贈者 {len(n_tp)}，各自時間點數：{n_tp.to_dict()}")
    print(f"保留 ≥ {MIN_TP} 個時間點的捐贈者 {len(keep)} 位；排除：{sorted(set(n_tp.index) - set(keep))}")
    info = info[info["donor"].isin(keep)].reset_index(drop=True)

    E = df[info["col"]]
    E.index = E.index.astype(str).str.split(".").str[0]
    gm = gtf_map(GTF)
    sym = E.index.map(lambda g: gm.get(g, np.nan))
    print(f"Ensembl 編號 {len(E):,}；對到 GENCODE v36 符號 {int(pd.notna(sym).sum()):,}")
    E = E[pd.notna(sym)].copy(); E["__sym"] = np.asarray(sym[pd.notna(sym)])
    E["__mean"] = E[info["col"]].mean(axis=1)
    E = E[E["__mean"] > EXPR_MIN].sort_values("__mean", ascending=False).drop_duplicates("__sym")
    E.index = E["__sym"].to_numpy(); E = E.drop(columns=["__sym", "__mean"])
    print(f"平均值 > {EXPR_MIN} 且合併重複符號後：{len(E):,} 個基因")

    sc = pd.read_csv(PIPE + "step_final_scores.csv"); sc["gene"] = sc["gene_symbol"].astype(str).str.upper()
    cand = set(pd.read_csv(PIPE + "step_final_candidates.csv")["gene_symbol"].astype(str).str.upper())
    stab = pd.read_csv(PIPE + "cand_stability_genes.tsv", sep="\t")
    core = set(stab.loc[stab["sel_total"] >= 36, "gene_symbol"].astype(str).str.upper())
    f75 = set(stab.loc[stab["sel_total"] >= 30, "gene_symbol"].astype(str).str.upper())
    anchors = set(sc.loc[sc["is_anchor"] == 1, "gene"])
    trpos = set(sc.loc[(sc["is_hi_amp"] == 1) & (sc["is_anchor"] == 0), "gene"])
    circ = set(sc.loc[sc["is_circadb"] == 1, "gene"])
    pool = set(sc.loc[sc["annotation_confidence_level"].isin(["LOW", "VERY_LOW"]) & (sc["is_anchor"] == 0) & (sc["is_circadb"] == 0), "gene"]) - cand
    groups = {"候選 399": cand, "穩定核心（≥ 36／40）": core, "入選頻率 ≥ 30／40": f75, "候選中的非核心": cand - core,
              "候選中的 IEG": cand & set(IEG), "候選中的非 IEG": cand - set(IEG), "對照池（低文獻、非 CircaDB、非候選）": pool,
              "19 個 anchor（正對照）": anchors, "CircaDB 訓練正例 164": trpos, "CircaDB 全部": circ}

    Y = E.T.to_numpy(float)
    t = info["clock"].to_numpy(float); subj = info["donor"].to_numpy(); el = info["elapsed"].to_numpy(float)
    summary = []
    for mname, trend in [("A 標準 cosinor", None), ("B 加切片經過時間趨勢", el)]:
        p, amp, df2 = ftest(Y, t, subj, trend)
        res = pd.DataFrame({"gene": E.index, "p": p, "q": bh(p), "amplitude": amp})
        res.to_csv(f"ts_GSE108539_gene_stats_{mname[0]}.tsv", sep="\t", index=False)
        nlp = pd.Series(-np.log10(np.clip(p, 1e-300, 1)), index=E.index)
        measured = set(E.index); pool_m = pool & measured
        print(f"\n=== 模型 {mname}（{len(t)} 個樣本、{len(np.unique(subj))} 位捐贈者，殘差自由度 {df2}）===")
        for crit, mask in [("q < 0.05", res["q"] < 0.05), ("p < 0.01", res["p"] < 0.01)]:
            rset = set(res.loc[mask, "gene"])
            print(f"  ── 門檻：{crit}（全體 {mask.mean():.1%}）")
            for g, gs in groups.items():
                gmz = gs & measured; k = len(gmz & rset); frac = k / len(gmz) if gmz else np.nan
                row = {"model": mname, "criterion": crit, "group": g, "n_measured": len(gmz), "n_rhythmic": k, "frac_rhythmic": frac}
                if g != "對照池（低文獻、非 CircaDB、非候選）" and gmz and pool_m:
                    c = len(pool_m & rset)
                    orr, fp = stats.fisher_exact([[k, len(gmz) - k], [c, len(pool_m) - c]], alternative="greater")
                    u = stats.mannwhitneyu(nlp[list(gmz)], nlp[list(pool_m)], alternative="greater")
                    row.update({"odds_ratio_vs_pool": orr, "fisher_p_vs_pool": fp, "mwu_p_vs_pool": u.pvalue})
                summary.append(row)
                extra = f"；OR {row['odds_ratio_vs_pool']:.2f}，Fisher p = {row['fisher_p_vs_pool']:.1e}，排名比較 p = {row['mwu_p_vs_pool']:.1e}" if "fisher_p_vs_pool" in row else ""
                print(f"    {g:<22} {k:>4}／{len(gmz):<5}（{frac:.1%}）{extra}")
    pd.DataFrame(summary).to_csv("ts_GSE108539_summary.tsv", sep="\t", index=False)
    print("\n已輸出 ts_GSE108539_gene_stats_<模型>.tsv、ts_GSE108539_summary.tsv")


if __name__ == "__main__":
    main()