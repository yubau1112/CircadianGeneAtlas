#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CGA 第 4 項延伸：候選穩定度分析
每種負例 × 種子 0–9，照 cga_ms3_final_candidates.py 的流程訓練、全基因打分、挑候選
（ai_score 前 5% ∩ 文獻 LOW／VERY_LOW ∩ 非 anchor ∩ 非 CircaDB），統計每個基因的入選頻率。

執行位置：/home/yubau/aging/20250519_circadian/Phase2/circadian_pipeline
輸出：cand_stability_genes.tsv（每個基因的入選次數）、cand_stability_runs.tsv（每次的候選數與重疊）
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer

COEXPR = "coexpr_features_adj.csv"; METH = "meth_features_adj.csv"; CIRCADB = "CircaDB_HA_query.csv"
FINAL_CAND = "step_final_candidates.csv"
LVL, FDR, AMP_FRAC, AI_PCT, SEEDS = "annotation_confidence_level", 0.05, 0.40, 0.95, list(range(10))
LOWLIT = {"VERY_LOW", "LOW"}
TABLE8 = ["TRIB1", "CPEB4", "CCN2", "RGS2", "PPP1R3B"]


def build_df():
    """與 cga_ms3_final_candidates.py 相同"""
    co = pd.read_csv(COEXPR); me = pd.read_csv(METH)
    coex = [c for c in co.columns if c.startswith("radj_")]
    meth = [c for c in me.columns if c.startswith("meth_radj_")]
    df = co.merge(me[["gene_symbol"] + meth], on="gene_symbol", how="left")
    R = df[coex].to_numpy(float)
    POS = ["radj_ARNTL", "radj_CLOCK", "radj_NPAS2", "radj_RORA", "radj_RORB", "radj_RORC", "radj_ARNTL2"]
    NEG = ["radj_PER1", "radj_PER2", "radj_PER3", "radj_CRY1", "radj_CRY2", "radj_NR1D1", "radj_NR1D2", "radj_BHLHE40", "radj_BHLHE41"]
    df["d_posneg"] = np.nanmean(df[[c for c in POS if c in coex]], axis=1) - np.nanmean(df[[c for c in NEG if c in coex]], axis=1)
    df["d_mean"] = np.nanmean(R, axis=1); df["d_max"] = np.nanmax(R, axis=1); df["d_sign"] = np.nanmean(np.sign(R), axis=1)
    return df, coex + ["d_posneg", "d_mean", "d_max", "d_sign"] + meth


def main():
    df, cols = build_df()
    X_all = df[cols].to_numpy(float)
    cd = pd.read_csv(CIRCADB); circ_all = cd[cd["FDR"] <= FDR].copy()
    circ_all["Gene_Symbol"] = circ_all["Gene_Symbol"].astype(str).str.strip()
    circ = set(circ_all["Gene_Symbol"])
    amp = circ_all.groupby("Gene_Symbol")["rAMP"].max()
    hi = set(amp[amp >= amp.quantile(1 - AMP_FRAC)].index)
    ia = df["is_anchor"] == 1
    not_circ = ~df["gene_symbol"].isin(circ)
    pos = df[df["gene_symbol"].isin(hi) & ~ia]
    if len(pos) != 164:
        raise SystemExit("正例不是 164 個，先停止")
    eligible = (df[LVL].isin(LOWLIT) & ~ia & not_circ).to_numpy()
    lvl = df[LVL]
    pools = {
        "N1": df[(lvl == "VERY_LOW") & ~ia & not_circ],
        "N2": df[~ia & not_circ],
        "N3": df[lvl.isin({"LOW", "MEDIUM", "HIGH"}) & ~ia & not_circ],
        "N4": df[lvl.isin({"MEDIUM", "HIGH"}) & ~ia & not_circ],
    }
    names = {"N1": "原本：文獻 VERY_LOW", "N2": "不看文獻：全部非 CircaDB",
             "N3": "文獻有提及：LOW／MEDIUM／HIGH", "N4": "文獻多：MEDIUM／HIGH"}

    orig = set(pd.read_csv(FINAL_CAND)["gene_symbol"])
    orig_df = pd.read_csv(FINAL_CAND)
    top20 = list(orig_df.sort_values("ai_score", ascending=False)["gene_symbol"].head(20))
    count = {p: np.zeros(len(df), int) for p in pools}
    runs = []

    def one_run(pname, seed):
        neg = pools[pname].sample(n=len(pos), random_state=seed)
        tr = pd.concat([pos, neg]); y = np.r_[np.ones(len(pos), int), np.zeros(len(neg), int)]
        pipe = Pipeline([("i", SimpleImputer(strategy="median")), ("s", StandardScaler()),
                         ("m", LogisticRegression(max_iter=3000, class_weight="balanced", penalty="l1", solver="liblinear"))])
        pipe.fit(tr[cols].to_numpy(float), y)
        score = pipe.predict_proba(X_all)[:, 1]
        return (score >= np.quantile(score, AI_PCT)) & eligible

    # 基準：N1 種子 42 應重現原 399
    chk = set(df.loc[one_run("N1", 42), "gene_symbol"])
    print(f"重現檢查：N1 種子 42 → {len(chk)} 個候選，與原 399 {'完全相同' if chk == orig else '不同，請留意'}")

    for p in pools:
        for s in SEEDS:
            sel = one_run(p, s)
            count[p] += sel
            cset = set(df.loc[sel, "gene_symbol"])
            runs.append({"negatives": p, "seed": s, "n_candidates": len(cset), "overlap_with_399": len(cset & orig)})
    runs = pd.DataFrame(runs)

    print("\n=== 各負例、種子 0–9：候選數與「與原 399 重疊」===")
    for p in pools:
        r = runs[runs.negatives == p]
        print(f"  {p} {names[p]:<24}：候選 {r.n_candidates.mean():.0f} 個（{r.n_candidates.min()}–{r.n_candidates.max()}）；"
              f"與原 399 重疊 平均 {r.overlap_with_399.mean():.0f} 個（{r.overlap_with_399.min()}–{r.overlap_with_399.max()}）")

    total = sum(count.values()); nrun = len(pools) * len(SEEDS)
    g = pd.DataFrame({"gene_symbol": df["gene_symbol"], "in_original_399": df["gene_symbol"].isin(orig).astype(int)})
    for p in pools:
        g[f"sel_{p}"] = count[p]
    g["sel_total"] = total
    g["sel_frac"] = total / nrun
    g = g[g.sel_total > 0].sort_values(["sel_total", "in_original_399"], ascending=False)
    g.to_csv("cand_stability_genes.tsv", sep="\t", index=False)
    runs.to_csv("cand_stability_runs.tsv", sep="\t", index=False)

    print(f"\n=== 入選頻率（共 {nrun} 次：4 種負例 × 10 個種子）===")
    for k, lab in [(nrun, "每次都入選"), (int(0.9 * nrun), "≥ 90%"), (int(0.75 * nrun), "≥ 75%"), (int(0.5 * nrun), "≥ 50%")]:
        sub = g[g.sel_total >= k]
        print(f"  {lab:<6}（≥ {k} 次）：{len(sub)} 個基因，其中屬於原 399 的 {int(sub.in_original_399.sum())} 個")
    n1 = g[g.sel_N1 >= 9]
    print(f"  只看 N1（原本負例），10 次中 ≥ 9 次入選：{len(n1)} 個，其中屬於原 399 的 {int(n1.in_original_399.sum())} 個")
    never = orig - set(g.gene_symbol)
    print(f"  原 399 中，40 次一次都沒入選的：{len(never)} 個（原 399 本身是種子 42，不在這 40 次裡）")

    def show(glist, title):
        print(f"\n=== {title} ===")
        for sym in glist:
            row = g[g.gene_symbol == sym]
            if len(row) == 0:
                print(f"  {sym:<10} 0 / {nrun}")
            else:
                r = row.iloc[0]
                print(f"  {sym:<10} {int(r.sel_total):>2} / {nrun}（N1 {r.sel_N1}、N2 {r.sel_N2}、N3 {r.sel_N3}、N4 {r.sel_N4}）")
    show(TABLE8, "Table 8 的 5 個基因")
    show(top20, "原本前 20 名候選（Figure 10 STRING 網路）")

    print("\n已輸出 cand_stability_genes.tsv、cand_stability_runs.tsv")


if __name__ == "__main__":
    main()