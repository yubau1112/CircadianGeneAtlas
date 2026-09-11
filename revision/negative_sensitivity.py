#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CGA 第 4 項：負樣本敏感度分析（另含第 5 項的特徵對照）
流程與 cga_ms3_robust.py 的 B 段完全相同：
  正例 = CircaDB（FDR ≤ 0.05）每基因最大 rAMP 前 40%，排除 anchor → 164 個
  模型 = 中位數補值 → 標準化 → L1 邏輯迴歸（liblinear, class_weight=balanced）
  五摺分層交叉驗證；種子 0–9；AUC 取平均 ± 標準差（與 robust.py 相同用 ddof=0）
只換負例來源；另外用種子 42 訓練最終模型、照原規則挑候選，與 step_final_candidates.csv 比重疊。

執行位置：/home/yubau/aging/20250519_circadian/Phase2/circadian_pipeline
輸出：neg_sensitivity_metrics.tsv、neg_sensitivity_aucs.tsv、neg_sensitivity_candidates.tsv、neg_train_seed42.txt
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score

COEXPR = "coexpr_features_adj.csv"; METH = "meth_features_adj.csv"; CIRCADB = "CircaDB_HA_query.csv"
FINAL_CAND = "step_final_candidates.csv"
LVL, FDR, AMP_FRAC, SEEDS, RS_FINAL, AI_PCT = "annotation_confidence_level", 0.05, 0.40, range(10), 42, 0.95
LOWLIT = {"VERY_LOW", "LOW"}


def build_df():
    """與 cga_ms3_robust.py 的 build_df 相同"""
    co = pd.read_csv(COEXPR); me = pd.read_csv(METH)
    coex = [c for c in co.columns if c.startswith("radj_")]
    meth = [c for c in me.columns if c.startswith("meth_radj_")]
    df = co.merge(me[["gene_symbol"] + meth], on="gene_symbol", how="left")
    R = df[coex].to_numpy(float)
    POS = ["radj_ARNTL", "radj_CLOCK", "radj_NPAS2", "radj_RORA", "radj_RORB", "radj_RORC", "radj_ARNTL2"]
    NEG = ["radj_PER1", "radj_PER2", "radj_PER3", "radj_CRY1", "radj_CRY2", "radj_NR1D1", "radj_NR1D2", "radj_BHLHE40", "radj_BHLHE41"]
    df["d_posneg"] = np.nanmean(df[[c for c in POS if c in coex]], axis=1) - np.nanmean(df[[c for c in NEG if c in coex]], axis=1)
    df["d_mean"] = np.nanmean(R, axis=1); df["d_max"] = np.nanmax(R, axis=1); df["d_sign"] = np.nanmean(np.sign(R), axis=1)
    deriv = ["d_posneg", "d_mean", "d_max", "d_sign"]
    return df, coex, deriv, meth


def make_pipe():
    return Pipeline([("i", SimpleImputer(strategy="median")), ("s", StandardScaler()),
                     ("m", LogisticRegression(max_iter=3000, class_weight="balanced", penalty="l1", solver="liblinear"))])


def cv_auc(df, cols, pos, negp, seed):
    neg = negp.sample(n=len(pos), random_state=seed)
    tr = pd.concat([pos, neg]); X = tr[cols].to_numpy(float)
    y = np.r_[np.ones(len(pos), int), np.zeros(len(neg), int)]
    cv = StratifiedKFold(5, shuffle=True, random_state=seed)
    oof = cross_val_predict(make_pipe(), X, y, cv=cv, method="predict_proba", n_jobs=-1)[:, 1]
    return roc_auc_score(y, oof)


def final_candidates(df, cols, pos, negp, ia, circ):
    """與 cga_ms3_final_candidates.py 相同：種子 42 抽負例、訓練、全基因打分、挑候選"""
    neg = negp.sample(n=len(pos), random_state=RS_FINAL)
    tr = pd.concat([pos, neg]); y = np.r_[np.ones(len(pos), int), np.zeros(len(neg), int)]
    pipe = make_pipe().fit(tr[cols].to_numpy(float), y)
    score = pipe.predict_proba(df[cols].to_numpy(float))[:, 1]
    cut = np.quantile(score, AI_PCT)
    keep = (score >= cut) & df[LVL].isin(LOWLIT).to_numpy() & (~ia).to_numpy() & (~df["gene_symbol"].isin(circ)).to_numpy()
    return set(df.loc[keep, "gene_symbol"]), set(neg["gene_symbol"])


def main():
    df, coex, deriv, meth = build_df()
    allf = coex + deriv + meth
    cd = pd.read_csv(CIRCADB); circ_all = cd[cd["FDR"] <= FDR].copy()
    circ_all["Gene_Symbol"] = circ_all["Gene_Symbol"].astype(str).str.strip()
    circ = set(circ_all["Gene_Symbol"])
    amp = circ_all.groupby("Gene_Symbol")["rAMP"].max()
    hi = set(amp[amp >= amp.quantile(1 - AMP_FRAC)].index)

    ia = df["is_anchor"] == 1
    not_circ = ~df["gene_symbol"].isin(circ)
    pos = df[df["gene_symbol"].isin(hi) & ~ia]
    print(f"基因數 {len(df):,}；CircaDB（FDR ≤ 0.05）{len(circ)}；前 40% {len(hi)}；正例 {len(pos)}")
    if len(pos) != 164:
        raise SystemExit("正例不是 164 個，與論文流程不同，先停止")

    lvl = df[LVL]
    pools = {
        "N1 原本：文獻 VERY_LOW": df[(lvl == "VERY_LOW") & ~ia & not_circ],
        "N2 不看文獻：全部非 CircaDB": df[~ia & not_circ],
        "N3 文獻有提及：LOW／MEDIUM／HIGH": df[lvl.isin({"LOW", "MEDIUM", "HIGH"}) & ~ia & not_circ],
        "N4 文獻多：MEDIUM／HIGH": df[lvl.isin({"MEDIUM", "HIGH"}) & ~ia & not_circ],
    }
    feats = {
        "全部特徵": allf,
        "只用共表達": coex,
        "共表達＋衍生": coex + deriv,
        "共表達＋甲基化": coex + meth,
    }
    runs = [(p, "全部特徵") for p in pools] + [("N1 原本：文獻 VERY_LOW", f) for f in list(feats)[1:]]

    per_seed, summary = [], []
    for pname, fname in runs:
        negp, cols = pools[pname], feats[fname]
        if len(negp) < len(pos):
            print(f"  {pname}：負例池只有 {len(negp)} 個，不足 164，跳過"); continue
        aucs = np.array([cv_auc(df, cols, pos, negp, s) for s in SEEDS])
        for s, a in zip(SEEDS, aucs):
            per_seed.append({"negatives": pname, "features": fname, "seed": s, "AUC": a})
        summary.append({"negatives": pname, "features": fname, "neg_pool": len(negp), "n_features": len(cols),
                        "AUC_mean": aucs.mean(), "AUC_sd": aucs.std(), "AUC_min": aucs.min(), "AUC_max": aucs.max()})
        print(f"  {pname:<24}｜{fname:<8}：AUC {aucs.mean():.3f} ± {aucs.std():.3f}（範圍 {aucs.min():.3f}–{aucs.max():.3f}，負例池 {len(negp):,}）")

    summ = pd.DataFrame(summary)
    base = summ[(summ.negatives.str.startswith("N1")) & (summ.features == "全部特徵")].iloc[0]
    ok = abs(base.AUC_mean - 0.805) < 0.0006 and abs(base.AUC_sd - 0.018) < 0.0006
    print(f"\n重現檢查：N1 全部特徵 = {base.AUC_mean:.3f} ± {base.AUC_sd:.3f}（論文 0.805 ± 0.018）→ {'相符' if ok else '不同，請留意'}")

    # ── 候選重疊（種子 42，與正式版相同）──
    orig = set(pd.read_csv(FINAL_CAND)["gene_symbol"])
    cand_rows = []
    print(f"\n=== 候選重疊（與 {FINAL_CAND} 的 {len(orig)} 個比）===")
    for pname, negp in pools.items():
        if len(negp) < len(pos):
            continue
        cset, negset = final_candidates(df, allf, pos, negp, ia, circ)
        both = cset & orig
        jac = len(both) / len(cset | orig) if (cset | orig) else np.nan
        cand_rows.append({"negatives": pname, "n_candidates": len(cset), "overlap_with_399": len(both), "jaccard": jac})
        print(f"  {pname:<24}：候選 {len(cset)} 個，與原 399 重疊 {len(both)} 個（Jaccard {jac:.3f}）")
        if pname.startswith("N1"):
            same = (cset == orig)
            print(f"  重現檢查：N1 種子 42 挑出的候選與原 399 {'完全相同' if same else '不同，請留意'}")
            with open("neg_train_seed42.txt", "w", encoding="utf-8") as f:
                f.write("\n".join(sorted(negset)) + "\n")
            print(f"  正式模型的 164 個訓練負例已存 neg_train_seed42.txt（{len(negset)} 個）")

    summ.to_csv("neg_sensitivity_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(per_seed).to_csv("neg_sensitivity_aucs.tsv", sep="\t", index=False)
    pd.DataFrame(cand_rows).to_csv("neg_sensitivity_candidates.tsv", sep="\t", index=False)
    print("\n已輸出 neg_sensitivity_metrics.tsv、neg_sensitivity_aucs.tsv、neg_sensitivity_candidates.tsv、neg_train_seed42.txt")


if __name__ == "__main__":
    main()