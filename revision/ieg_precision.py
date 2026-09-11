#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CGA 第 12 項（IEG 驗證）＋ 第 7 項（候選 precision 估計）
流程沿用 cga_ms3_final_candidates.py / cga_ms3_robust.py：
  正例 = CircaDB 前 40% 高振幅（164）；負例 = 文獻 VERY_LOW、非 CircaDB；L1 邏輯迴歸；五摺交叉驗證

A. 候選與各 anchor 的共表達輪廓（radj_*）：候選、訓練正例、全體基因各自的平均
B. 拿掉 BHLHE40 的所有特徵（radj_BHLHE40、meth_radj_BHLHE40，衍生特徵也不含 BHLHE40）後重跑：
   AUC（種子 0–9）、種子 42 候選與原 399 的重疊、IEG 在候選中的比例
C. precision 估計：以種子 42 正式模型的全基因組第 95 百分位分數為門檻，
   用種子 0–9 的交叉驗證（out-of-fold）機率算 TPR 與 FPR，再依不同的節律基因比例假設換算 precision

IEG 清單為本分析暫定的常見立即早期基因，寫進論文前要換成有引用出處的清單。
執行位置：/home/yubau/aging/20250519_circadian/Phase2/circadian_pipeline
輸出：ieg_anchor_profile.tsv、ieg_drop_bhlhe40.tsv、precision_estimate.tsv
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
LVL, FDR, AMP_FRAC, AI_PCT, RS_FINAL, SEEDS = "annotation_confidence_level", 0.05, 0.40, 0.95, 42, range(10)
LOWLIT = {"VERY_LOW", "LOW"}
IEG = ["FOS", "FOSB", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND", "EGR1", "EGR2", "EGR3", "EGR4",
       "IER2", "IER3", "IER5", "NR4A1", "NR4A2", "NR4A3", "ATF3", "DUSP1", "DUSP2", "DUSP5", "DUSP6",
       "KLF2", "KLF4", "KLF6", "KLF10", "ZFP36", "ZFP36L1", "ZFP36L2", "CCN1", "CCN2", "BTG2",
       "PPP1R15A", "GADD45B", "GADD45G", "MYC", "ARC", "NPAS4", "SGK1", "CEBPB", "CEBPD", "SOCS3",
       "HES1", "PLK2", "CSRNP1", "NFKBIA", "IRF1", "TNFAIP3", "KDM6B", "SIK1"]
POS_ARM = ["ARNTL", "CLOCK", "NPAS2", "RORA", "RORB", "RORC", "ARNTL2"]
NEG_ARM = ["PER1", "PER2", "PER3", "CRY1", "CRY2", "NR1D1", "NR1D2", "BHLHE40", "BHLHE41"]


def build(exclude=()):
    """與原腳本相同的特徵；exclude 裡的 anchor 其共表達、甲基化與衍生特徵全部不用"""
    co = pd.read_csv(COEXPR); me = pd.read_csv(METH)
    coex = [c for c in co.columns if c.startswith("radj_") and c[5:] not in exclude]
    meth = [c for c in me.columns if c.startswith("meth_radj_") and c[10:] not in exclude]
    df = co.merge(me[["gene_symbol"] + [c for c in me.columns if c.startswith("meth_radj_")]], on="gene_symbol", how="left")
    R = df[coex].to_numpy(float)
    pa = [f"radj_{g}" for g in POS_ARM if f"radj_{g}" in coex]
    na = [f"radj_{g}" for g in NEG_ARM if f"radj_{g}" in coex]
    df["d_posneg"] = np.nanmean(df[pa], axis=1) - np.nanmean(df[na], axis=1)
    df["d_mean"] = np.nanmean(R, axis=1); df["d_max"] = np.nanmax(R, axis=1); df["d_sign"] = np.nanmean(np.sign(R), axis=1)
    return df, coex + ["d_posneg", "d_mean", "d_max", "d_sign"] + meth


def pipe():
    return Pipeline([("i", SimpleImputer(strategy="median")), ("s", StandardScaler()),
                     ("m", LogisticRegression(max_iter=3000, class_weight="balanced", penalty="l1", solver="liblinear"))])


def main():
    cd = pd.read_csv(CIRCADB); ca = cd[cd["FDR"] <= FDR].copy()
    ca["Gene_Symbol"] = ca["Gene_Symbol"].astype(str).str.strip()
    circ = set(ca["Gene_Symbol"]); amp = ca.groupby("Gene_Symbol")["rAMP"].max()
    hi = set(amp[amp >= amp.quantile(1 - AMP_FRAC)].index)
    orig = set(pd.read_csv(FINAL_CAND)["gene_symbol"])

    def setup(df):
        ia = df["is_anchor"] == 1
        pos = df[df["gene_symbol"].isin(hi) & ~ia]
        negp = df[(df[LVL] == "VERY_LOW") & ~ia & ~df["gene_symbol"].isin(circ)]
        elig = (df[LVL].isin(LOWLIT) & ~ia & ~df["gene_symbol"].isin(circ)).to_numpy()
        return pos, negp, elig

    def final_run(df, cols):
        pos, negp, elig = setup(df)
        neg = negp.sample(n=len(pos), random_state=RS_FINAL)
        tr = pd.concat([pos, neg]); y = np.r_[np.ones(len(pos), int), np.zeros(len(neg), int)]
        m = pipe().fit(tr[cols].to_numpy(float), y)
        score = m.predict_proba(df[cols].to_numpy(float))[:, 1]
        cut = np.quantile(score, AI_PCT)
        return set(df.loc[(score >= cut) & elig, "gene_symbol"]), cut

    def oof(df, cols, seed):
        pos, negp, _ = setup(df)
        neg = negp.sample(n=len(pos), random_state=seed)
        tr = pd.concat([pos, neg]); y = np.r_[np.ones(len(pos), int), np.zeros(len(neg), int)]
        p = cross_val_predict(pipe(), tr[cols].to_numpy(float), y, cv=StratifiedKFold(5, shuffle=True, random_state=seed),
                              method="predict_proba", n_jobs=-1)[:, 1]
        return y, p

    df, cols = build()
    pos, _, _ = setup(df)
    if len(pos) != 164:
        raise SystemExit("正例不是 164 個，先停止")
    cand, cut = final_run(df, cols)
    print(f"重現檢查：種子 42 候選 {len(cand)} 個，與原 399 {'完全相同' if cand == orig else '不同，請留意'}；第 95 百分位門檻 = {cut:.4f}")

    # ── A. 共表達輪廓 ──
    anchors = [c[5:] for c in df.columns if c.startswith("radj_")]
    grp = {"候選 399": df["gene_symbol"].isin(orig), "訓練正例 164": df["gene_symbol"].isin(set(pos["gene_symbol"])),
           "全體基因": pd.Series(True, index=df.index)}
    prof = pd.DataFrame({g: [df.loc[m, f"radj_{a}"].mean() for a in anchors] for g, m in grp.items()}, index=anchors)
    prof["候選減全體"] = prof["候選 399"] - prof["全體基因"]
    prof = prof.sort_values("候選 399", ascending=False)
    prof.to_csv("ieg_anchor_profile.tsv", sep="\t")
    print("\n=== A. 各 anchor 的平均共表達（radj，已校正）：候選 vs 訓練正例 vs 全體 ===")
    print(prof.round(3).to_string())

    ieg_in = sorted(set(IEG) & orig)
    print(f"\n原 399 候選中的 IEG（暫定清單 {len(IEG)} 個）：{len(ieg_in)} 個 → {', '.join(ieg_in)}")

    # ── B. 拿掉 BHLHE40 ──
    df2, cols2 = build(exclude=("BHLHE40",))
    a_full = [roc_auc_score(*oof(df, cols, s)) for s in SEEDS]
    a_drop = [roc_auc_score(*oof(df2, cols2, s)) for s in SEEDS]
    cand2, _ = final_run(df2, cols2)
    ieg2 = sorted(set(IEG) & cand2)
    rowsB = [{"model": "全部特徵（正式）", "n_features": len(cols), "AUC_mean": np.mean(a_full), "AUC_sd": np.std(a_full),
              "n_candidates": len(cand), "overlap_with_399": len(cand & orig), "n_IEG_in_candidates": len(ieg_in)},
             {"model": "拿掉 BHLHE40", "n_features": len(cols2), "AUC_mean": np.mean(a_drop), "AUC_sd": np.std(a_drop),
              "n_candidates": len(cand2), "overlap_with_399": len(cand2 & orig), "n_IEG_in_candidates": len(ieg2)}]
    pd.DataFrame(rowsB).to_csv("ieg_drop_bhlhe40.tsv", sep="\t", index=False)
    print("\n=== B. 拿掉 BHLHE40 的所有特徵後重跑 ===")
    print(f"  全部特徵：AUC {np.mean(a_full):.3f} ± {np.std(a_full):.3f}（{len(cols)} 個特徵）；候選 {len(cand)}；IEG {len(ieg_in)} 個")
    print(f"  拿掉 BHLHE40：AUC {np.mean(a_drop):.3f} ± {np.std(a_drop):.3f}（{len(cols2)} 個特徵）；候選 {len(cand2)}，"
          f"與原 399 重疊 {len(cand2 & orig)}；IEG {len(ieg2)} 個 → {', '.join(ieg2)}")

    # ── C. precision 估計 ──
    tprs, fprs = [], []
    for s in SEEDS:
        y, p = oof(df, cols, s)
        tprs.append(np.mean(p[y == 1] >= cut)); fprs.append(np.mean(p[y == 0] >= cut))
    tpr, fpr = np.mean(tprs), np.mean(fprs)
    print(f"\n=== C. 候選 precision 估計（門檻 = 第 95 百分位 {cut:.4f}；種子 0–9 的 out-of-fold 機率）===")
    print(f"  正例超過門檻的比例（TPR）= {tpr:.3f}（範圍 {min(tprs):.3f}–{max(tprs):.3f}）")
    print(f"  負例超過門檻的比例（FPR）= {fpr:.3f}（範圍 {min(fprs):.3f}–{max(fprs):.3f}）")
    rowsC = []
    for pi in (0.02, 0.05, 0.10, 0.20, 0.30):
        prec = pi * tpr / (pi * tpr + (1 - pi) * fpr) if (pi * tpr + (1 - pi) * fpr) > 0 else np.nan
        rowsC.append({"assumed_prevalence": pi, "TPR": tpr, "FPR": fpr, "precision": prec,
                      "expected_false_in_399": 399 * (1 - prec)})
        print(f"  若背景基因中有 {pi:>4.0%} 屬於高振幅節律基因：precision ≈ {prec:.3f}，399 個中預期假陽性 ≈ {399 * (1 - prec):.0f} 個")
    pd.DataFrame(rowsC).to_csv("precision_estimate.tsv", sep="\t", index=False)
    print("\n已輸出 ieg_anchor_profile.tsv、ieg_drop_bhlhe40.tsv、precision_estimate.tsv")


if __name__ == "__main__":
    main()