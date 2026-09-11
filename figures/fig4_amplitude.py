#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CGA MS3 — 圖：正例震盪幅度純度 vs AUC 單調曲線（突破證據）
================================================================
穩健驗證 B 的數據：rAMP 門檻越嚴（正例越純），AUC 越高，且變異小。
此單調關係證明提升由真實節律訊號驅動，非抽樣運氣。適合作論文新圖。
輸出：figures/fig12_amplitude_auc.png
================================================================
"""
import os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTDIR="figures"; os.makedirs(OUTDIR,exist_ok=True)
PURPLE,ORANGE="#3b0764","#f59e0b"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":12,"axes.spines.top":False,
                     "axes.spines.right":False,"figure.dpi":300})

# 穩健驗證 B 實測數據（各門檻：正例數、AUC 平均、標準差）
data=[("Top 40%",164,0.805,0.018),
      ("Top 50%",209,0.800,0.017),
      ("Top 60%",243,0.792,0.011),
      ("Top 70%",284,0.781,0.012),
      ("All",419,0.764,0.009)]
# x 軸用「正例純度」：從嚴(左)到寬(右)，即 40→100%
xpct=[40,50,60,70,100]
labels=[d[0] for d in data]
npos=[d[1] for d in data]
auc=[d[2] for d in data]
err=[d[3] for d in data]

fig,ax=plt.subplots(figsize=(7,5))
ax.errorbar(xpct,auc,yerr=err,fmt="o-",color=PURPLE,mfc=ORANGE,mec=PURPLE,
            mew=1.5,ms=9,lw=2,capsize=4,ecolor=PURPLE,label="Mean ± SD (15 seeds)")
# 基準線
ax.axhline(0.764,ls="--",color="#bbbbbb",lw=1)
ax.text(43,0.767,"baseline (all positives)",fontsize=8.5,va="bottom",ha="left",color="#999999")
# 正例數標註：全部放圖頂端一排，各自對齊 x，完全離開資料點與 error bar
for x,n in zip(xpct,npos):
    ax.text(x,0.828,f"n={n}",ha="center",va="top",fontsize=8.5,color="#555555")
ax.set_xlabel("Positive set: rhythmic amplitude percentile (stricter → looser)")
ax.set_ylabel("Cross-validated ROC-AUC")
ax.set_title("Purer high-amplitude positives yield higher AUC")
ax.set_xticks(xpct); ax.set_xticklabels(labels)
ax.set_ylim(0.74,0.832)
ax.legend(loc="lower left",frameon=True,framealpha=0.9,fontsize=9)
plt.tight_layout(); plt.savefig(f"{OUTDIR}/fig12_amplitude_auc.png"); plt.close()
print("  fig12_amplitude_auc.png")
print("完成。")