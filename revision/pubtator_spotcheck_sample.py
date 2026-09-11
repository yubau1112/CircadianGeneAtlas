#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CGA 第 10 項：PubTator3 基因抽取的人工抽查樣本（回應審稿 R2 缺點 5、建議 3）
抽樣單位：一篇文章 × 一個基因（只含進入評分的 12,629 個人類 Gene ID）；同一配對取第一次出現的上下文
抽樣：簡單隨機抽 100 個配對，固定種子 20260911
輸入：pubtator_all_tsv_results/20250912_Results_gene_descriptions.tsv
      circadian_genes_results_dual_optimized.csv（Gene_ID、Gene_Symbol、Confidence_Level）
輸出：pubtator_spotcheck_sample_20260911.xlsx（抽查表＋說明）、pubtator_spotcheck_sample_20260911.tsv
執行位置：/home/yubau/aging/20250519_circadian
v2（9/11）：改用不處理引號的模式讀取（v1 只讀到 2,530,024 筆，SQL 表為 2,570,640 筆）；筆數不符即停止
"""
import csv
import pandas as pd

DESC = "pubtator_all_tsv_results/20250912_Results_gene_descriptions.tsv"
SCORED = "circadian_genes_results_dual_optimized.csv"
N_SAMPLE, SEED = 100, 20260911
OUT = "pubtator_spotcheck_sample_20260911_v2"
EXPECTED_ROWS = 2570640

sc = pd.read_csv(SCORED, encoding="utf-8-sig", escapechar="\\", dtype=str, keep_default_na=False,
                 usecols=["Gene_ID", "Gene_Symbol", "Confidence_Level"])
if len(sc) != 12629:
    raise SystemExit(f"評分表列數不是 12,629，而是 {len(sc):,}，先停止")
human = set(sc["Gene_ID"])

d = pd.read_csv(DESC, sep="\t", dtype=str, keep_default_na=False, quoting=csv.QUOTE_NONE, low_memory=False,
                usecols=["pmid", "gene_id", "gene_name", "gene_text", "species_name", "context", "passage_type"])
print(f"gene_descriptions 全部 {len(d):,} 筆")
if len(d) != EXPECTED_ROWS:
    raise SystemExit(f"讀到 {len(d):,} 筆，與 SQL 表 Results_gene_descriptions_20250912 的 {EXPECTED_ROWS:,} 筆不符，先停止")
d = d[d["gene_id"].isin(human)]
print(f"屬於 12,629 個人類評分基因的提及 {len(d):,} 筆")
pairs = d.drop_duplicates(["pmid", "gene_id"], keep="first")
print(f"不重複的「文章 × 基因」配對 {len(pairs):,} 個，涵蓋 {pairs['gene_id'].nunique():,} 個基因、{pairs['pmid'].nunique():,} 篇文章")
if len(pairs) < N_SAMPLE:
    raise SystemExit("配對數少於抽樣數，先停止")

s = pairs.sample(n=N_SAMPLE, random_state=SEED).merge(sc, left_on="gene_id", right_on="Gene_ID", how="left")
s = s.reset_index(drop=True)
s.insert(0, "sample_no", range(1, len(s) + 1))
s["pubmed_url"] = "https://pubmed.ncbi.nlm.nih.gov/" + s["pmid"] + "/"
s["correct_gene（Y／N／?）"] = ""
s["circadian_context（Y／N／?）"] = ""
s["note"] = ""
cols = ["sample_no", "pmid", "pubmed_url", "gene_id", "Gene_Symbol", "Confidence_Level", "gene_text", "species_name",
        "passage_type", "context", "correct_gene（Y／N／?）", "circadian_context（Y／N／?）", "note"]
s = s[cols]
s.to_csv(f"{OUT}.tsv", sep="\t", index=False)

print("\n樣本中各文獻等級的配對數：")
print(s["Confidence_Level"].value_counts().to_string())
print("\n母體（全部配對）中各文獻等級的比例：")
print((pairs.merge(sc, left_on="gene_id", right_on="Gene_ID")["Confidence_Level"].value_counts(normalize=True) * 100).round(1).to_string())

guide = pd.DataFrame({"項目": [
    "目的", "抽樣", "correct_gene", "circadian_context", "判讀方式", "回報"],
    "說明": [
    "估計 PubTator3 把文中提及對應到人類基因的正確率，以及這些配對是否真的在講生理時鐘（審稿 R2 缺點 5、建議 3）",
    f"從 {len(pairs):,} 個「文章 × 基因」配對中，以固定種子 {SEED} 簡單隨機抽 {N_SAMPLE} 個；每個配對取第一次出現的上下文（前後各 50 個字元，** 標出提及）",
    "Y：上下文中的提及確實指這個人類基因（Gene_Symbol）；N：指的是別的東西（例如縮寫撞名、別的物種、細胞類型、試劑、人名）；?：上下文不足以判斷，可點 pubmed_url 看摘要",
    "Y：這段上下文或該篇摘要在討論生理時鐘、晝夜節律或相關現象；N：只是同一篇文章裡剛好提到；?：無法判斷",
    "每列填 Y、N 或 ?；有疑問的在 note 簡單記原因，例如「縮寫撞名」「小鼠基因」",
    "填完傳回，會計算兩項的正確率與 95% 信賴區間（Wilson），寫進回覆信與補充資料"]})
with pd.ExcelWriter(f"{OUT}.xlsx", engine="openpyxl") as w:
    s.to_excel(w, sheet_name="抽查表", index=False)
    guide.to_excel(w, sheet_name="說明", index=False)
    ws = w.sheets["抽查表"]
    widths = {"A": 9, "B": 11, "C": 38, "D": 11, "E": 13, "F": 12, "G": 16, "H": 14, "I": 10, "J": 90, "K": 16, "L": 20, "M": 24}
    for col, wd in widths.items():
        ws.column_dimensions[col].width = wd
    ws.freeze_panes = "B2"
    for row in range(2, len(s) + 2):
        c = ws[f"C{row}"]; c.hyperlink = c.value; c.style = "Hyperlink"
    g = w.sheets["說明"]; g.column_dimensions["A"].width = 18; g.column_dimensions["B"].width = 120
print(f"\n已輸出 {OUT}.xlsx（抽查表＋說明）與 {OUT}.tsv")