# Circadian Gene Atlas (CGA)

Code and processed data for the manuscript
**“A computational framework for literature-guided discovery of circadian-associated genes using pan-cancer multi-omics and machine learning”** (APBC 2026, submission 98).

Web platform: https://bal.asia.edu.tw/CircadianGeneAtlas/

## Overview

| Stage | Folder | Description |
|---|---|---|
| 1. Literature mining and scoring | `stage1_literature/` | PubTator3 gene extraction from 149,757 PubMed records, human-gene filtering, dual literature/annotation scoring (Section 2.4, Supplementary Table S5), article-level relevance scoring |
| 2. Machine learning | `stage2_ml/` | TCGA co-expression and promoter-methylation features, background correction, L1-regularized logistic regression, candidate prioritization, alternative strategies (Figure 5), enrichment and STRING network |
| Figures | `figures/` | One script per manuscript figure (see table below) |
| Revision analyses | `revision/` | Negative-set sensitivity, candidate stability and precision, benchmark recheck, re-ranking without the curated prior, PubTator3 spot-check sampling, background-correction check, time-series validation (GSE39445, GSE108539) |
| Processed data | `data/` | Feature matrices, final scores, candidate list, positive labels and reference gene sets |

## Figures and scripts

| Figure | Script |
|---|---|
| Figure 1 | `figures/fig1_overview.py` |
| Figure 2 | `figures/fig2_workflow.py` |
| Figures 3, 5, S1 | `figures/fig3_fig5_figS1_roc_strategy_importance.py` |
| Figure 4 | `figures/fig4_amplitude.py` |
| Figures 6, 7 | `figures/fig6_fig7_prioritization_heatmap.py` |
| Figure 8 | `figures/fig8_benchmark_data.py` (data) and `figures/fig8_benchmark_plot.py` (plot) |
| Figures S2, S6 | `figures/figS2_figS6_scores_threshold.py` |
| Figure S3 | `figures/figS3_background_correction.py` |
| Figures S4, S5 | `figures/figS4_figS5_pr_confusion.py` |
| Figure S7 | `figures/figS7_case_scatter.py` |
| Figures 9, 10 | `stage2_ml/evaluation/enrichment_network.py` |

## Data sources (not redistributed)

- TCGA gene expression, DNA methylation, copy number, somatic mutation, microRNA expression and clinical data: NCI Genomic Data Commons, Data Release 38.0 (primary tumor and adjacent normal samples, 24 cancer types); gene-level mutation frequencies in Table 3 from the GDC somatic mutation snapshot of 24 May 2024 (all programs)
- GTEx v8 expression (used only for the excluded leakage-control branch and Figure S7)
- PubMed records (searched 2 August 2025) and PubTator3 annotations (retrieved 4 August 2025)
- CircaDB human dataset (rhythmic genes inferred from GTEx by CYCLOPS; Ruben et al., Sci Transl Med 2018), accessed 5 July 2026; the query result is included as `data/CircaDB_HA_query.csv`
- Time-series transcriptomes: GEO GSE39445 (blood) and GSE108539 (skeletal muscle)
- GENCODE v36 annotation (Ensembl-to-symbol mapping)
- GENCODE v48 annotation and MyGene.info (gene annotation integration for the Stage 1 gene table; Ensembl identifiers retrieved 28 July 2025)

In the original analysis, TCGA data were read from a local SQL Server database. Scripts that query the database read the connection settings from environment variables (see `.env.example`); the processed feature matrices in `data/` allow the model, candidates and figures to be reproduced without database access.

## Environment

Python 3.11.9 with the packages in `requirements.txt`. Database access additionally requires `pyodbc` and Microsoft ODBC Driver 18 for SQL Server.

```bash
pip install -r requirements.txt
export CGA_SQL_SERVER=<host,port>      # only for scripts that query the database
export CGA_SQL_USER=<user>
export CGA_SQL_PASSWORD=<your password>
export NCBI_API_KEY=<your key>         # optional; faster NCBI E-utilities
export NCBI_EMAIL=you@example.org      # required by NCBI for E-utilities
```

Most scripts expect to be run from the folder that contains their input files; input and output file names are listed in the header of each script.

## Citation

Please cite the APBC 2026 paper (citation to be added upon publication).

## License

MIT License (see `LICENSE`).
