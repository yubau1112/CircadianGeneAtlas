#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
論文晝夜節律關聯性分析器 - 雙重分析版本

基於基因分析版本改為論文雙重分析：
1. 保持原有的JSON快取格式兼容性
2. 保持原有的評分算法不變
3. 雙重分析：標題分析 + 摘要分析
4. 分別計算標題和摘要的晝夜節律相關性評分
5. 輸出詳細的雙重評分結果和置信度
"""

import os
import sys
import json
import time
import random
import argparse
import traceback
import csv
import glob
import logging
import pandas as pd
from datetime import datetime
from collections import defaultdict
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import lru_cache
import threading

# 設置基本目錄
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_PAPER_CACHE = os.path.join(SCRIPT_DIR, "paper_analysis_cache.json")
DEFAULT_OUTPUT_TSV = os.path.join(SCRIPT_DIR, "paper_circadian_results.tsv")

# 確保日誌目錄存在
LOG_DIR = os.path.join(SCRIPT_DIR, "logs")
os.makedirs(LOG_DIR, exist_ok=True)

# 設置日誌
LOG_FILE = os.path.join(LOG_DIR, f"paper_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

# 全局編譯模式實例

cache_lock = threading.Lock()

# 已知的核心晝夜節律基因列表
VERIFIED_CIRCADIAN_GENES = [
    "CLOCK", "ARNTL", "BMAL1", "ARNTL2", "BMAL2",
    "PER1", "PER2", "PER3", 
    "CRY1", "CRY2",
    "NR1D1", "REV-ERBA", "NR1D2", "REV-ERBB",
    "RORA", "RORB", "RORC",
    "NPAS2", "CSNK1D", "CSNK1E",
    "TIMELESS", "DEC1", "DEC2", "BHLHE40", "BHLHE41"
]

# 轉換為set以提高查找速度
VERIFIED_CIRCADIAN_GENES_SET = set(VERIFIED_CIRCADIAN_GENES)

# 晝夜節律關鍵字分類
CORE_CIRCADIAN_KEYWORDS = [
    # 核心術語 - 使用詞邊界匹配
    r"\bcircadian\b", r"\bclock\b", r"\brhythm\b", r"\bdiurnal\b", r"\bnocturnal\b",
    
    # 擴展術語 - 精確片語匹配
    "circadian rhythm", "biological clock", "clock genes", "circadian gene", 
    "circadian regulation", "circadian oscillation", "circadian oscillator",
    "circadian timing", "clock function", "clock protein", "clock-controlled",
    "day/night cycle", "light/dark cycle", "light cycle", "dark cycle"
]

PATHWAY_CIRCADIAN_KEYWORDS = [
    # 核心生理時鐘基因
    r"\bCLOCK\b(?:\s+gene|\s+protein|$|\s)", 
    r"\bBMAL1\b", r"\bARNTL\b",
    r"\bPER1\b", r"\bPER2\b", r"\bPER3\b",
    r"\bCRY1\b", r"\bCRY2\b",
    r"\bREV-ERB\b", r"\bROR\b",
    "period gene", r"\bcryptochrome\b", r"\btimeless\b", 
    "suprachiasmatic", "suprachiasmatic nucleus",
    
    # 生理時鐘調節相關
    "sleep/wake", "sleep-wake", r"\bmelatonin\b", "SCN",
    "phase shift", "phase response", r"\bentrainment\b", "free-running"
]

MECHANISM_CIRCADIAN_KEYWORDS = [
    r"\bzeitgeber\b", r"\bphotoperiod\b", r"\bchronobiology\b", r"\bchronobiological\b",
    
    # 時間相關術語
    "daily variation", "daily oscillation", "daily rhythm", "daily cycle",
    "24-hour", "24 hour", "24h", "24 h"
]

DISEASE_CIRCADIAN_KEYWORDS = [
    "circadian disorder", "sleep disorder", "jet lag", "shift work disorder"
]

RESEARCH_CIRCADIAN_KEYWORDS = [
    "chronotherapy", "circadian medicine", "sleep research"
]

# 需要上下文驗證的關鍵字
CONTEXT_DEPENDENT_KEYWORDS = [
    r"\btime\b", r"\bdaily\b", r"\bnight\b", r"\bday\b"
]

# === 優化2: 加快關鍵字查找 ===
    
def find_circadian_keywords_in_text(text, keywords_list):
    """在文本中查找晝夜節律關鍵字"""
    if not text:
        return []
    
    found_keywords = []
    text_lower = text.lower()
    
    for keyword in keywords_list:
        try:
            if keyword.startswith(r'\b') or '\\b' in keyword:
                # 正則表達式匹配
                if re.search(keyword, text_lower, re.IGNORECASE):
                    match = re.search(keyword, text_lower, re.IGNORECASE)
                    if match:
                        found_keywords.append(match.group())
            else:
                # 簡單字串匹配
                if keyword.lower() in text_lower:
                    found_keywords.append(keyword)
        except re.error:
            # 正則表達式錯誤時回退到簡單匹配
            if keyword.lower() in text_lower:
                found_keywords.append(keyword)
    
    return found_keywords

def analyze_text_content(text, text_type):
    """分析單一文本內容（標題或摘要）- 完全按照原始邏輯 + 返回關鍵字"""
    score = 0
    detailed_reasons = []
    found_keywords = []  # 新增：收集找到的關鍵字
    
    if not text or not text.strip():
        return 0, [f"無{text_type}內容"], []
    
    text_lower = text.lower()
    
    # 1. 已知晝夜節律基因加分 (20分)
    circadian_gene_mentioned = False
    for gene in VERIFIED_CIRCADIAN_GENES:
        if re.search(rf'\b{re.escape(gene)}\b', text, re.IGNORECASE):
            circadian_gene_mentioned = True
            found_keywords.append(gene)  # 添加基因到關鍵字列表
            break

    if circadian_gene_mentioned:
        score += 20
        detailed_reasons.append(f"{text_type}提及已知晝夜節律基因 (+20)")

    # 2. 核心關鍵字匹配 (每個3分，最多24分)
    core_matches = 0
    for keyword in CORE_CIRCADIAN_KEYWORDS:
        if re.search(keyword, text_lower, re.IGNORECASE):
            core_matches += 1
            # 收集實際匹配的文本
            match = re.search(keyword, text_lower, re.IGNORECASE)
            if match:
                found_keywords.append(match.group())
            if core_matches <= 8:
                score += 3
                detailed_reasons.append(f"{text_type}核心關鍵字: {keyword} (+3)")
    
    # 3. 通路關鍵字匹配 (每個2分，最多16分)
    pathway_matches = 0
    for keyword in PATHWAY_CIRCADIAN_KEYWORDS:
        if re.search(keyword, text_lower, re.IGNORECASE):
            pathway_matches += 1
            # 收集實際匹配的文本
            match = re.search(keyword, text_lower, re.IGNORECASE)
            if match:
                found_keywords.append(match.group())
            if pathway_matches <= 8:
                score += 2
                detailed_reasons.append(f"{text_type}通路關鍵字: {keyword} (+2)")
    
    # 4. 機制關鍵字匹配 (每個1.5分，最多12分)
    mechanism_matches = 0
    for keyword in MECHANISM_CIRCADIAN_KEYWORDS:
        if re.search(keyword, text_lower, re.IGNORECASE):
            mechanism_matches += 1
            # 收集實際匹配的文本
            match = re.search(keyword, text_lower, re.IGNORECASE)
            if match:
                found_keywords.append(match.group())
            if mechanism_matches <= 8:
                score += 1.5
                detailed_reasons.append(f"{text_type}機制關鍵字: {keyword} (+1.5)")
    
    # 5. 疾病關鍵字匹配 (每個1分，最多8分)
    disease_matches = 0
    for keyword in DISEASE_CIRCADIAN_KEYWORDS:
        if re.search(keyword, text_lower, re.IGNORECASE):
            disease_matches += 1
            # 收集實際匹配的文本
            match = re.search(keyword, text_lower, re.IGNORECASE)
            if match:
                found_keywords.append(match.group())
            if disease_matches <= 8:
                score += 1
                detailed_reasons.append(f"{text_type}疾病關鍵字: {keyword} (+1)")
    
    # 6. 上下文驗證獎勵 (最多10分)
    context_bonus = 0
    for context_kw in CONTEXT_DEPENDENT_KEYWORDS:
        if re.search(context_kw, text_lower, re.IGNORECASE):
            # 檢查是否有晝夜節律相關上下文
            circadian_context_patterns = [
                r"circadian", r"clock", r"rhythm", r"diurnal",
                r"24-hour", r"daily", r"sleep-wake"
            ]
            has_circadian_context = any(
                re.search(pattern, text_lower, re.IGNORECASE) 
                for pattern in circadian_context_patterns
            )
            if has_circadian_context:
                context_bonus += 2
                # 收集上下文關鍵字
                match = re.search(context_kw, text_lower, re.IGNORECASE)
                if match:
                    found_keywords.append(match.group())
                detailed_reasons.append(f"{text_type}時間上下文: {context_kw} (+2)")  # ✅ 改為"時間"
                if context_bonus >= 10:
                    break
    
    score += context_bonus
    
    # 7. 多樣性獎勵 (5分) - 對單一文本來說就是有晝夜節律關鍵字就給
    if score > 10:  # 只有當已經有一定分數時才給多樣性獎勵
        has_circadian_keywords = any(re.search(kw, text, re.IGNORECASE) for kw in CORE_CIRCADIAN_KEYWORDS[:5])
        if has_circadian_keywords:
            score += 5
            detailed_reasons.append(f"{text_type}內容豐富獎勵 (+5)")
    
    # 限制總分不超過100
    score = min(100, score)
    
    # 去重關鍵字並保持順序
    unique_keywords = []
    seen = set()
    for kw in found_keywords:
        if kw.lower() not in seen:
            unique_keywords.append(kw)
            seen.add(kw.lower())
    
    return round(score, 1), detailed_reasons, unique_keywords

def process_paper_circadian_analysis_dual(pmid, title, abstract, paper_cache=None):
    """雙重分析：標題分析 + 摘要分析 - 修改版包含關鍵字"""

    # 1. 標題分析
    title_score, title_reasons, title_keywords = analyze_text_content(title, "標題")

    # 2. 摘要分析  
    abstract_score, abstract_reasons, abstract_keywords = analyze_text_content(abstract, "摘要")

    # 3. 綜合評分
    combined_score = title_score + abstract_score
   
    # 4. 判斷置信度和晝夜節律相關性（修改版）
    # 新邏輯：只要標題或摘要其中一個 >0分 就算晝夜節律相關
    title_has_circadian = title_score > 0
    abstract_has_circadian = abstract_score > 0
    is_circadian_related = title_has_circadian or abstract_has_circadian

    # 置信度還是基於綜合評分
    if combined_score >= 80:
        confidence = "HIGH"
    elif combined_score >= 40:
        confidence = "MEDIUM"
    elif combined_score >= 16:
        confidence = "LOW"
    else:
        confidence = "VERY_LOW"

    # 如果標題或摘要有晝夜節律相關但綜合評分很低，至少給LOW
    if is_circadian_related and confidence == "VERY_LOW":
        confidence = "LOW"
   
    # 5. 收集找到的關鍵字 - 合併標題和摘要的關鍵字（去重）
    all_keywords = title_keywords + abstract_keywords
    unique_keywords = []
    seen = set()
    for kw in all_keywords:
        if kw.lower() not in seen:
            unique_keywords.append(kw)
            seen.add(kw.lower())

    return {
       "pmid": pmid,
       "title": title if title else "",
       "abstract": abstract if abstract else "",
       "title_score": title_score,
       "abstract_score": abstract_score,
       "combined_score": combined_score,
       "circadian_related": is_circadian_related,
       "confidence_level": confidence,
       "title_reasons": title_reasons,
       "abstract_reasons": abstract_reasons,
       "title_keywords": title_keywords,      # 新增：標題關鍵字
       "abstract_keywords": abstract_keywords,  # 新增：摘要關鍵字
       "circadian_keywords": unique_keywords,       # 保持原有：合併關鍵字
       "analysis_version": "paper_dual_v1",
       "last_updated": datetime.now().isoformat()
    }

# === 保持原始快取格式不變 ===
def load_cache(cache_file):
    """保持原版本的快取載入邏輯"""
    if os.path.exists(cache_file):
        try:
            with open(cache_file, 'r', encoding='utf-8') as f:
                return json.load(f)
        except Exception as e:
            logger.error(f"讀取快取檔案時出錯: {e}")
            return {}
    return {}

def save_cache(cache, cache_file):
    """改進的快取保存邏輯，防止檔案損壞"""
    with cache_lock:  # 添加線程鎖
        try:
            # 先保存到臨時檔案
            temp_file = cache_file + '.tmp'
            with open(temp_file, 'w', encoding='utf-8') as f:
                json.dump(cache, f, ensure_ascii=False, indent=2)
            
            # 檢查臨時檔案是否正常
            if os.path.getsize(temp_file) > 0:
                # 原子性替換
                if os.name == 'nt':  # Windows
                    if os.path.exists(cache_file):
                        os.remove(cache_file)
                os.rename(temp_file, cache_file)
                logger.debug(f"快取已保存: {len(cache)} 條目")
            else:
                logger.error("臨時檔案為空，取消保存")
                os.remove(temp_file)
                
        except Exception as e:
            logger.error(f"保存快取檔案時出錯: {e}")
            # 清理臨時檔案
            temp_file = cache_file + '.tmp'
            if os.path.exists(temp_file):
                try:
                    os.remove(temp_file)
                except:
                    pass

# === 優化5: 智能批量處理 ===
def process_paper_batch_optimized(papers_df, paper_cache, cache_file, 
                                save_interval=100, max_workers=3):
    """優化版批量處理論文，支援並行處理"""
    total = len(papers_df)
    circadian_related_count = 0
    
    logger.info(f"開始優化處理 {total:,} 篇論文")
    
    # 處理單篇論文的函數
    def process_single_paper(paper_data):
        pmid, title, abstract, row_index = paper_data
        
        # 檢查快取
        cache_key = f"{pmid}_paper_analysis"
        if cache_key in paper_cache:
            cached_result = paper_cache[cache_key]
            if cached_result.get('analysis_version') == 'paper_dual_v1':
                return cached_result, cache_key
        
        # 進行分析
        paper_result = process_paper_circadian_analysis_dual(pmid, title, abstract)
        
        return paper_result, cache_key
    
    # 準備論文數據 - 分離快取和需要處理的論文
    paper_data_list = []
    papers_need_processing = []
    cached_results = []

    for i, (index, row) in enumerate(papers_df.iterrows()):
        pmid = str(row['pmid'])
        cache_key = f"{pmid}_paper_analysis"
        
        if cache_key in paper_cache and paper_cache[cache_key].get('analysis_version') == 'paper_dual_v1':
            # 已有快取，直接使用
            cached_results.append(paper_cache[cache_key])
        else:
            # 需要重新處理
            paper_data_list.append((
                pmid,
                str(row['title']) if pd.notna(row['title']) else "",
                str(row['abstract']) if pd.notna(row['abstract']) else "",
                i
            ))
            papers_need_processing.append(pmid)

    logger.info(f"從快取載入: {len(cached_results):,} 篇論文")
    logger.info(f"需要重新處理: {len(paper_data_list):,} 篇論文")
    if papers_need_processing:
        logger.info(f"需要處理的PMID前10個: {papers_need_processing[:10]}")
    
    # 統計快取中的晝夜節律相關論文
    cached_circadian_count = len([r for r in cached_results if r.get('circadian_related', False)])
    circadian_related_count += cached_circadian_count
    
    # 初始化新處理的結果列表
    new_results = []
    
    # 如果有需要處理的論文，使用線程池並行處理
    if paper_data_list:
        logger.info(f"開始並行處理，使用 {max_workers} 個線程")
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_paper = {executor.submit(process_single_paper, paper_data): paper_data for paper_data in paper_data_list}
            
            for i, future in enumerate(as_completed(future_to_paper)):
                try:
                    paper_result, cache_key = future.result()
                    
                    # 驗證結果格式
                    if isinstance(paper_result, dict) and 'pmid' in paper_result and 'analysis_version' in paper_result:
                        new_results.append(paper_result)
                    else:
                        logger.error(f"跳過無效並行結果: {paper_result}")
                        continue
                    
                    # 更新快取
                    paper_cache[cache_key] = paper_result
                    
                    # 統計新處理的晝夜節律相關論文
                    if paper_result.get('circadian_related', False):
                        circadian_related_count += 1
                    
                    # 計算總進度（包含快取的論文）
                    total_processed = len(cached_results) + (i + 1)
                    progress = total_processed / total * 100
                    
                    # 進度顯示
                    if (i + 1) % 1000 == 0 or (i + 1) == len(paper_data_list):
                        circadian_percentage = (circadian_related_count / total_processed) * 100
                        logger.info(f"進度: {total_processed:,}/{total:,} ({progress:.1f}%), 晝夜節律相關: {circadian_related_count:,} ({circadian_percentage:.1f}%)")
                    
                    # 定期保存快取
                    if (i + 1) % save_interval == 0 or (i + 1) == len(paper_data_list):
                        try:
                            save_cache(paper_cache, cache_file)
                            # 只在較大間隔時才顯示保存日誌
                            if (i + 1) % (save_interval * 10) == 0 or (i + 1) == len(paper_data_list):
                                logger.info(f"已保存快取 (新處理: {i+1:,}/{len(paper_data_list):,})")
                        except Exception as e:
                            logger.error(f"定期保存快取失敗: {e}")
                        
                except Exception as e:
                    logger.error(f"處理論文時出錯: {e}")
                    continue
    else:
        logger.info("所有論文都已在快取中，無需重新處理")
    
    # 合併快取結果和新處理結果
    all_results = cached_results + new_results
    
    # 最終統計
    final_circadian_count = len([r for r in all_results if r.get('circadian_related', False)])
    
    logger.info(f"處理完成！")
    logger.info(f"   - 總論文數: {len(all_results):,}")
    logger.info(f"   - 從快取載入: {len(cached_results):,} 篇")
    logger.info(f"   - 新處理: {len(new_results):,} 篇")
    logger.info(f"   - 晝夜節律相關: {final_circadian_count:,} ({final_circadian_count/len(all_results)*100:.1f}%)")
    
    # 驗證結果完整性
    if len(all_results) != total:
        logger.warning(f"結果數量不匹配: 預期 {total:,}, 實際 {len(all_results):,}")
    
    # 檢查是否有重複的PMID
    pmids = [r.get('pmid') for r in all_results]
    unique_pmids = len(set(pmids))
    if unique_pmids != len(all_results):
        logger.warning(f"發現重複PMID: {len(all_results) - unique_pmids} 個重複")
    
    return all_results

def main():
    # 設置命令行參數 - 改為論文分析參數
    parser = argparse.ArgumentParser(description='論文晝夜節律關聯性分析 - 優化版本')
    parser.add_argument('--paper-cache-file', type=str, 
                      default=os.path.join(SCRIPT_DIR, "papers_cache.json"),
                      help='論文快取檔案路徑（PubMed數據，包含標題和摘要）')
    parser.add_argument('--analysis-cache-file', type=str, default='paper_analysis_cache.json',
                      help='分析結果快取檔案路徑')
    parser.add_argument('--output-tsv', type=str, default='paper_circadian_results.tsv',
                      help='輸出TSV檔案路徑')
    parser.add_argument('--pmid-file', type=str, 
                      help='PMID列表文件路徑（可選，不指定則處理快取中所有論文）')
    parser.add_argument('--max-papers', type=int, default=0,
                      help='處理的最大論文數量(0表示不限制)')
    parser.add_argument('--start-index', type=int, default=0,
                      help='從論文列表中的指定索引開始處理')
    parser.add_argument('--save-interval', type=int, default=100,
                      help='每處理多少篇論文保存一次快取')
    parser.add_argument('--max-workers', type=int, default=4,
                      help='並行處理線程數')
    parser.add_argument('--min-circadian-score', type=float, default=8,
                      help='最低晝夜節律評分要求')
    args = parser.parse_args()
    
    # 載入分析結果快取
    analysis_cache = load_cache(args.analysis_cache_file)
    logger.info(f"已載入分析快取，包含{len(analysis_cache):,}篇論文")
    
    # 兼容性檢查
    if analysis_cache:
        sample_key = next(iter(analysis_cache.keys()))
        sample_version = analysis_cache[sample_key].get('analysis_version', '未知')
        logger.info(f"快取版本: {sample_version}")
        if sample_version == 'paper_dual_v1':
            logger.info("快取格式兼容，可直接使用")
        else:
            logger.warning("快取版本較舊，建議重新分析以獲得最佳效果")
    
    # 讀取目標PMID列表（如果指定）
    target_pmids = None
    if args.pmid_file:
        if not os.path.exists(args.pmid_file):
            logger.error(f"PMID文件不存在: {args.pmid_file}")
            sys.exit(1)
        
        with open(args.pmid_file, 'r') as f:
            target_pmids = [line.strip() for line in f if line.strip()]
        logger.info(f"讀取目標PMID: {len(target_pmids):,} 個")
    
    # 載入論文數據
    logger.info("載入論文數據...")
    start_time = time.time()
    papers_df = load_paper_data_from_cache(args.paper_cache_file, target_pmids)
    load_time = time.time() - start_time
    
    if papers_df.empty:
        logger.error("無法載入論文資料")
        sys.exit(1)
    
    logger.info(f"數據載入完成 (耗時 {load_time:.1f}s)")
    logger.info(f"載入了 {len(papers_df):,} 篇論文的資料")
    
    # 設置處理範圍
    if args.max_papers > 0 and args.max_papers < len(papers_df):
        logger.info(f"將僅處理前 {args.max_papers:,} 篇論文")
        papers_df = papers_df.iloc[:args.max_papers]
    
    if args.start_index > 0:
        if args.start_index < len(papers_df):
            papers_df = papers_df.iloc[args.start_index:]
            logger.info(f"從索引 {args.start_index} 處開始處理")
        else:
            logger.warning(f"指定的起始索引 {args.start_index} 超出了論文列表範圍")
    
    logger.info(f"即將處理 {len(papers_df):,} 篇論文")
    logger.info(f"並行線程數: {args.max_workers}")
    
    # 主處理流程
    if not papers_df.empty:
        main_start_time = time.time()
        
        results = process_paper_batch_optimized(
            papers_df, 
            analysis_cache, 
            args.analysis_cache_file,
            save_interval=args.save_interval,
            max_workers=args.max_workers
        )
        
        main_processing_time = time.time() - main_start_time
        
        # 根據最低評分要求篩選結果
        original_circadian_count = len([r for r in results if r.get('circadian_related', False)])
        
        if args.min_circadian_score > 0:
            for result in results:
                if result.get('combined_score', 0) < args.min_circadian_score:
                    result['circadian_related'] = False
            
            filtered_circadian_count = len([r for r in results if r.get('circadian_related', False)])
            logger.info(f"應用最低評分篩選({args.min_circadian_score}): {original_circadian_count:,} -> {filtered_circadian_count:,}")
        
        # 保存結果到TSV - 修正為正確的函數名
        tsv_file, csv_file = save_results_to_tsv_paper(results, args.output_tsv)
        
        # 在這裡添加文件路徑信息
        if tsv_file and csv_file:
            logger.info(f"輸出文件:")
            logger.info(f"   - TSV格式: {tsv_file}")
            logger.info(f"   - Excel專用CSV: {csv_file}")
        else:
            logger.error("文件保存失敗")        
        
        # 性能統計
        total_time = load_time + main_processing_time
        papers_per_second = len(results) / main_processing_time if main_processing_time > 0 else 0

        # 詳細統計
        unique_pmids = len(set(r['pmid'] for r in results))
        expected_papers = len(papers_df)

        logger.info(f"性能統計:")
        logger.info(f"   - 數據載入: {load_time:.1f}s")
        logger.info(f"   - 主處理: {main_processing_time:.1f}s")
        logger.info(f"   - 總耗時: {total_time:.1f}s")
        logger.info(f"   - 處理速度: {papers_per_second:.1f} 論文/秒")

        logger.info(f"數據完整性檢查:")
        logger.info(f"   - 預期論文數: {expected_papers:,}")
        logger.info(f"   - 實際處理: {len(results):,}")
        logger.info(f"   - 唯一PMID數: {unique_pmids:,}")

        if len(results) != expected_papers:
            logger.warning(f"   - 處理數量不符: 差異 {abs(len(results) - expected_papers):,} 篇")

        if unique_pmids != len(results):
            logger.warning(f"   - 發現重複處理: {len(results) - unique_pmids:,} 篇重複")
        
        # 生成並顯示摘要報告
        summary_report = generate_paper_summary_report(results)
        print(summary_report)
        logger.info("摘要報告已生成")
        
    else:
        logger.info("沒有需要處理的論文")
    
    # 保存最終快取
    save_cache(analysis_cache, args.analysis_cache_file)
    logger.info(f"已將所有 {len(analysis_cache):,} 個記錄保存到快取 {args.analysis_cache_file}")
    
    logger.info("處理完成!")

def load_paper_data_from_cache(cache_file, target_pmids=None):
   """從快取中載入論文數據"""
   if not os.path.exists(cache_file):
       logger.error(f"論文快取檔案不存在: {cache_file}")
       return pd.DataFrame()
   
   paper_cache = load_cache(cache_file)
   logger.info(f"從快取載入了 {len(paper_cache):,} 篇論文")
   
   # 在篩選之前，先分析所有失敗記錄
   deleted_pmids = []
   book_pmids = []
   other_failed_pmids = []

   for pmid, data in paper_cache.items():
       if not data.get('success', True):
           error = data.get('error', '')
           if 'deleted' in error.lower() or 'redirected' in error.lower():
               deleted_pmids.append(pmid)
           elif 'book' in error.lower() or 'report' in error.lower():
               book_pmids.append(pmid)
           else:
               other_failed_pmids.append(pmid)

   # 顯示排除信息
   if deleted_pmids:
       logger.info(f"排除已刪除記錄: {len(deleted_pmids)} 篇 (PMID已被刪除/重定向) - PMID: {deleted_pmids}")
   if book_pmids:
       logger.info(f"排除書籍/報告記錄: {len(book_pmids)} 篇 (非期刊論文類型) - PMID: {book_pmids}")
   if other_failed_pmids:
       logger.info(f"排除其他失敗記錄: {len(other_failed_pmids)} 篇 (其他錯誤) - PMID: {other_failed_pmids}")
   
   # 然後進行原有的篩選邏輯
   if target_pmids:
       target_pmids_set = set(map(str, target_pmids))
       filtered_cache = {pmid: data for pmid, data in paper_cache.items() 
                        if pmid in target_pmids_set and data.get('success', True)}
       logger.info(f"根據PMID列表篩選後: {len(filtered_cache):,} 篇論文")
   else:
       filtered_cache = {pmid: data for pmid, data in paper_cache.items() 
                        if data.get('success', True)}
       logger.info(f"篩選成功記錄後: {len(filtered_cache):,} 篇論文")
   
   # 轉換為 DataFrame 格式
   papers_data = []
   for pmid, data in filtered_cache.items():
       papers_data.append({
           'pmid': pmid,
           'title': data.get('title', ''),
           'abstract': data.get('abstract', '')
       })
   
   papers_df = pd.DataFrame(papers_data)
   
   # 標準的內容檢查（標題和摘要都為空但success=True的情況）
   initial_count = len(papers_df)
   empty_content_mask = (
       (papers_df['title'].str.strip() == '') & 
       (papers_df['abstract'].str.strip() == '')
   )
   empty_papers = papers_df[empty_content_mask]
   empty_pmids = empty_papers['pmid'].astype(str).tolist()

   # 如果有success=True但內容為空的記錄，也顯示出來
   if empty_pmids:
       logger.info(f"排除無內容記錄: {len(empty_pmids)} 篇 (success=True但無標題或摘要) - PMID: {empty_pmids}")

   # 過濾掉無內容論文
   papers_df = papers_df[~empty_content_mask]
   
   return papers_df

def save_results_to_tsv_paper(results, tsv_file):
    """保存論文分析結果到TSV文件 - Excel兼容版本 + 關鍵字欄位"""
    try:
        # 生成帶時間戳的檔名
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        name, ext = os.path.splitext(tsv_file)
        if ext.lower() != '.tsv':
            ext = '.tsv'
        timestamped_file = f"{name}_{timestamp}{ext}"
        
        # 同時生成CSV版本給Excel使用
        csv_file = f"{name}_{timestamp}.csv"
        
        # Excel友好的清理函數
        def clean_for_excel(text):
            if isinstance(text, str):
                # 移除所有控制字符和可能導致問題的字符
                text = text.replace('\t', ' ')      # TAB → 空格
                text = text.replace('\n', ' ')      # 換行 → 空格  
                text = text.replace('\r', ' ')      # 回車 → 空格
                text = text.replace('\0', '')       # NULL字符
                text = text.replace('"', '""')      # 雙引號轉義
                text = text.replace('\\', '')       # 反斜線
                # 規範化空白
                text = ' '.join(text.split())
                
                # 如果包含逗號、換行或雙引號，用雙引號包圍（CSV格式）
                if ',' in text or '\n' in text or '"' in text:
                    text = f'"{text}"'
                    
                return text
            return str(text) if text is not None else ''
        
        # 新增：處理關鍵字列表的函數
        def format_keywords(keywords_list):
            if isinstance(keywords_list, list):
                # 用分號分隔關鍵字
                return ';'.join(keywords_list)
            return ''
        
        # 方法1: 保存標準TSV（UTF-8）
        with open(timestamped_file, 'w', encoding='utf-8', newline='') as f:
            # 寫入標頭 - 新增關鍵字欄位
            headers = [
                'PMID', 'Title', 'Abstract', 'Title_Score', 'Abstract_Score', 'Combined_Score',
                'Is_Circadian_Related', 'Confidence_Level', 'Title_Keywords', 'Abstract_Keywords',
                'Title_Reasons', 'Abstract_Reasons'
            ]
            f.write('\t'.join(headers) + '\n')
            
            # 寫入數據
            for result in results:
                if not isinstance(result, dict) or 'pmid' not in result:
                    logger.warning(f"跳過無效記錄: {type(result)}")
                    continue
                
                try:
                    row_data = [
                        clean_for_excel(result['pmid']),
                        clean_for_excel(result.get('title', '')),
                        clean_for_excel(result.get('abstract', '')),
                        str(result.get('title_score', 0)),
                        str(result.get('abstract_score', 0)),
                        str(result.get('combined_score', 0)),
                        'Yes' if result.get('circadian_related', False) else 'No',
                        clean_for_excel(result.get('confidence_level', 'VERY_LOW')),
                        clean_for_excel(format_keywords(result.get('title_keywords', []))),
                        clean_for_excel(format_keywords(result.get('abstract_keywords', []))),
                        clean_for_excel(' | '.join(result.get('title_reasons', []))),
                        clean_for_excel(' | '.join(result.get('abstract_reasons', [])))
                    ]
                    
                    f.write('\t'.join(row_data) + '\n')
                    
                except Exception as e:
                    logger.error(f"處理PMID {result.get('pmid', 'Unknown')} 時出錯: {e}")
                    continue
        
        # 方法2: 保存Excel專用CSV（UTF-8 with BOM）
        import csv
        with open(csv_file, 'w', encoding='utf-8-sig', newline='') as f:
            # UTF-8-sig 會自動添加BOM，讓Excel正確識別編碼
            writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            
            # 寫入標頭 - 新增關鍵字欄位
            headers = [
                'PMID', 'Title', 'Abstract', 'Title_Score', 'Abstract_Score', 'Combined_Score',
                'Is_Circadian_Related', 'Confidence_Level', 'Title_Keywords', 'Abstract_Keywords',
                'Title_Reasons', 'Abstract_Reasons'
            ]
            writer.writerow(headers)
            
            # 寫入數據
            for result in results:
                if not isinstance(result, dict) or 'pmid' not in result:
                    continue
                
                try:
                    # 簡化的清理函數用於CSV
                    def clean_for_csv(text):
                        if isinstance(text, str):
                            text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
                            text = ' '.join(text.split())
                            return text
                        return str(text) if text is not None else ''
                    
                    row_data = [
                        clean_for_csv(result['pmid']),
                        clean_for_csv(result.get('title', '')),
                        clean_for_csv(result.get('abstract', '')),
                        result.get('title_score', 0),
                        result.get('abstract_score', 0),
                        result.get('combined_score', 0),
                        'Yes' if result.get('circadian_related', False) else 'No',
                        clean_for_csv(result.get('confidence_level', 'VERY_LOW')),
                        format_keywords(result.get('title_keywords', [])),
                        format_keywords(result.get('abstract_keywords', [])),
                        clean_for_csv(' | '.join(result.get('title_reasons', []))),
                        clean_for_csv(' | '.join(result.get('abstract_reasons', [])))
                    ]
                    
                    writer.writerow(row_data)
                    
                except Exception as e:
                    logger.error(f"處理PMID {result.get('pmid', 'Unknown')} 時出錯: {e}")
                    continue
        
        logger.info(f"TSV結果已保存至 {timestamped_file}")
        logger.info(f"Excel專用CSV已保存至 {csv_file}")
        logger.info(f"成功寫入 {len(results):,} 篇論文記錄")
        logger.info(f"")
        logger.info(f"新增欄位:")
        logger.info(f"  - Title_Keywords: 標題中發現的晝夜節律關鍵字 (分號分隔)")
        logger.info(f"  - Abstract_Keywords: 摘要中發現的晝夜節律關鍵字 (分號分隔)")
        logger.info(f"")
        logger.info(f"使用建議:")
        logger.info(f"  - 用一般文本編輯器/數據分析工具: 使用 {timestamped_file}")
        logger.info(f"  - 用Excel打開: 使用 {csv_file} (含BOM，解決中文亂碼)")
        
        return timestamped_file, csv_file
        
    except Exception as e:
        logger.error(f"保存文件時出錯: {e}")
        import traceback
        traceback.print_exc()
        return None, None

def generate_paper_summary_report(results):
    """生成論文分析摘要報告"""
    total_papers = len(results)
    circadian_related = [r for r in results if r.get('circadian_related', False)]
    
    # 按置信度分類
    high_conf = [r for r in circadian_related if r.get('confidence_level') == 'HIGH']
    medium_conf = [r for r in circadian_related if r.get('confidence_level') == 'MEDIUM']
    low_conf = [r for r in circadian_related if r.get('confidence_level') == 'LOW']
    
    # 按評分排序
    circadian_related_sorted = sorted(circadian_related, key=lambda x: x.get('combined_score', 0), reverse=True)
    
    # 評分分布統計
    if circadian_related:
        scores = [r.get('combined_score', 0) for r in circadian_related]
        avg_score = sum(scores) / len(scores)
        max_score = max(scores)
        min_score = min(scores)
    else:
        avg_score = max_score = min_score = 0
    
    # 生成報告
    report = f"""
    
{'='*80}
論文晝夜節律關聯性分析報告
{'='*80}

總體統計:
- 總處理論文數: {total_papers:,}
- 晝夜節律相關論文數: {len(circadian_related):,} ({len(circadian_related)/total_papers*100:.1f}%)

置信度分布:
- 高置信度 (評分≥80): {len(high_conf):,} 篇 ({len(high_conf)/len(circadian_related)*100:.1f}% of circadian papers)
- 中等置信度 (評分40-79): {len(medium_conf):,} 篇 ({len(medium_conf)/len(circadian_related)*100:.1f}% of circadian papers)  
- 低置信度 (評分16-39): {len(low_conf):,} 篇 ({len(low_conf)/len(circadian_related)*100:.1f}% of circadian papers)

評分統計 (晝夜節律相關論文):
- 平均評分: {avg_score:.1f}
- 最高評分: {max_score:.1f}
- 最低評分: {min_score:.1f}

Top 15 晝夜節律相關論文 (按評分排序):
"""
    
    for i, paper in enumerate(circadian_related_sorted[:15], 1):
        pmid = paper.get('pmid', 'Unknown')
        title = paper.get('title', 'No title')[:60] + "..." if len(paper.get('title', '')) > 60 else paper.get('title', 'No title')
        score = paper.get('combined_score', 0)
        confidence = paper.get('confidence_level', 'Unknown')
        report += f"{i:2d}. PMID:{pmid} (評分:{score:4.1f}, {confidence:<6}) {title}\n"
    
    # 分析洞察
    report += f"\n分析洞察:\n"
    
    # 高分論文統計 - 調整門檻
    very_high_score_papers = [r for r in circadian_related if r.get('combined_score', 0) >= 90]
    high_score_papers = [r for r in circadian_related if 70 <= r.get('combined_score', 0) < 90]
    medium_score_papers = [r for r in circadian_related if 40 <= r.get('combined_score', 0) < 70]

    if very_high_score_papers:
        report += f"- 極高分論文 (≥90分): {len(very_high_score_papers):,} 篇，核心晝夜節律研究\n"
    if high_score_papers:
        report += f"- 高分論文 (70-89分): {len(high_score_papers):,} 篇，重要晝夜節律研究\n"
    if medium_score_papers:
        report += f"- 中等分數論文 (40-69分): {len(medium_score_papers):,} 篇，相關晝夜節律研究\n"
    
    # 關鍵字分析
    all_keywords = []
    for paper in circadian_related:
        all_keywords.extend(paper.get('circadian_keywords', []))
    
    if all_keywords:
        from collections import Counter
        keyword_counts = Counter(all_keywords)
        top_keywords = keyword_counts.most_common(10)
        report += f"- 最常見晝夜節律關鍵字: {', '.join([f'{kw}({count:,})' for kw, count in top_keywords[:5]])}\n"
    
    report += f"\n建議:\n"
    report += f"- 高置信度論文可作為晝夜節律研究的核心文獻\n"
    report += f"- 中等置信度論文可用於擴展研究範圍\n"
    report += f"- 低置信度論文可能包含間接相關的有用信息\n"
    
    report += f"\n{'='*80}\n"
    
    return report

if __name__ == "__main__":
    try:
        print("\n" + "="*80)
        print("論文晝夜節律關聯性分析腳本 - 雙重分析版本")
        print("特色: 標題和摘要分別分析PubMed論文的晝夜節律相關性")
        print("用法: 基於論文標題和摘要進行獨立的晝夜節律關鍵字分析")
        print("輸出: 雙重評分、綜合置信度和詳細分析報告")
        print("="*80 + "\n")
        
        main()
    except KeyboardInterrupt:
        print("\n程序被用戶中斷")
        logger.warning("程序被用戶中斷")
    except Exception as e:
        print(f"發生未預期的錯誤: {e}")
        logger.exception(f"發生未預期的錯誤: {e}")
        traceback.print_exc()
        sys.exit(1)