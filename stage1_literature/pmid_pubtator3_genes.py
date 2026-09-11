#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
高效能PMID批量處理程式
- 單線程設計，避免鎖競爭
- 記憶體內緩存，快速基因物種查詢
- 合併輸出檔案，減少檔案I/O
- 詳細日誌記錄處理過程
"""

import json
import requests
import sys
import datetime
import os
import time
import logging
from collections import defaultdict
import pandas as pd

# ===== NCBI API 設定 =====
NCBI_EMAIL = os.environ.get("NCBI_EMAIL", "")  # 您的郵箱
NCBI_API_KEY = os.environ.get("NCBI_API_KEY", "")  # 您的API Key

# 設定檔案路徑
GENE_CACHE_FILE = "gene_species_cache.json"
PMID_FILE = "circadian_all_pmid/circadian_complete_pmids_20250802_195005.tsv"
JSON_BASE_DIR = "circadian_all_pmid/json_files"
OUTPUT_DIR = "pubtator_all_tsv_results"

# 確保輸出目錄存在
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 設定日誌
def setup_logging():
    """設定日誌系統"""
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = "pmid_batch_processing_{}.log".format(timestamp)
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, encoding='utf-8'),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)

logger = setup_logging()

class GeneCache:
    """基因緩存管理類"""
    
    def __init__(self, cache_file=GENE_CACHE_FILE):
        self.cache_file = cache_file
        self.cache = {}
        self.new_entries = 0
        self.load_cache()     
    
    def load_cache(self):
        """載入基因緩存到記憶體"""
        if os.path.exists(self.cache_file):
            try:
                with open(self.cache_file, 'r', encoding='utf-8') as f:
                    self.cache = json.load(f)
                logger.info("載入基因緩存: {} 個基因".format(len(self.cache)))
            except Exception as e:
                logger.error("載入緩存失敗: {}".format(e))
                self.cache = {}
        else:
            logger.info("緩存檔案不存在，將創建新緩存")
            self.cache = {}
    
    def get_gene_species(self, gene_id):
        """從緩存獲取基因物種資訊"""
        gene_id = str(gene_id)
        if gene_id in self.cache:
            cached_info = self.cache[gene_id]
            return cached_info.get('species_name', 'unknown')
        return None
    
    def query_gene_from_api(self, gene_id):
        """從NCBI API查詢基因資訊"""
        gene_id = str(gene_id)
        
        try:
            # 添加 email 和 api_key 參數
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            params = {
                'db': 'gene',
                'id': gene_id,
                'retmode': 'json',
                'email': NCBI_EMAIL,
                'api_key': NCBI_API_KEY
            }
            
            response = requests.get(url, params=params, timeout=30)
            
            if response.status_code == 200:
                gene_info = response.json()
                if 'result' in gene_info and gene_id in gene_info['result']:
                    organism = gene_info['result'][gene_id].get('organism', {})
                    taxid = str(organism.get('taxid', ''))
                    species_name = organism.get('scientificname', 'unknown')
                    
                    # 加入緩存
                    cache_entry = {
                        'species': taxid,
                        'species_name': species_name,
                        'last_updated': datetime.datetime.now().isoformat()
                    }
                    
                    self.cache[gene_id] = cache_entry
                    self.new_entries += 1
                    
                    return species_name
            
            logger.warning("API查詢失敗: {}".format(gene_id))
            return 'unknown'
            
        except Exception as e:
            logger.error("API查詢異常 ({}): {}".format(gene_id, e))
            return 'unknown'
    
    def get_gene_species_name(self, gene_id):
        """獲取基因物種名稱（支援多基因ID拆分）"""
        original_gene_id = str(gene_id)
        
        # 檢查是否為空值或None
        if not gene_id or gene_id in ['None', 'none', '', 'null']:
            return 'unknown'
        
        # 檢查是否為多基因ID
        if ';' in original_gene_id:
            gene_ids = [gid.strip() for gid in original_gene_id.split(';') if gid.strip() and gid.strip() not in ['None', 'none', 'null']]
            
            if not gene_ids:
                return 'unknown'
                
            # logger.info("多基因ID檢測: {}，拆分為: {}".format(original_gene_id, gene_ids))
            
            # 分別查詢每個基因ID，返回結果列表
            results = []
            for single_gene_id in gene_ids:
                # 先檢查緩存
                species_name = self.get_gene_species(single_gene_id)
                if species_name is not None:
                    results.append((single_gene_id, species_name))
                else:
                    # 緩存中沒有，查詢API
                    # logger.info("API查詢基因: {}".format(single_gene_id))
                    species_name = self.query_gene_from_api(single_gene_id)
                    results.append((single_gene_id, species_name))
            
            return results  # 返回 [(gene_id, species_name), ...] 列表
        else:
            # 單一基因ID的處理
            # 先檢查緩存
            species_name = self.get_gene_species(gene_id)
            if species_name is not None:
                return species_name
            
            # 緩存中沒有，查詢API
            # logger.info("API查詢基因: {}".format(gene_id))
            return self.query_gene_from_api(gene_id)
    
    def save_cache(self):
        """保存緩存到檔案"""
        if self.new_entries > 0:
            try:
                with open(self.cache_file, 'w', encoding='utf-8') as f:
                    json.dump(self.cache, f, ensure_ascii=False, indent=2)
                logger.info("保存緩存: 新增 {} 個基因，總計 {} 個".format(self.new_entries, len(self.cache)))
            except Exception as e:
                logger.error("保存緩存失敗: {}".format(e))

class PMIDProcessor:
    """PMID處理主類"""
    
    def __init__(self):
        self.gene_cache = GeneCache()
        self.timestamp = datetime.datetime.now().strftime("%Y%m%d")
        self.output_files = self._init_output_files()
        self.json_file_map = {}  # 新增這行
        self._build_json_file_map()  # 新增這行
        self.stats = {
            'total_pmids': 0,
            'processed_pmids': 0,
            'failed_pmids': 0,  # 新增失敗計數
            'json_reads': 0,
            'api_queries': 0,
            'total_gene_mentions': 0,
            'start_time': time.time()
        }
    
    def _init_output_files(self):
        """初始化輸出檔案"""
        files = {
            'all_genes': os.path.join(OUTPUT_DIR, "{}_Results_all_genes.tsv".format(self.timestamp)),
            'chemicals': os.path.join(OUTPUT_DIR, "{}_Results_chemicals.tsv".format(self.timestamp)),
            'gene_descriptions': os.path.join(OUTPUT_DIR, "{}_Results_gene_descriptions.tsv".format(self.timestamp)),
            'species': os.path.join(OUTPUT_DIR, "{}_Results_species.tsv".format(self.timestamp)),
            'relations': os.path.join(OUTPUT_DIR, "{}_Results_gene_relations.tsv".format(self.timestamp))
        }
        
        # 寫入標題行
        headers = {
            'all_genes': ['pmid', 'gene_id', 'gene_name', 'gene_text', 'gene_position', 'passage_type', 
                         'homologene_id', 'species_name', 'original_gene_ids'],
            'chemicals': ['pmid', 'chemical_id', 'chemical_name', 'chemical_text', 'passage_type', 
                         'database', 'context'],
            'gene_descriptions': ['pmid', 'gene_id', 'gene_name', 'gene_text', 'species_name', 'context', 'passage_type', 'original_gene_ids'],
            'species': ['pmid', 'species_id', 'species_name', 'species_text', 'passage_type', 'mention_location'],
            'relations': ['pmid', 'relation_id', 'relation_type', 'score', 'role1_id', 'role1_type', 
                         'role1_name', 'role1_species_name', 'role2_id', 'role2_type', 'role2_name', 'role2_species_name']
        }
        
        for file_type, filepath in files.items():
            with open(filepath, 'w', encoding='utf-8') as f:
                f.write('\t'.join(headers[file_type]) + '\n')
        
        logger.info("初始化輸出檔案: {} 個TSV檔案".format(len(files)))
        return files
        
    def _build_json_file_map(self):
        """建立PMID到JSON檔案路徑的映射"""
        logger.info("建立JSON檔案映射...")
        
        if not os.path.exists(JSON_BASE_DIR):
            logger.error("JSON基礎目錄不存在: {}".format(JSON_BASE_DIR))
            return
        
        total_files = 0
        
        # 遍歷所有子目錄
        for year_dir in os.listdir(JSON_BASE_DIR):
            year_path = os.path.join(JSON_BASE_DIR, year_dir)
            if not os.path.isdir(year_path):
                continue
            
            # 遍歷年份目錄下的所有JSON檔案
            for filename in os.listdir(year_path):
                if filename.endswith('.json'):
                    # 從檔名提取PMID（格式：timestamp_pmid.json）
                    parts = filename.split('_')
                    if len(parts) >= 2:
                        pmid = parts[-1].replace('.json', '')
                        full_path = os.path.join(year_path, filename)
                        
                        # 如果同一個PMID有多個檔案，選擇最新的
                        if pmid in self.json_file_map:
                            existing_path = self.json_file_map[pmid]
                            if os.path.getmtime(full_path) > os.path.getmtime(existing_path):
                                self.json_file_map[pmid] = full_path
                        else:
                            self.json_file_map[pmid] = full_path
                        
                        total_files += 1
        
        logger.info("JSON檔案映射完成: {} 個檔案，{} 個唯一PMID".format(total_files, len(self.json_file_map)))        
    
    def read_pmid_list(self, pmid_file=PMID_FILE):
        """讀取PMID列表（從TSV檔案）"""
        try:
            pmids = []
            with open(pmid_file, 'r', encoding='utf-8') as f:
                # 跳過標題行
                next(f)
                for line in f:
                    parts = line.strip().split('\t')
                    if parts and parts[0].strip():
                        pmids.append(parts[0].strip())
            
            logger.info("讀取PMID列表: {} 個".format(len(pmids)))
            self.stats['total_pmids'] = len(pmids)
            return pmids
        except Exception as e:
            logger.error("讀取PMID列表失敗: {}".format(e))
            return []
    
    def find_json_file(self, pmid):
        """從映射中尋找對應的JSON檔案"""
        return self.json_file_map.get(str(pmid))
    
    def load_json_data(self, json_file):
        """載入JSON資料"""
        try:
            with open(json_file, 'r', encoding='utf-8') as f:
                data = json.load(f)
            self.stats['json_reads'] += 1
            return data
        except Exception as e:
            logger.error("讀取JSON失敗 ({}): {}".format(json_file, e))
            return None
    
    # 修改 process_pmid 方法
    def process_pmid(self, pmid):
        """處理單個PMID"""
        start_time = time.time()
        
        # 尋找JSON檔案
        json_file = self.find_json_file(pmid)
        if not json_file:
            # logger.warning("找不到PMID {} 的JSON檔案".format(pmid))  # 註解掉這行
            self.stats['failed_pmids'] += 1  # 新增
            return False
        
        # 載入JSON資料
        data = self.load_json_data(json_file)
        if not data:
            self.stats['failed_pmids'] += 1  # 新增
            return False
        
        # logger.info("讀取JSON: {}".format(os.path.basename(json_file)))  # 也可以註解掉這行
        
        # 處理資料
        api_queries_before = self.gene_cache.new_entries
        results = self.process_pubtator_data(data, pmid)
        api_queries_after = self.gene_cache.new_entries
        
        new_api_queries = api_queries_after - api_queries_before
        if new_api_queries > 0:
            self.stats['api_queries'] += new_api_queries
            # logger.info("API查詢: {} 個新基因".format(new_api_queries))  # 也可以註解掉這行
        
        # 寫入結果
        self.write_results(results)
        
        processing_time = time.time() - start_time
        self.stats['processed_pmids'] += 1
        
        # logger.info("PMID {} 處理完成 (耗時: {:.2f}秒)".format(pmid, processing_time))  # 註解掉這行
        
        # 每100個PMID顯示進度
        if self.stats['processed_pmids'] % 10000 == 0:
            self.show_progress()
        
        return True
    
    def process_pubtator_data(self, data, pmid):
        """處理PubTator3資料"""
        articles = data.get('PubTator3', [])
        
        gene_data = []
        species_data = []
        gene_descriptions = []
        relations_data = []
        chemicals_data = []
        
        for article in articles:
            article_pmid = article.get('pmid', '')
            if str(article_pmid) != str(pmid):
                continue
            
            # 處理物種
            for passage in article.get('passages', []):
                for annotation in passage.get('annotations', []):
                    info = annotation.get('infons', {})
                    
                    if info.get('type') == 'Species':
                        species_id = str(info.get('identifier', ''))
                        
                        species_info = {
                            'pmid': article_pmid,
                            'species_id': species_id,
                            'species_name': info.get('name', ''),
                            'species_text': annotation.get('text', ''),
                            'passage_type': passage.get('infons', {}).get('type', ''),
                            'mention_location': "{}:{}".format(
                                passage.get('offset', 0),
                                annotation.get('locations', [{}])[0].get('offset', 0) if annotation.get('locations') else 0
                            )
                        }
                        species_data.append(species_info)
            
            # 處理化學物質
            for passage in article.get('passages', []):
                passage_text = passage.get('text', '')
                passage_type = passage.get('infons', {}).get('type', '')
                offset = passage.get('offset', 0)
                
                for annotation in passage.get('annotations', []):
                    info = annotation.get('infons', {})
                    
                    if info.get('type') == 'Chemical':
                        chemical_id = str(info.get('identifier', ''))
                        chemical_name = info.get('name', '')
                        chemical_text = annotation.get('text', '')
                        
                        # 提取上下文
                        loc = annotation.get('locations', [{}])[0] if annotation.get('locations') else {}
                        start = loc.get('offset', 0) - offset
                        length = loc.get('length', len(chemical_text))
                        
                        context_start = max(0, start - 50)
                        context_end = min(len(passage_text), start + length + 50)
                        
                        context = passage_text[context_start:context_end]
                        context = context.replace(chemical_text, "**{}**".format(chemical_text))
                        
                        chemical_info = {
                            'pmid': article_pmid,
                            'chemical_id': chemical_id,
                            'chemical_name': chemical_name,
                            'chemical_text': chemical_text,
                            'passage_type': passage_type,
                            'database': info.get('database', ''),
                            'context': context
                        }
                        chemicals_data.append(chemical_info)
            
            # 處理基因
            for passage in article.get('passages', []):
                passage_text = passage.get('text', '')
                passage_type = passage.get('infons', {}).get('type', '')
                offset = passage.get('offset', 0)
                
                for annotation in passage.get('annotations', []):
                    info = annotation.get('infons', {})
                    
                    if info.get('type') == 'Gene':
                        original_gene_id = str(info.get('identifier', ''))
                        gene_name = info.get('name', '')
                        gene_text = annotation.get('text', '')
                        homologene_id = info.get('ncbi_homologene', '')
                        
                        # 跳過空的基因ID
                        if not original_gene_id or original_gene_id in ['None', 'none', '', 'null']:
                            continue
                        
                        # 查詢基因物種（可能返回單一結果或多基因結果）
                        species_result = self.gene_cache.get_gene_species_name(original_gene_id)
                        
                        # 提取基因上下文
                        loc = annotation.get('locations', [{}])[0] if annotation.get('locations') else {}
                        start = loc.get('offset', 0) - offset
                        length = loc.get('length', len(gene_text))
                        
                        context_start = max(0, start - 50)
                        context_end = min(len(passage_text), start + length + 50)
                        
                        context = passage_text[context_start:context_end]
                        context = context.replace(gene_text, "**{}**".format(gene_text))
                        
                        # 處理多基因ID的情況
                        if isinstance(species_result, list):
                            # 多基因ID，拆分成多筆記錄
                            for gene_id, species_name in species_result:
                                self.stats['total_gene_mentions'] += 1
                                
                                gene_info = {
                                    'pmid': article_pmid,
                                    'gene_id': gene_id,
                                    'gene_name': gene_name,
                                    'gene_text': gene_text,
                                    'gene_position': "{}:{}".format(
                                        offset,
                                        annotation.get('locations', [{}])[0].get('offset', 0) if annotation.get('locations') else 0
                                    ),
                                    'passage_type': passage_type,
                                    'homologene_id': homologene_id,
                                    'species_name': species_name,
                                    'original_gene_ids': original_gene_id  # 原始的多基因ID
                                }
                                gene_data.append(gene_info)
                                
                                desc = {
                                    'pmid': article_pmid,
                                    'gene_id': gene_id,
                                    'gene_name': gene_name,
                                    'gene_text': gene_text,
                                    'species_name': species_name,
                                    'context': context,
                                    'passage_type': passage_type,
                                    'original_gene_ids': original_gene_id  # 原始的多基因ID
                                }
                                gene_descriptions.append(desc)
                        else:
                            # 單基因ID
                            self.stats['total_gene_mentions'] += 1
                            
                            gene_info = {
                                'pmid': article_pmid,
                                'gene_id': original_gene_id,
                                'gene_name': gene_name,
                                'gene_text': gene_text,
                                'gene_position': "{}:{}".format(
                                    offset,
                                    annotation.get('locations', [{}])[0].get('offset', 0) if annotation.get('locations') else 0
                                ),
                                'passage_type': passage_type,
                                'homologene_id': homologene_id,
                                'species_name': species_result,
                                'original_gene_ids': original_gene_id  # 與gene_id相同
                            }
                            gene_data.append(gene_info)
                            
                            desc = {
                                'pmid': article_pmid,
                                'gene_id': original_gene_id,
                                'gene_name': gene_name,
                                'gene_text': gene_text,
                                'species_name': species_result,
                                'context': context,
                                'passage_type': passage_type,
                                'original_gene_ids': original_gene_id  # 與gene_id相同
                            }
                            gene_descriptions.append(desc)
            
            # 處理關係
            for relation in article.get('relations', []):
                rel_id = relation.get('id', '')
                infons = relation.get('infons', {})
                rel_type = infons.get('type', '')
                
                role1 = infons.get('role1', {})
                role2 = infons.get('role2', {})
                
                # 檢查角色是否是基因，並獲取物種資訊
                role1_species_name = ''
                if role1.get('type') == 'Gene':
                    role1_id = str(role1.get('identifier', ''))
                    role1_species_name = self.gene_cache.get_gene_species_name(role1_id)
                
                role2_species_name = ''
                if role2.get('type') == 'Gene':
                    role2_id = str(role2.get('identifier', ''))
                    role2_species_name = self.gene_cache.get_gene_species_name(role2_id)
                
                rel_info = {
                    'pmid': article_pmid,
                    'relation_id': rel_id,
                    'relation_type': rel_type,
                    'score': infons.get('score', ''),
                    'role1_id': role1.get('identifier', ''),
                    'role1_type': role1.get('type', ''),
                    'role1_name': role1.get('name', ''),
                    'role1_species_name': role1_species_name,
                    'role2_id': role2.get('identifier', ''),
                    'role2_type': role2.get('type', ''),
                    'role2_name': role2.get('name', ''),
                    'role2_species_name': role2_species_name
                }
                relations_data.append(rel_info)
        
        return {
            'all_genes': gene_data,
            'chemicals': chemicals_data,
            'gene_descriptions': gene_descriptions,
            'species': species_data,
            'relations': relations_data
        }
    
    def write_results(self, results):
        """寫入結果到TSV檔案"""
        for data_type, data_list in results.items():
            if data_list and data_type in self.output_files:
                filepath = self.output_files[data_type]
                
                with open(filepath, 'a', encoding='utf-8') as f:
                    for item in data_list:
                        # 根據不同的資料類型寫入不同的欄位，確保所有值都轉為字串並清理特殊字符
                        if data_type == 'all_genes':
                            row = [
                                str(item.get('pmid', '')),
                                str(item.get('gene_id', '')),
                                str(item.get('gene_name', '')).replace('\t', ' ').replace('\n', ' ').replace('\r', ' '),
                                str(item.get('gene_text', '')).replace('\t', ' ').replace('\n', ' ').replace('\r', ' '),
                                str(item.get('gene_position', '')),
                                str(item.get('passage_type', '')).replace('\t', ' ').replace('\n', ' ').replace('\r', ' '),
                                str(item.get('homologene_id', '')),
                                str(item.get('species_name', '')).replace('\t', ' ').replace('\n', ' ').replace('\r', ' '),
                                str(item.get('original_gene_ids', ''))
                            ]
                        elif data_type == 'chemicals':
                            row = [
                                str(item.get('pmid', '')),
                                str(item.get('chemical_id', '')),
                                str(item.get('chemical_name', '')).replace('\t', ' ').replace('\n', ' ').replace('\r', ' '),
                                str(item.get('chemical_text', '')).replace('\t', ' ').replace('\n', ' ').replace('\r', ' '),
                                str(item.get('passage_type', '')).replace('\t', ' ').replace('\n', ' ').replace('\r', ' '),
                                str(item.get('database', '')),
                                str(item.get('context', '')).replace('\t', ' ').replace('\n', ' ').replace('\r', ' ')
                            ]
                        elif data_type == 'gene_descriptions':
                            row = [
                                str(item.get('pmid', '')),
                                str(item.get('gene_id', '')),
                                str(item.get('gene_name', '')).replace('\t', ' ').replace('\n', ' ').replace('\r', ' '),
                                str(item.get('gene_text', '')).replace('\t', ' ').replace('\n', ' ').replace('\r', ' '),
                                str(item.get('species_name', '')).replace('\t', ' ').replace('\n', ' ').replace('\r', ' '),
                                str(item.get('context', '')).replace('\t', ' ').replace('\n', ' ').replace('\r', ' '),
                                str(item.get('passage_type', '')).replace('\t', ' ').replace('\n', ' ').replace('\r', ' '),
                                str(item.get('original_gene_ids', ''))
                            ]
                        elif data_type == 'species':
                            row = [
                                str(item.get('pmid', '')),
                                str(item.get('species_id', '')),
                                str(item.get('species_name', '')).replace('\t', ' ').replace('\n', ' ').replace('\r', ' '),
                                str(item.get('species_text', '')).replace('\t', ' ').replace('\n', ' ').replace('\r', ' '),
                                str(item.get('passage_type', '')).replace('\t', ' ').replace('\n', ' ').replace('\r', ' '),
                                str(item.get('mention_location', ''))
                            ]
                        elif data_type == 'relations':
                            row = [
                                str(item.get('pmid', '')),
                                str(item.get('relation_id', '')),
                                str(item.get('relation_type', '')).replace('\t', ' ').replace('\n', ' ').replace('\r', ' '),
                                str(item.get('score', '')),
                                str(item.get('role1_id', '')),
                                str(item.get('role1_type', '')).replace('\t', ' ').replace('\n', ' ').replace('\r', ' '),
                                str(item.get('role1_name', '')).replace('\t', ' ').replace('\n', ' ').replace('\r', ' '),
                                str(item.get('role1_species_name', '')).replace('\t', ' ').replace('\n', ' ').replace('\r', ' '),
                                str(item.get('role2_id', '')),
                                str(item.get('role2_type', '')).replace('\t', ' ').replace('\n', ' ').replace('\r', ' '),
                                str(item.get('role2_name', '')).replace('\t', ' ').replace('\n', ' ').replace('\r', ' '),
                                str(item.get('role2_species_name', '')).replace('\t', ' ').replace('\n', ' ').replace('\r', ' ')
                            ]
                        
                        f.write('\t'.join(row) + '\n')
    
    def show_progress(self):
        """顯示處理進度"""
        elapsed = time.time() - self.stats['start_time']
        rate = self.stats['processed_pmids'] / elapsed if elapsed > 0 else 0
        remaining = self.stats['total_pmids'] - self.stats['processed_pmids'] - self.stats['failed_pmids']
        eta = remaining / rate if rate > 0 else 0
        
        logger.info("進度: {:,}/{:,} ({:.1f}%), 速度: {:.1f} PMID/秒, 預計剩餘: {:.1f}小時, 失敗: {:,}".format(
                   self.stats['processed_pmids'], self.stats['total_pmids'], 
                   (self.stats['processed_pmids'] + self.stats['failed_pmids'])/self.stats['total_pmids']*100, 
                   rate, eta/3600, self.stats['failed_pmids']))
        logger.info("統計: JSON讀取 {:,}, API查詢 {:,}, 基因註釋筆數 {:,}".format(
                   self.stats['json_reads'], self.stats['api_queries'], self.stats['total_gene_mentions']))
    
    def create_statistics_summary(self):
        """創建統計摘要Excel"""
        try:
            stats_file = os.path.join(OUTPUT_DIR, "{}_Statistics.xlsx".format(self.timestamp))
            
            # 統計各個資料集的基本信息
            stats_data = []
            for data_type, filepath in self.output_files.items():
                if os.path.exists(filepath):
                    try:
                        # 只讀取前幾行來獲取基本統計
                        df_sample = pd.read_csv(filepath, sep='\t', nrows=1000, 
                                              on_bad_lines='skip', engine='python', quoting=3)
                        
                        # 計算總行數（不讀取全部資料）
                        with open(filepath, 'r', encoding='utf-8') as f:
                            total_lines = sum(1 for _ in f) - 1  # 減去標題行
                        
                        stats_data.append({
                            'Dataset': data_type,
                            'Total_Rows': total_lines,
                            'Columns': len(df_sample.columns),
                            'File_Size_MB': round(os.path.getsize(filepath) / 1024 / 1024, 2)
                        })
                        
                    except Exception as e:
                        logger.warning("統計 {} 時出現問題: {}".format(data_type, e))
            
            if stats_data:
                stats_df = pd.DataFrame(stats_data)
                
                # 添加處理統計
                process_stats = pd.DataFrame([{
                    'Metric': 'Total PMIDs',
                    'Value': self.stats['total_pmids']
                }, {
                    'Metric': 'Processed PMIDs', 
                    'Value': self.stats['processed_pmids']
                }, {
                    'Metric': 'JSON Reads',
                    'Value': self.stats['json_reads']
                }, {
                    'Metric': 'API Queries',
                    'Value': self.stats['api_queries']
                }, {
                    'Metric': 'Total Genes',
                    'Value': self.stats['total_gene_mentions']
                }, {
                    'Metric': 'Cache Size',
                    'Value': len(self.gene_cache.cache)
                }])
                
                with pd.ExcelWriter(stats_file) as writer:
                    stats_df.to_excel(writer, sheet_name='Dataset_Statistics', index=False)
                    process_stats.to_excel(writer, sheet_name='Process_Statistics', index=False)
                
                logger.info("創建統計摘要: {}".format(stats_file))
                
        except Exception as e:
            logger.error("創建統計摘要失敗: {}".format(e))
    
    def run(self):
        """主執行函數"""
        logger.info("開始PMID批量處理")
        
        # 讀取PMID列表
        pmids = self.read_pmid_list()
        if not pmids:
            logger.error("沒有PMID需要處理")
            return
        
        # 處理每個PMID
        success_count = 0
        for i, pmid in enumerate(pmids, 1):
            # logger.info("處理進度: {}/{} - PMID: {}".format(i, len(pmids), pmid))
            
            if self.process_pmid(pmid):
                success_count += 1
            # else:
                # logger.warning("PMID {} 處理失敗".format(pmid))
        
        # 保存緩存
        self.gene_cache.save_cache()
        
        # 創建TSV統計摘要
        self.create_tsv_statistics()
        
        # 最終統計
        elapsed = time.time() - self.stats['start_time']
        logger.info("處理完成！")
        logger.info("最終統計:")
        logger.info("   - 總PMID: {}".format(self.stats['total_pmids']))
        logger.info("   - 成功處理: {}".format(success_count))
        logger.info("   - 失敗: {}".format(self.stats['failed_pmids']))  # 修改這行
        logger.info("   - 總耗時: {:.2f}秒 ({:.2f}小時)".format(elapsed, elapsed/3600))
        logger.info("   - 平均速度: {:.2f} PMID/秒".format(self.stats['processed_pmids']/elapsed))
        logger.info("   - JSON讀取: {}".format(self.stats['json_reads']))
        logger.info("   - API查詢: {}".format(self.stats['api_queries']))
        logger.info("   - 基因註釋筆數: {:,}".format(self.stats['total_gene_mentions']))
        
    def create_tsv_statistics(self):
        """創建TSV格式的統計摘要"""
        try:
            stats_file = os.path.join(OUTPUT_DIR, "{}_Statistics.tsv".format(self.timestamp))
            
            with open(stats_file, 'w', encoding='utf-8') as f:
                # 寫入標題
                f.write("Dataset\tTotal_Rows\tColumns\tFile_Size_MB\n")
                
                # 統計各個資料集
                for data_type, filepath in self.output_files.items():
                    if os.path.exists(filepath):
                        try:
                            # 計算總行數
                            with open(filepath, 'r', encoding='utf-8') as data_f:
                                total_lines = sum(1 for _ in data_f) - 1  # 減去標題行
                            
                            # 讀取標題行計算欄位數
                            with open(filepath, 'r', encoding='utf-8') as data_f:
                                header = data_f.readline().strip()
                                columns = len(header.split('\t'))
                            
                            file_size_mb = round(os.path.getsize(filepath) / 1024 / 1024, 2)
                            
                            f.write("{}\t{}\t{}\t{}\n".format(
                                data_type, total_lines, columns, file_size_mb))
                            
                        except Exception as e:
                            logger.warning("統計 {} 時出現問題: {}".format(data_type, e))
                
                # 添加處理統計
                f.write("\n# Process Statistics\n")
                f.write("Metric\tValue\n")
                f.write("Total_PMIDs\t{}\n".format(self.stats['total_pmids']))
                f.write("Processed_PMIDs\t{}\n".format(self.stats['processed_pmids']))
                f.write("JSON_Reads\t{}\n".format(self.stats['json_reads']))
                f.write("API_Queries\t{}\n".format(self.stats['api_queries']))
                f.write("Total_Gene_Mentions\t{}\n".format(self.stats['total_gene_mentions']))
                f.write("Cache_Size\t{}\n".format(len(self.gene_cache.cache)))
            
            logger.info("創建TSV統計摘要: {}".format(stats_file))
            
        except Exception as e:
            logger.error("創建TSV統計摘要失敗: {}".format(e))            

def main():
    """主函數"""
    if len(sys.argv) > 1:
        print("用法: python pmid_batch_processor.py")
        print("此程式會自動讀取 circadian_all_pmid/circadian_complete_pmids_20250802_195005.tsv 中的所有PMID進行批量處理")
        return
    
    processor = PMIDProcessor()
    processor.run()

if __name__ == "__main__":
    main()