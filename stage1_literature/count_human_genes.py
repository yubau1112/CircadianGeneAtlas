import os
import glob
import pandas as pd
from collections import defaultdict
import sys
import locale
import csv  # 加在其他 import 之後
import json  # 加在其他 import 之後

def load_gene_cache(cache_file="gene_species_cache.json"):
    """載入基因物種緩存"""
    if os.path.exists(cache_file):
        try:
            with open(cache_file, 'r', encoding='utf-8') as f:
                cache = json.load(f)
            print(f"✓ 載入基因緩存: {len(cache):,} 個基因")
            return cache
        except Exception as e:
            print(f"✗ 載入緩存失敗: {e}")
            return {}
    else:
        print(f"✗ 找不到緩存檔案: {cache_file}")
        return {}

def get_species_from_cache(gene_id, gene_cache):
    """從緩存中獲取基因物種，優先使用物種ID判斷"""
    gene_id = str(gene_id)
    if gene_id in gene_cache:
        cached_info = gene_cache[gene_id]
        # 如果是人類基因 (NCBI ID: 9606)，直接返回標準名稱
        if cached_info.get('species') == '9606':
            return 'Homo sapiens'
        return cached_info.get('species_name', 'unknown')
    return 'unknown'

def count_human_genes_per_article():
    # 設定千位分隔符
    locale.setlocale(locale.LC_ALL, '')
    
    # 載入基因緩存
    gene_cache = load_gene_cache()
    
    print("=" * 60)
    print("開始執行人類基因統計分析")
    print("=" * 60)
    
    # 步驟1: 讀取 pmid-circadianTitle-set.txt 並統計 PMID 數量
    print("\n步驟1: 讀取 PMID 清單")
    print("-" * 30)
    pmid_file = 'pmid-circadianTitle-set.txt'
    
    if os.path.exists(pmid_file):
        try:
            with open(pmid_file, 'r', encoding='utf-8') as f:
                pmid_lines = f.readlines()
            # 移除空行和只有空白字元的行
            pmid_list = [line.strip() for line in pmid_lines if line.strip()]
            total_pmids = len(pmid_list)
            print(f"✓ PMID 清單檔案: {pmid_file}")
            print(f"✓ 總 PMID 數: {format(total_pmids, ',')}")
        except Exception as e:
            print(f"✗ 無法讀取 PMID 清單檔案: {e}")
            total_pmids = 0
            pmid_list = []
    else:
        print(f"✗ 找不到 PMID 清單檔案: {pmid_file}")
        total_pmids = 0
        pmid_list = []
    
    # 步驟2: 統計 JSON 檔案數並找出缺失的 PMID
    print("\n步驟2: 統計 JSON 檔案並檢查缺失")
    print("-" * 30)
    
    # 修改：搜尋分年份資料夾中的所有 JSON 檔案
    json_base_path = 'circadian_all_pmid/json_files'
    
    # 方法1: 使用遞迴搜尋所有年份資料夾
    json_files = glob.glob(f'{json_base_path}/*/*.json')
    
    # 方法2: 如果方法1沒找到檔案，嘗試其他可能的路徑
    if not json_files:
        # 嘗試不同的路徑結構
        alternative_paths = [
            'json_files/*/*.json',  # 原始路徑下的年份資料夾
            '*/json_files/*/*.json'  # 更上層的結構
        ]
        
        for path_pattern in alternative_paths:
            json_files = glob.glob(path_pattern)
            if json_files:
                print(f"✓ 在路徑找到檔案: {path_pattern}")
                break
    
    total_json_files = len(json_files)
    print(f"✓ JSON 檔案搜尋路徑: {json_base_path}")
    print(f"✓ JSON 檔案數: {format(total_json_files, ',')}")
    
    # 如果找不到檔案，顯示診斷資訊
    if total_json_files == 0:
        print("⚠ 找不到 JSON 檔案，檢查資料夾結構...")
        
        # 檢查基礎路徑是否存在
        if os.path.exists(json_base_path):
            print(f"✓ 基礎路徑存在: {json_base_path}")
            
            # 列出年份資料夾
            year_folders = [d for d in os.listdir(json_base_path) 
                          if os.path.isdir(os.path.join(json_base_path, d))]
            print(f"✓ 找到年份資料夾: {len(year_folders)} 個")
            print(f"✓ 年份資料夾範例: {sorted(year_folders)[:10]}")
            
            # 檢查第一個年份資料夾中的檔案
            if year_folders:
                first_year = sorted(year_folders)[0]
                year_path = os.path.join(json_base_path, first_year)
                json_files_in_year = glob.glob(os.path.join(year_path, '*.json'))
                print(f"✓ {first_year} 年資料夾中的 JSON 檔案數: {len(json_files_in_year)}")
                
                if json_files_in_year:
                    print(f"✓ 範例檔案: {os.path.basename(json_files_in_year[0])}")
        else:
            print(f"✗ 基礎路徑不存在: {json_base_path}")
            
            # 嘗試尋找可能的路徑
            possible_paths = [
                'json_files',
                'circadian_all_pmid',
                './json_files',
                '../json_files'
            ]
            
            for path in possible_paths:
                if os.path.exists(path):
                    print(f"✓ 找到可能的路徑: {path}")
        
        # 如果沒有找到檔案就停止
        if total_json_files == 0:
            print("✗ 請確認 JSON 檔案路徑是否正確")
            return
    
    # 顯示按年份分布的統計
    if total_json_files > 0:
        year_distribution = defaultdict(int)
        for json_file in json_files:
            # 從路徑中提取年份
            path_parts = json_file.split(os.sep)
            if len(path_parts) >= 2:
                year = path_parts[-2]  # 倒數第二個部分應該是年份
                year_distribution[year] += 1
        
        year_dist_list = [f"{year}: {format(count, ',')}" for year, count in sorted(year_distribution.items())]
        print(f"\n✓ JSON 檔案年份分布: {', '.join(year_dist_list)}")
    
    # 提取所有 JSON 檔案的 PMID (檔名去掉.json，並處理時間戳記前綴)
    json_pmids = set()
    for json_file in json_files:
        filename = os.path.basename(json_file)
        # 移除 .json 副檔名
        pmid_part = filename.replace('.json', '')
        
        # 處理時間戳記前綴 (格式: YYYYMMDDHHMMSS_PMID)
        if '_' in pmid_part:
            # 取最後一個底線後的部分作為 PMID
            pmid = pmid_part.split('_')[-1]
        else:
            # 如果沒有底線，直接使用整個檔名作為 PMID
            pmid = pmid_part
        
        json_pmids.add(pmid)
    
    print(f"✓ 從檔名提取的 PMID 數: {format(len(json_pmids), ',')}")
    
    # 調試：檢查 PMID 格式
    if total_pmids > 0 and len(json_pmids) > 0:
        print(f"✓ 原始 PMID 範例 (前5個): {pmid_list[:5]}")
        print(f"✓ 提取的 JSON PMID 範例 (前5個): {list(sorted(json_pmids))[:5]}")
    
    # 比較 PMID 清單和 JSON 檔案
    if total_pmids > 0 and total_json_files > 0:
        completion_rate = round(total_json_files / total_pmids * 100, 2)
        print(f"✓ 完成率: {completion_rate}% ({total_json_files}/{total_pmids})")
        
        # 確保 PMID 格式一致 (都轉為字串並去除空白)
        pmid_set = set(str(pmid).strip() for pmid in pmid_list) if total_pmids > 0 else set()
        json_pmids_cleaned = set(str(pmid).strip() for pmid in json_pmids)
        
        print(f"✓ 清理後原始 PMID 數: {format(len(pmid_set), ',')}")
        print(f"✓ 清理後 JSON PMID 數: {format(len(json_pmids_cleaned), ',')}")
        
        # 找出缺失的 PMID
        missing_pmids = pmid_set - json_pmids_cleaned
        extra_pmids = json_pmids_cleaned - pmid_set
        
        # 調試：檢查交集
        intersection = pmid_set & json_pmids_cleaned
        print(f"✓ 匹配的 PMID 數: {format(len(intersection), ',')}")
        print(f"✓ 匹配率: {round(len(intersection) / len(pmid_set) * 100, 2) if len(pmid_set) > 0 else 0}%")
        
        if missing_pmids:
            missing_count = len(missing_pmids)
            print(f"⚠ 缺少 JSON 檔案的 PMID 數: {format(missing_count, ',')}")
            
            # 將缺失的 PMID 寫入檔案
            missing_file = 'missing_pmids.txt'
            try:
                with open(missing_file, 'w', encoding='utf-8') as f:
                    f.write("# 缺少 JSON 檔案的 PMID (可能是書籍而非論文)\n")
                    f.write(f"# 總數: {missing_count}\n")
                    f.write(f"# 生成時間: {pd.Timestamp.now()}\n\n")
                    for pmid in sorted(missing_pmids):
                        f.write(f"{pmid}\n")
                print(f"✓ 缺失的 PMID 已寫入: {missing_file}")
                
                # 顯示前10個缺失的 PMID 作為範例
                sample_missing = list(sorted(missing_pmids))[:10]
                print(f"✓ 缺失 PMID 範例 (前10個): {', '.join(sample_missing)}")
                
            except Exception as e:
                print(f"✗ 無法寫入缺失 PMID 檔案: {e}")
        
        if extra_pmids:
            extra_count = len(extra_pmids)
            print(f"⚠ 多出的 JSON 檔案 PMID 數: {format(extra_count, ',')}")
            
            # 將多出的 PMID 寫入檔案
            extra_file = 'extra_pmids.txt'
            try:
                with open(extra_file, 'w', encoding='utf-8') as f:
                    f.write("# 多出的 JSON 檔案 PMID (不在原始清單中)\n")
                    f.write(f"# 總數: {extra_count}\n")
                    f.write(f"# 生成時間: {pd.Timestamp.now()}\n\n")
                    for pmid in sorted(extra_pmids):
                        f.write(f"{pmid}\n")
                print(f"✓ 多出的 PMID 已寫入: {extra_file}")
                
            except Exception as e:
                print(f"✗ 無法寫入多出 PMID 檔案: {e}")
        
        if not missing_pmids and not extra_pmids:
            print("✓ JSON 檔案與 PMID 清單完全匹配")
    
    # 步驟3: 尋找基因資料檔案
    print("\n步驟3: 尋找基因資料檔案")
    print("-" * 30)
    gene_files = glob.glob('pubtator_all_tsv_results/*_Results_all_genes.tsv')
    
    if not gene_files:
        print("✗ 找不到 pmid_batch_processor.py 生成的基因檔案")
        print("✗ 請確認檔案格式為: YYYYMMDD_Results_all_genes.tsv")
        return
    
    # 使用最新的檔案（如果有多個）
    gene_file = max(gene_files, key=os.path.getmtime)
    print(f"✓ 找到 {len(gene_files)} 個基因資料檔案")
    print(f"✓ 使用最新檔案: {gene_file}")
    
    # 步驟4: 讀取並分析基因資料
    print("\n步驟4: 讀取基因資料")
    print("-" * 30)
    
    # 先檢查檔案大小和行數
    try:
        import subprocess
        result = subprocess.run(['wc', '-l', gene_file], capture_output=True, text=True)
        if result.returncode == 0:
            total_lines = int(result.stdout.split()[0])
            print(f"✓ 檔案總行數 (含標題): {format(total_lines, ',')}")
        else:
            print("⚠ 無法取得檔案行數")
            total_lines = 0
    except Exception as e:
        print(f"⚠ 檔案行數檢查失敗: {e}")
        total_lines = 0
    
    try:
        print("✓ 嘗試讀取檔案...")
        
        # 方法1: 預設方式（修復引號問題）
        df = pd.read_csv(gene_file, sep='\t', low_memory=False, quoting=csv.QUOTE_NONE)
        
        # 方法2: 檢查是否有引號問題 (無引號模式)
        try:
            df2 = pd.read_csv(gene_file, sep='\t', low_memory=False, quoting=csv.QUOTE_NONE)  # 改成這個
            print(f"✓ 無引號模式記錄數: {format(len(df2), ',')}")
            if len(df2) > len(df):
                df = df2
                print("✓ 使用無引號模式")
        except Exception as e:
            print(f"⚠ 無引號模式失敗: {e}")
        
        # 方法3: 檢查是否有錯誤行 (跳過壞行模式)
        try:
            df3 = pd.read_csv(gene_file, sep='\t', low_memory=False, on_bad_lines='skip', quoting=csv.QUOTE_NONE)
            if len(df3) > len(df):
                df = df3
                print("✓ 使用跳過壞行模式")
        except Exception as e:
            print(f"⚠ 跳過壞行模式失敗: {e}")
        
        print(f"✓ 最終使用記錄數: {format(len(df), ',')}")
        
        # 計算遺失比例
        if total_lines > 0:
            expected_records = total_lines - 1  # 扣除標題行
            missing_records = expected_records - len(df)
            missing_percentage = round(missing_records / expected_records * 100, 2) if expected_records > 0 else 0
            print(f"⚠ 遺失記錄數: {format(missing_records, ',')} ({missing_percentage}%)")
        
        # 檢查檔案格式
        print(f"✓ 檔案欄位數: {len(df.columns)}")
        print(f"✓ 檔案欄位: {list(df.columns)}")
        
        # 檢查各欄位的空值情況
        null_counts = df.isnull().sum()
        if null_counts.sum() > 0:
            print("⚠ 發現空值:")
            for col, count in null_counts.items():
                if count > 0:
                    print(f"   {col}: {format(count, ',')} 個空值")
        else:
            print("✓ 無空值")
            
    except Exception as e:
        print(f"✗ 無法讀取檔案 {gene_file}: {e}")
        return
    
    # 檢查必要的欄位
    required_columns = ['pmid', 'gene_id', 'gene_name', 'species_name']
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        print(f"✗ 缺少必要欄位: {missing_columns}")
        print(f"✓ 可用欄位: {list(df.columns)}")
        return
    else:
        print(f"✓ 所有必要欄位都存在: {required_columns}")
    
    # 步驟5: 檢查API錯誤記錄
    print("\n步驟5: 檢查API錯誤記錄")
    print("-" * 30)
    
    # 檢查基因名稱欄位中的錯誤
    error_keywords = ['Error occurred', 'cannot get document summary', 'API error']
    error_mask = df['gene_name'].str.contains('|'.join(error_keywords), case=False, na=False)
    error_records = df[error_mask]
    
    print(f"✓ 總基因記錄數: {format(len(df), ',')}")
    print(f"⚠ API錯誤記錄數: {format(len(error_records), ',')}")
    if len(df) > 0:
        print(f"⚠ 錯誤比例: {round(len(error_records)/len(df)*100, 2)}%")
    
    # 統計錯誤記錄涉及的文章數
    if len(error_records) > 0:
        error_pmids = error_records['pmid'].nunique()
        print(f"⚠ 包含錯誤記錄的文章數: {format(error_pmids, ',')}")
        
        # 將錯誤記錄保存到檔案
        error_file = 'api_error_records.csv'
        try:
            error_records.to_csv(error_file, index=False)
            print(f"✓ 錯誤記錄已保存到: {error_file}")
        except Exception as e:
            print(f"⚠ 無法保存錯誤記錄: {e}")
    
    # 步驟6: 統計基因資料中的文章數
    print("\n步驟6: 分析基因資料覆蓋範圍")
    print("-" * 30)
    unique_pmids_in_gene_data = df['pmid'].nunique()
    print(f"✓ 基因資料中的獨特 PMID 數: {format(unique_pmids_in_gene_data, ',')}")
    
    if total_json_files > 0:
        coverage_rate = round(unique_pmids_in_gene_data / total_json_files * 100, 2)
        print(f"✓ 基因資料覆蓋率 (vs JSON檔案): {coverage_rate}%")
    
    if total_pmids > 0:
        original_coverage = round(unique_pmids_in_gene_data / total_pmids * 100, 2)
        print(f"✓ 基因資料覆蓋率 (vs 原始PMID): {original_coverage}%")
    
    # 步驟7: 過濾人類基因（使用基因緩存 - 優化版本）
    print("\n步驟7: 過濾人類基因（使用基因緩存）")
    print("-" * 30)

    # 優化版本 - 使用向量化操作
    if gene_cache:
        print("✓ 使用基因緩存補強物種資訊...")
        
        # 方法1: 使用 map() 函數進行向量化操作
        def get_cached_species(gene_id):
            gene_id_str = str(gene_id)
            if gene_id_str in gene_cache:
                cached_info = gene_cache[gene_id_str]
                if cached_info.get('species') == '9606':
                    return 'Homo sapiens'
                return cached_info.get('species_name', 'unknown')
            return None
        
        # 獲取緩存中的物種資訊
        print("✓ 正在查找緩存中的物種資訊...")
        cached_species = df['gene_id'].map(get_cached_species)
        
        # 只更新有緩存資訊的記錄
        update_mask = cached_species.notna()
        
        # 統計更新前的狀態
        before_update = len(df)
        human_before = (df['species_name'] == 'Homo sapiens').sum()
        
        # 執行批量更新
        print("✓ 正在批量更新物種資訊...")
        df.loc[update_mask, 'species_name'] = cached_species[update_mask]
        
        # 統計更新結果
        updated_count = update_mask.sum()
        human_after = (df['species_name'] == 'Homo sapiens').sum()
        new_human_genes = human_after - human_before
        
        print(f"✓ 更新了 {updated_count:,} 筆記錄的物種資訊")
        print(f"✓ 新增人類基因記錄: {new_human_genes:,} 筆")
        print(f"✓ 總人類基因記錄: {human_after:,} 筆")

    # 顯示物種分布
    species_counts = df['species_name'].value_counts().head(10)
    print("✓ 前10個物種分布:")
    for species, count in species_counts.items():
        print(f"   {species}: {format(count, ',')}")

    # 創建人類基因篩選條件 - 直接使用物種名稱
    print("✓ 正在篩選人類基因...")
    human_mask = (df['species_name'] == 'Homo sapiens')

    # 排除API錯誤的記錄
    error_keywords = ['Error occurred', 'cannot get document summary', 'API error']
    error_mask = df['gene_name'].str.contains('|'.join(error_keywords), case=False, na=False)
    valid_mask = ~error_mask

    # 原始人類基因（包含錯誤）
    all_human_genes_df = df[human_mask].copy()

    # 有效人類基因（排除錯誤）
    human_genes_df = df[human_mask & valid_mask].copy()

    print(f"\n✓ 人類基因記錄數 (含錯誤): {format(len(all_human_genes_df), ',')}")
    print(f"✓ 有效人類基因記錄數: {format(len(human_genes_df), ',')}")
    print(f"⚠ 人類基因中的錯誤記錄數: {format(len(all_human_genes_df) - len(human_genes_df), ',')}")

    if len(df) > 0:
        print(f"✓ 有效人類基因比例: {round(len(human_genes_df)/len(df)*100, 2)}%")
    
    # 統計包含人類基因的文章數
    articles_with_human_genes = all_human_genes_df['pmid'].nunique()
    articles_with_valid_human_genes = human_genes_df['pmid'].nunique()
    
    print(f"✓ 包含人類基因的文章數 (含錯誤): {format(articles_with_human_genes, ',')}")
    print(f"✓ 包含有效人類基因的文章數: {format(articles_with_valid_human_genes, ',')}")
    
    if unique_pmids_in_gene_data > 0:
        valid_human_gene_article_rate = round(articles_with_valid_human_genes / unique_pmids_in_gene_data * 100, 2)
        print(f"✓ 有效人類基因文章比例: {valid_human_gene_article_rate}%")
    
    # 步驟8: 處理基因頻率統計
    print("\n步驟8: 統計有效基因頻率")
    print("-" * 30)
    
    if len(human_genes_df) == 0:
        print("✗ 沒有有效的人類基因記錄可供統計")
        return
    
    # 創建字典來存儲結果
    result_dict = defaultdict(lambda: defaultdict(int))
    gene_symbols = {}  # 用於存儲 gene_id 對應的 gene_name
    
    # 處理每筆記錄，並顯示進度條
    total_records = len(human_genes_df)
    print(f"✓ 開始處理 {format(total_records, ',')} 筆有效人類基因記錄...")
    
    for i, (_, row) in enumerate(human_genes_df.iterrows()):
        # 計算並顯示進度
        if i % 10000 == 0 or i == total_records - 1:  # 每10000筆或最後一筆顯示進度
            progress = (i + 1) / total_records
            bar_length = 20
            bar = '=' * int(bar_length * progress) + ' ' * (bar_length - int(bar_length * progress))
            percent = round(progress * 100, 2)
            sys.stdout.write(f'\r   處理進度: [{bar}] {format(i+1, ",")}/{format(total_records, ",")} = {percent}%')
            sys.stdout.flush()
        
        pmid = str(row['pmid'])
        gene_id = str(row['gene_id'])
        gene_name = str(row['gene_name']) if pd.notna(row['gene_name']) else ""
        
        # 計算每個基因ID在每篇文章中的出現次數
        result_dict[pmid][gene_id] += 1
        
        # 存儲基因ID對應的基因符號 (gene_name)
        if gene_id not in gene_symbols and gene_name:
            gene_symbols[gene_id] = gene_name
    
    # 完成處理後換行
    print("\n✓ 基因頻率統計完成")
    
    # 步驟9: 準備輸出結果
    print("\n步驟9: 準備輸出結果")
    print("-" * 30)
    
    output_data = []
    for pmid, gene_counts in result_dict.items():
        for gene_id, count in gene_counts.items():
            # 獲取基因符號，如果沒有則用空字串
            gene_symbol = gene_symbols.get(gene_id, "")
            
            output_data.append({
                'pmid': pmid,
                'gene_id': gene_id,
                'gene_symbol': gene_symbol,
                'count': count
            })
    
    # 創建DataFrame並排序
    result_df = pd.DataFrame(output_data)
    
    if result_df.empty:
        print("✗ 沒有找到有效的人類基因資料")
        return
    
    # 按照count從大到小排序
    result_df = result_df.sort_values(by=['count'], ascending=False)
    
    # 最終統計
    print("\n" + "=" * 60)
    print("最終統計結果")
    print("=" * 60)
    print(f"✓ 總基因-文章對數: {format(len(result_df), ',')}")
    print(f"✓ 涉及文章數: {format(result_df['pmid'].nunique(), ',')}")
    print(f"✓ 涉及基因數: {format(result_df['gene_id'].nunique(), ',')}")
    print(f"✓ 平均每篇文章的基因數: {round(len(result_df) / result_df['pmid'].nunique(), 2)}")
    print(f"✓ 平均每個基因出現在文章數: {round(len(result_df) / result_df['gene_id'].nunique(), 2)}")
    
    # 在輸出到CSV之前，顯示前20筆資料
    print(f"\n前20筆結果 (按出現次數排序):")
    print("-" * 50)
    print(result_df.head(20).to_string(index=False))
    
    # 輸出到CSV檔案
    output_file = 'pmid_gene_frequency.csv'
    result_df.to_csv(output_file, index=False)
    print(f"\n✓ 結果已儲存到: {output_file}")
    
    # 額外統計：最常出現的基因
    print(f"\n最常出現的前10個基因 (跨所有文章):")
    print("-" * 50)
    gene_total_counts = result_df.groupby(['gene_id', 'gene_symbol'])['count'].sum().reset_index()
    gene_total_counts = gene_total_counts.sort_values(by='count', ascending=False)
    print(gene_total_counts.head(10).to_string(index=False))
    
    print("\n" + "=" * 60)
    print("分析完成！")
    print("API錯誤記錄已被自動排除，最終結果不包含錯誤資料")
    print("=" * 60)

if __name__ == "__main__":
    count_human_genes_per_article()