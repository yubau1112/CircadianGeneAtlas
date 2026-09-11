#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Cache to TSV Tool - English Version (Updated)
Extract results from circadian_genes_cache_dual.json and convert to TSV format
All output in English with ASCII-only characters for maximum compatibility
Default: Full NCBI summary output (no truncation)
"""

import json
import csv
import argparse
import os
import sys
from datetime import datetime

def load_cache_file(cache_file):
    """Load cache file"""
    try:
        with open(cache_file, 'r', encoding='utf-8') as f:
            cache_data = json.load(f)
        print(f"Successfully loaded cache file with {len(cache_data)} gene records")
        return cache_data
    except FileNotFoundError:
        print(f"Error: Cache file not found - {cache_file}")
        return None
    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON format in cache file - {e}")
        return None
    except Exception as e:
        print(f"Error: Problem reading cache file - {e}")
        return None

def clean_text_for_tsv(text, max_length=None):
    """Clean text for TSV format - keep only ASCII characters, optional truncation"""
    if not text:
        return ""
    
    # Convert to string
    text = str(text)
    
    # Remove or replace problematic characters
    import re
    # Keep only English letters, numbers, basic punctuation and spaces
    cleaned = re.sub(r'[^\x00-\x7F]+', '', text)  # Remove all non-ASCII characters
    
    # Clean special characters
    cleaned = cleaned.replace('\n', ' ').replace('\t', ' ').replace('\r', ' ')
    cleaned = re.sub(r'\s+', ' ', cleaned)  # Merge multiple spaces into one
    
    # Optional truncation
    if max_length and len(cleaned) > max_length:
        cleaned = cleaned[:max_length] + "..."
    
    return cleaned.strip()

def extract_score_info_english(score_reasons):
    """Extract scoring information in English from Chinese score reasons - GENE VERSION"""
    if not score_reasons:
        return "No_scoring_info"
    
    extracted_info = []
    
    for reason in score_reasons:
        reason_str = str(reason)
        
        # Extract score from brackets
        import re
        score_match = re.search(r'\(\+(\d+(?:\.\d+)?)\)', reason_str)
        score = score_match.group(1) if score_match else "0"
        
        # 基因版本的匹配模式
        if "已知晝夜節律基因" in reason_str or "Known_circadian_gene" in reason_str:
            extracted_info.append(f"Known_circadian_gene_+{score}")
        elif "文獻核心關鍵字" in reason_str or "Literature_core" in reason_str:
            # Extract keyword after colon
            keyword_match = re.search(r'文獻核心關鍵字:\s*([^(]+)', reason_str)
            if keyword_match:
                keyword = keyword_match.group(1).strip()
                keyword = re.sub(r'\\b|\\\\b', '', keyword)  # Remove regex markers
                keyword = clean_text_for_tsv(keyword)
                extracted_info.append(f"Lit_core_{keyword}_+{score}")
            else:
                extracted_info.append(f"Lit_core_+{score}")
        elif "文獻通路關鍵字" in reason_str or "Literature_pathway" in reason_str:
            keyword_match = re.search(r'文獻通路關鍵字:\s*([^(]+)', reason_str)
            if keyword_match:
                keyword = keyword_match.group(1).strip()
                keyword = re.sub(r'\\b|\\\\b', '', keyword)
                keyword = clean_text_for_tsv(keyword)
                extracted_info.append(f"Lit_pathway_{keyword}_+{score}")
            else:
                extracted_info.append(f"Lit_pathway_+{score}")
        elif "文獻機制關鍵字" in reason_str or "Literature_mechanism" in reason_str:
            keyword_match = re.search(r'文獻機制關鍵字:\s*([^(]+)', reason_str)
            if keyword_match:
                keyword = keyword_match.group(1).strip()
                keyword = clean_text_for_tsv(keyword)
                extracted_info.append(f"Lit_mechanism_{keyword}_+{score}")
            else:
                extracted_info.append(f"Lit_mechanism_+{score}")
        elif "文獻疾病關鍵字" in reason_str or "Literature_disease" in reason_str:
            keyword_match = re.search(r'文獻疾病關鍵字:\s*([^(]+)', reason_str)
            if keyword_match:
                keyword = keyword_match.group(1).strip()
                keyword = clean_text_for_tsv(keyword)
                extracted_info.append(f"Lit_disease_{keyword}_+{score}")
            else:
                extracted_info.append(f"Lit_disease_+{score}")
        elif "NCBI核心關鍵字" in reason_str or "NCBI_core" in reason_str:
            keyword_match = re.search(r'NCBI核心關鍵字:\s*([^(]+)', reason_str)
            if keyword_match:
                keyword = keyword_match.group(1).strip()
                keyword = re.sub(r'\\b|\\\\b', '', keyword)
                keyword = clean_text_for_tsv(keyword)
                extracted_info.append(f"NCBI_core_{keyword}_+{score}")
            else:
                extracted_info.append(f"NCBI_core_+{score}")
        elif "NCBI通路關鍵字" in reason_str or "NCBI_pathway" in reason_str:
            keyword_match = re.search(r'NCBI通路關鍵字:\s*([^(]+)', reason_str)
            if keyword_match:
                keyword = keyword_match.group(1).strip()
                keyword = re.sub(r'\\b|\\\\b', '', keyword)
                keyword = clean_text_for_tsv(keyword)
                extracted_info.append(f"NCBI_pathway_{keyword}_+{score}")
            else:
                extracted_info.append(f"NCBI_pathway_+{score}")
        elif "NCBI機制關鍵字" in reason_str or "NCBI_mechanism" in reason_str:
            keyword_match = re.search(r'NCBI機制關鍵字:\s*([^(]+)', reason_str)
            if keyword_match:
                keyword = keyword_match.group(1).strip()
                keyword = clean_text_for_tsv(keyword)
                extracted_info.append(f"NCBI_mechanism_{keyword}_+{score}")
            else:
                extracted_info.append(f"NCBI_mechanism_+{score}")
        elif "NCBI基因名稱相關" in reason_str or "NCBI_gene_name" in reason_str:
            extracted_info.append(f"NCBI_gene_name_+{score}")
        elif "多文獻提及獎勵" in reason_str or "Multi_literature" in reason_str:
            extracted_info.append(f"Multi_lit_bonus_+{score}")
        elif "文獻時間上下文" in reason_str or "Time_context" in reason_str:
            keyword_match = re.search(r'文獻時間上下文:\s*([^(]+)', reason_str)
            if keyword_match:
                keyword = keyword_match.group(1).strip()
                keyword = re.sub(r'\\b|\\\\b', '', keyword)
                keyword = clean_text_for_tsv(keyword)
                extracted_info.append(f"Time_context_{keyword}_+{score}")
            else:
                extracted_info.append(f"Time_context_+{score}")
        else:
            # Generic handling for other factors
            cleaned_reason = clean_text_for_tsv(reason_str)
            if cleaned_reason and len(cleaned_reason) > 5:
                if "文獻" in reason_str:
                    extracted_info.append(f"Lit_other_factor_+{score}")
                elif "NCBI" in reason_str:
                    extracted_info.append(f"NCBI_other_factor_+{score}")
                else:
                    extracted_info.append(f"Other_factor_+{score}")
    
    return " | ".join(extracted_info) if extracted_info else "No_scoring_info"

def translate_status_to_english(ncbi_status):
    """Translate NCBI status to English"""
    if not ncbi_status:
        return "Unknown"
    
    status_mapping = {
        '成功': 'Success',
        '跳過NCBI查詢': 'Skip_NCBI_Query',
        '找不到基因': 'Gene_Not_Found',
        '查詢失敗': 'Query_Failed',
        '獲取失敗': 'Fetch_Failed',
        '達到最大重試次數': 'Max_Retries_Reached'
    }
    
    # Check for exact matches first
    for chinese, english in status_mapping.items():
        if chinese in ncbi_status:
            return english
    
    # Clean any remaining Chinese characters
    return clean_text_for_tsv(ncbi_status)

def convert_cache_to_tsv(cache_data, output_file, min_score=0, confidence_filter=None, circadian_only=False, truncate_summary=False):
    """Convert cache data to TSV format with English-only content"""
    
    if not cache_data:
        print("No data to convert")
        return
    
    # Filter records based on criteria
    filtered_results = []
    
    for gene_key, gene_data in cache_data.items():
        # Check if it's dual analysis version data
        if gene_data.get('analysis_version') not in ['dual_v1', 'dual_v1_optimized']:
            continue
            
        # Apply circadian-only filter
        if circadian_only and not gene_data.get('circadian_related', False):
            continue
            
        # Apply score filter
        combined_score = gene_data.get('combined_score', 0)
        if combined_score < min_score:
            continue
            
        # Apply confidence filter
        confidence = gene_data.get('confidence_level', 'VERY_LOW')
        if confidence_filter:
            confidence_order = {'HIGH': 4, 'MEDIUM': 3, 'LOW': 2, 'VERY_LOW': 1}
            min_confidence_level = confidence_order.get(confidence_filter, 1)
            current_confidence_level = confidence_order.get(confidence, 1)
            
            if current_confidence_level < min_confidence_level:
                continue
        
        filtered_results.append(gene_data)
    
    print(f"Filtered to {len(filtered_results)} genes meeting criteria")
    
    if not filtered_results:
        print("No genes meet the filtering criteria")
        return
    
    # Sort by combined score
    filtered_results.sort(key=lambda x: x.get('combined_score', 0), reverse=True)
    
    # Write TSV file
    try:
        with open(output_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
            
            # Write header - updated to match main program structure
            header = [
                'Gene_ID', 'Gene_Symbol', 'Gene_Type', 'Article_Count', 'Total_Occurrences',
                'Literature_Score', 'NCBI_Score', 'Combined_Score',
                'NCBI_Summary', 'Aliases', 'Literature_Contexts', 'Circadian_Keywords', 
                'Is_Circadian_Related', 'Confidence_Level', 'Score_Reasons', 'NCBI_Status'
            ]
            writer.writerow(header)
            
            # Write data
            for gene_data in filtered_results:
                # Validate record format
                if not isinstance(gene_data, dict) or 'gene_id' not in gene_data or 'analysis_version' not in gene_data:
                    print(f"Skipping invalid record: {type(gene_data)} - {str(gene_data)[:100]}")
                    continue
                
                # Process NCBI status
                ncbi_status = translate_status_to_english(gene_data.get('ncbi_status', 'Unknown'))
                
                # Process score reasons
                score_reasons = gene_data.get('score_reasons', [])
                reasons_text = extract_score_info_english(score_reasons)
                
                # Process literature contexts
                contexts = gene_data.get('literature_contexts', [])
                if contexts:
                    cleaned_contexts = []
                    for ctx in contexts[:3]:  # Take only first 3 contexts
                        if ctx:
                            cleaned_ctx = clean_text_for_tsv(ctx, 150)  # Limit to 150 chars
                            if cleaned_ctx:
                                cleaned_contexts.append(cleaned_ctx)
                    contexts_text = " | ".join(cleaned_contexts) if cleaned_contexts else "No_context"
                else:
                    contexts_text = "No_context"
                
                # 🎯 Process NCBI summary - default to full summary, optional truncation
                summary_max_length = 400 if truncate_summary else None
                ncbi_summary = clean_text_for_tsv(gene_data.get('ncbi_summary', ''), summary_max_length)
                if not ncbi_summary:
                    ncbi_summary = "No_summary"
                
                # Process circadian keywords
                circadian_keywords = gene_data.get('circadian_keywords', [])
                cleaned_keywords = []
                for kw in circadian_keywords:
                    cleaned_kw = clean_text_for_tsv(kw)
                    if cleaned_kw:
                        cleaned_keywords.append(cleaned_kw)
                keywords_text = ";".join(cleaned_keywords) if cleaned_keywords else "No_keywords"
                
                # Process aliases
                aliases = gene_data.get('aliases', [])
                cleaned_aliases = []
                for alias in aliases:
                    cleaned_alias = clean_text_for_tsv(alias)
                    if cleaned_alias:
                        cleaned_aliases.append(cleaned_alias)
                aliases_text = ";".join(cleaned_aliases) if cleaned_aliases else "No_aliases"
                
                # Prepare data row
                row = [
                    clean_text_for_tsv(gene_data.get('gene_id', '')),
                    clean_text_for_tsv(gene_data.get('gene_symbol', '')),
                    clean_text_for_tsv(gene_data.get('gene_type', 'Unknown')),
                    gene_data.get('article_count', 0),
                    gene_data.get('total_occurrences', 0),
                    gene_data.get('literature_score', 0),
                    gene_data.get('ncbi_score', 0),
                    gene_data.get('combined_score', 0),
                    ncbi_summary,  # 🎯 Full summary by default
                    aliases_text,
                    contexts_text,
                    keywords_text,
                    'Yes' if gene_data.get('circadian_related', False) else 'No',
                    gene_data.get('confidence_level', 'VERY_LOW'),
                    reasons_text,
                    ncbi_status
                ]
                
                writer.writerow(row)
        
        print(f"Successfully saved results to {output_file}")
        print(f"Successfully wrote {len(filtered_results)} gene records")
        if truncate_summary:
            print("⚠️ NCBI summary truncated to 400 characters")
        else:
            print("✅ Using full NCBI summary (no truncation)")
        print("Note: All non-ASCII characters removed, English content only")
        
    except Exception as e:
        print(f"Error: Problem saving TSV file - {e}")
        import traceback
        traceback.print_exc()

def generate_summary_stats(cache_data):
    """Generate summary statistics"""
    if not cache_data:
        return
    
    dual_analysis_data = [data for data in cache_data.values() 
                         if data.get('analysis_version') in ['dual_v1', 'dual_v1_optimized']]
    
    if not dual_analysis_data:
        print("No dual analysis version data found in cache")
        return
    
    print(f"\n{'='*60}")
    print(f"Cache File Summary Statistics")
    print(f"{'='*60}")
    
    total_genes = len(dual_analysis_data)
    circadian_related = [data for data in dual_analysis_data if data.get('circadian_related', False)]
    
    print(f"Total genes: {total_genes:,}")
    print(f"Circadian-related genes: {len(circadian_related):,} ({len(circadian_related)/total_genes*100:.1f}%)")
    
    # 🎯 New: Circadian_Keywords statistics analysis
    no_keywords_count = 0
    has_keywords_count = 0
    
    for data in dual_analysis_data:
        circadian_keywords = data.get('circadian_keywords', [])
        if not circadian_keywords or circadian_keywords == []:
            no_keywords_count += 1
        else:
            has_keywords_count += 1
    
    print(f"\nCircadian_Keywords field statistics:")
    print(f"  No_keywords: {no_keywords_count:,} genes ({no_keywords_count/total_genes*100:.1f}%)")
    print(f"  Has_keywords: {has_keywords_count:,} genes ({has_keywords_count/total_genes*100:.1f}%)")
    if no_keywords_count > total_genes * 0.8:
        print(f"  ⚠️ Over 80% genes show No_keywords, recommend updating main program keyword extraction logic")
    
    # Version distribution
    version_counts = {}
    for data in dual_analysis_data:
        version = data.get('analysis_version', 'Unknown')
        version_counts[version] = version_counts.get(version, 0) + 1
    
    print(f"\nVersion distribution:")
    for version, count in version_counts.items():
        print(f"  {version}: {count:,} genes")
    
    # Confidence distribution
    confidence_counts = {}
    for data in circadian_related:
        conf = data.get('confidence_level', 'VERY_LOW')
        confidence_counts[conf] = confidence_counts.get(conf, 0) + 1
    
    print(f"\nConfidence distribution (circadian-related genes):")
    for conf in ['HIGH', 'MEDIUM', 'LOW', 'VERY_LOW']:
        count = confidence_counts.get(conf, 0)
        if count > 0:
            percentage = count / len(circadian_related) * 100
            print(f"  {conf}: {count:,} genes ({percentage:.1f}%)")
    
    # Score distribution
    if circadian_related:
        lit_scores = [data.get('literature_score', 0) for data in circadian_related]
        ncbi_scores = [data.get('ncbi_score', 0) for data in circadian_related]
        combined_scores = [data.get('combined_score', 0) for data in circadian_related]
        
        print(f"\nScore statistics (circadian-related genes):")
        print(f"  Literature scores - Max: {max(lit_scores):.1f}, Avg: {sum(lit_scores)/len(lit_scores):.1f}, Min: {min(lit_scores):.1f}")
        print(f"  NCBI scores - Max: {max(ncbi_scores):.1f}, Avg: {sum(ncbi_scores)/len(ncbi_scores):.1f}, Min: {min(ncbi_scores):.1f}")
        print(f"  Combined scores - Max: {max(combined_scores):.1f}, Avg: {sum(combined_scores)/len(combined_scores):.1f}, Min: {min(combined_scores):.1f}")
    
    # Top 15 genes (consistent with main program report)
    top_genes = sorted(circadian_related, key=lambda x: x.get('combined_score', 0), reverse=True)[:15]
    print(f"\nTop 15 circadian-related genes (by combined score):")
    for i, gene in enumerate(top_genes, 1):
        symbol = gene.get('gene_symbol', 'Unknown')
        lit_score = gene.get('literature_score', 0)
        ncbi_score = gene.get('ncbi_score', 0)
        combined_score = gene.get('combined_score', 0)
        confidence = gene.get('confidence_level', 'Unknown')
        article_count = gene.get('article_count', 0)
        version_mark = "✓" if gene.get('analysis_version') == 'dual_v1_optimized' else "○"
        
        # 🎯 New: Show circadian_keywords status
        circadian_keywords = gene.get('circadian_keywords', [])
        keywords_status = f"Keywords:{len(circadian_keywords)}" if circadian_keywords else "No_keywords"
        
        print(f"  {i:2d}. {symbol:<12} (Lit:{lit_score:4.1f} + NCBI:{ncbi_score:4.1f} = {combined_score:4.1f}, {confidence:<6}, Articles:{article_count:4d}, {keywords_status}) {version_mark}")
    
    # NCBI query statistics (enhanced)
    ncbi_attempted = [data for data in dual_analysis_data if data.get('ncbi_status') not in ['跳過NCBI查詢', 'Skip_NCBI_Query']]
    ncbi_success = [data for data in dual_analysis_data if data.get('ncbi_status') == '成功']
    ncbi_failed = [data for data in dual_analysis_data if data.get('ncbi_status', '').startswith('查詢失敗')]
    ncbi_not_found = [data for data in dual_analysis_data if data.get('ncbi_status') == '找不到基因']
    ncbi_no_summary = [data for data in dual_analysis_data if data.get('ncbi_status') == '成功' and not data.get('ncbi_summary')]
    
    print(f"\nNCBI query statistics:")
    print(f"  Attempted queries: {len(ncbi_attempted):,} genes")
    print(f"  Successful queries: {len(ncbi_success):,} genes")
    print(f"  Failed queries: {len(ncbi_failed):,} genes")
    print(f"  Genes not found: {len(ncbi_not_found):,} genes")
    print(f"  Success but no summary: {len(ncbi_no_summary):,} genes")
    if len(ncbi_attempted) > 0:
        print(f"  Query success rate: {len(ncbi_success)/len(ncbi_attempted)*100:.2f}%")
    
    # Gene type distribution (new)
    gene_types = {}
    for data in dual_analysis_data:
        gene_type = data.get('gene_type', 'Unknown')
        gene_types[gene_type] = gene_types.get(gene_type, 0) + 1
    
    print(f"\nGene type distribution (Top 10):")
    sorted_types = sorted(gene_types.items(), key=lambda x: x[1], reverse=True)[:10]
    for gene_type, count in sorted_types:
        percentage = count / total_genes * 100
        print(f"  {gene_type:<15}: {count:,} genes ({percentage:.1f}%)")
    
    print(f"{'='*60}\n")

def main():
    parser = argparse.ArgumentParser(description='Convert cache file to TSV results - English Version (Updated)')
    parser.add_argument('--cache-file', type=str, 
                       default='circadian_genes_cache_dual.json',
                       help='Cache file path')
    parser.add_argument('--output', type=str, 
                       default=f'circadian_genes_cache_results_english_{datetime.now().strftime("%Y%m%d_%H%M%S")}.tsv',
                       help='Output TSV file path')
    parser.add_argument('--min-score', type=float, default=0,
                       help='Minimum combined score filter (default: 0)')
    parser.add_argument('--min-confidence', type=str, 
                       choices=['HIGH', 'MEDIUM', 'LOW', 'VERY_LOW'],
                       help='Minimum confidence level filter')
    parser.add_argument('--circadian-only', action='store_true',
                       help='Output only circadian-related genes')
    parser.add_argument('--show-stats', action='store_true',
                       help='Display summary statistics')
    parser.add_argument('--high-confidence', action='store_true',
                       help='Quick option: output only high-confidence circadian genes (equivalent to --circadian-only --min-confidence HIGH)')
    # 🎯 New parameter: Control NCBI summary truncation (default to full output)
    parser.add_argument('--truncate-summary', action='store_true',
                       help='Truncate NCBI summary to 400 characters (default is full output)')
    
    args = parser.parse_args()
    
    # Handle quick option
    if args.high_confidence:
        args.circadian_only = True
        args.min_confidence = 'HIGH'
        print("Using high-confidence quick option: only high-confidence circadian-related genes")
    
    # Check if cache file exists
    if not os.path.exists(args.cache_file):
        print(f"Error: Cache file does not exist - {args.cache_file}")
        print("Please confirm the cache file path is correct")
        sys.exit(1)
    
    # Load cache file
    print(f"Loading cache file: {args.cache_file}")
    cache_data = load_cache_file(args.cache_file)
    
    if cache_data is None:
        sys.exit(1)
    
    # Display summary statistics
    if args.show_stats:
        generate_summary_stats(cache_data)
    
    # Convert to TSV
    print(f"Converting to TSV format...")
    print(f"Filter criteria:")
    print(f"  - Minimum score: {args.min_score}")
    if args.min_confidence:
        print(f"  - Minimum confidence: {args.min_confidence}")
    if args.circadian_only:
        print(f"  - Include only circadian-related genes")
    if args.truncate_summary:
        print(f"  - NCBI summary truncated to 400 characters")
    else:
        print(f"  - Using full NCBI summary (no truncation)")
    
    convert_cache_to_tsv(cache_data, args.output, args.min_score, args.min_confidence, args.circadian_only, args.truncate_summary)
    
    print(f"Conversion complete! File saved at: {args.output}")

if __name__ == "__main__":
    try:
        print("\n" + "="*60)
        print("  Cache to TSV Tool - English Version (Updated)")
        print("  - Extract results from circadian_genes_cache_dual.json")
        print("  - Convert to TSV format with English-only content")
        print("  - Default: Full NCBI summary output (no truncation)")
        print("  - Enhanced Circadian_Keywords statistics analysis")
        print("  - Maximum compatibility with Excel and analysis tools")
        print("="*60 + "\n")
        
        main()
        
    except KeyboardInterrupt:
        print("\nProgram interrupted by user")
    except Exception as e:
        print(f"Unexpected error occurred: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)