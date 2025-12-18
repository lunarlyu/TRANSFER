#!/usr/bin/env python3
"""
Compare Results Script

Compares output files between two result directories to identify differences.
Useful for verifying reproducibility of the pipeline.

Usage:
    python compare_results.py <dir1> <dir2>
    python compare_results.py Results Check_data  # Default paths

Example:
    python compare_results.py Results_Orig Check_data/intermediate_results_for_check
"""

import argparse
import csv
import ast
import os
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional

try:
    import pandas as pd
except ImportError:
    print("Error: pandas is required. Install with: pip install pandas")
    sys.exit(1)


def load_labels(filepath: str) -> Dict[Tuple[str, str], Set[str]]:
    """Load labels file into a dictionary."""
    df = pd.read_csv(filepath, sep='\t')
    result = {}
    for _, row in df.iterrows():
        tissue = row['Tissue']
        cell_type = row['Cell type']
        # Handle both column name variants
        gene_col = [c for c in df.columns if 'Gene Names' in c][0]
        genes = set(ast.literal_eval(row[gene_col]))
        result[(tissue, cell_type)] = genes
    return result


def load_recommendations(filepath: str) -> Dict[str, List[str]]:
    """Load recommendations file into a dictionary."""
    result = {}
    with open(filepath, 'r') as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        for row in reader:
            if len(row) >= 2:
                result[row[0]] = ast.literal_eval(row[1])
    return result


def compare_text_files(file1: str, file2: str) -> bool:
    """Compare two text files for exact equality."""
    with open(file1) as f1, open(file2) as f2:
        return f1.read() == f2.read()


def compare_expression_matrix(file1: str, file2: str, tolerance: float = 0.001) -> Dict:
    """Compare expression matrices numerically."""
    df1 = pd.read_csv(file1, sep='\t')
    df2 = pd.read_csv(file2, sep='\t')
    
    format_diffs = 0
    numeric_diffs = 0
    diff_details = []
    
    for col in df1.columns:
        if col in ['Tissue', 'Cell type']:
            continue
        if col not in df2.columns:
            continue
        for i in range(min(len(df1), len(df2))):
            v1 = df1[col].iloc[i]
            v2 = df2[col].iloc[i]
            if v1 != v2:
                try:
                    diff = abs(float(v1) - float(v2))
                    if diff < tolerance:
                        format_diffs += 1
                    else:
                        numeric_diffs += 1
                        if len(diff_details) < 5:
                            diff_details.append({
                                'row': i,
                                'gene': col,
                                'val1': v1,
                                'val2': v2,
                                'diff': diff
                            })
                except:
                    numeric_diffs += 1
    
    return {
        'format_diffs': format_diffs,
        'numeric_diffs': numeric_diffs,
        'details': diff_details
    }


def compare_labels(file1: str, file2: str) -> Dict:
    """Compare label files."""
    labels1 = load_labels(file1)
    labels2 = load_labels(file2)
    
    only_in_1 = set(labels1.keys()) - set(labels2.keys())
    only_in_2 = set(labels2.keys()) - set(labels1.keys())
    common = set(labels1.keys()) & set(labels2.keys())
    
    same_genes = 0
    diff_genes = []
    
    for cell in common:
        if labels1[cell] == labels2[cell]:
            same_genes += 1
        else:
            diff_genes.append({
                'cell': cell,
                'only_in_1': labels1[cell] - labels2[cell],
                'only_in_2': labels2[cell] - labels1[cell]
            })
    
    return {
        'rows1': len(labels1),
        'rows2': len(labels2),
        'only_in_1': list(only_in_1),
        'only_in_2': list(only_in_2),
        'common': len(common),
        'same_genes': same_genes,
        'diff_genes': diff_genes
    }


def compare_recommendations(file1: str, file2: str) -> Dict:
    """Compare recommendation files."""
    rec1 = load_recommendations(file1)
    rec2 = load_recommendations(file2)
    
    exact_match = 0
    same_set = 0
    different = []
    
    all_cells = set(rec1.keys()) | set(rec2.keys())
    
    for cell in all_cells:
        markers1 = rec1.get(cell, [])
        markers2 = rec2.get(cell, [])
        
        if markers1 == markers2:
            exact_match += 1
            same_set += 1
        elif set(markers1) == set(markers2):
            same_set += 1
        else:
            if len(different) < 10:
                different.append({
                    'cell': cell,
                    'markers1': markers1[:5],
                    'markers2': markers2[:5],
                    'only_in_1': set(markers1) - set(markers2),
                    'only_in_2': set(markers2) - set(markers1)
                })
    
    return {
        'total': len(all_cells),
        'exact_match': exact_match,
        'same_set': same_set,
        'different': different
    }


def find_file(directory: str, patterns: List[str]) -> Optional[str]:
    """Find a file matching one of the patterns in the directory tree."""
    dir_path = Path(directory)
    for pattern in patterns:
        # Direct match
        if (dir_path / pattern).exists():
            return str(dir_path / pattern)
        # Search subdirectories
        matches = list(dir_path.rglob(pattern))
        if matches:
            return str(matches[0])
        # Try with -2 suffix (Check_data format)
        pattern_2 = pattern.replace('.csv', '-2.csv').replace('.tsv', '-2.tsv')
        matches = list(dir_path.rglob(pattern_2))
        if matches:
            return str(matches[0])
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Compare result directories",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('dir1', nargs='?', default='Results',
                        help='First directory to compare (default: Results)')
    parser.add_argument('dir2', nargs='?', default='Check_data',
                        help='Second directory to compare (default: Check_data)')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Show detailed differences')
    
    args = parser.parse_args()
    
    print("=" * 70)
    print(f"COMPARING: {args.dir1} vs {args.dir2}")
    print("=" * 70)
    
    # Define files to compare
    files_to_compare = [
        {
            'name': 'tissue_cell_pairs.tsv',
            'type': 'text',
            'patterns': ['tissue_cell_pairs.tsv', 'tissue_cell_pairs-2.tsv']
        },
        {
            'name': 'gene_list.csv',
            'type': 'text',
            'patterns': ['gene_list.csv']
        },
        {
            'name': 'gene_expression_matrix_high.csv',
            'type': 'expression',
            'patterns': ['gene_expression_matrix_high.csv', 'gene_expression_matrix_high-2.csv']
        },
        {
            'name': 'gene_expression_matrix_median.csv',
            'type': 'expression',
            'patterns': ['gene_expression_matrix_median.csv', 'gene_expression_matrix_median-2.csv']
        },
        {
            'name': 'positives_labels.csv',
            'type': 'labels',
            'patterns': ['positives_labels.csv', 'positives_labels-2.csv']
        },
        {
            'name': 'negative_labels.csv',
            'type': 'labels',
            'patterns': ['negative_labels.csv', 'negative_labels-2.csv']
        },
        {
            'name': 'recommended_whole_body_markers_high.csv',
            'type': 'recommendations',
            'patterns': ['recommended_whole_body_markers_high.csv', 'recommended_whole_body_markers_high-2.csv']
        },
        {
            'name': 'recommended_whole_body_markers_median.csv',
            'type': 'recommendations',
            'patterns': ['recommended_whole_body_markers_median.csv', 'recommended_whole_body_markers_median-2.csv']
        }
    ]
    
    summary = {'identical': [], 'different': [], 'missing': []}
    
    for file_info in files_to_compare:
        print(f"\n{'-' * 50}")
        print(f"{file_info['name']}")
        print(f"{'-' * 50}")
        
        file1 = find_file(args.dir1, file_info['patterns'])
        file2 = find_file(args.dir2, file_info['patterns'])
        
        if not file1:
            print(f"  NOT FOUND in {args.dir1}")
            summary['missing'].append(f"{file_info['name']} (dir1)")
            continue
        if not file2:
            print(f"  NOT FOUND in {args.dir2}")
            summary['missing'].append(f"{file_info['name']} (dir2)")
            continue
        
        print(f"  File 1: {file1}")
        print(f"  File 2: {file2}")
        
        try:
            if file_info['type'] == 'text':
                identical = compare_text_files(file1, file2)
                print(f"  Identical: {identical}")
                if identical:
                    summary['identical'].append(file_info['name'])
                else:
                    summary['different'].append(file_info['name'])
            
            elif file_info['type'] == 'expression':
                result = compare_expression_matrix(file1, file2)
                print(f"  Format-only diffs: {result['format_diffs']}")
                print(f"  Numeric diffs: {result['numeric_diffs']}")
                if result['numeric_diffs'] == 0:
                    summary['identical'].append(file_info['name'])
                else:
                    summary['different'].append(file_info['name'])
                    if args.verbose and result['details']:
                        print("  Sample differences:")
                        for d in result['details'][:3]:
                            print(f"    Row {d['row']}, {d['gene']}: {d['val1']} vs {d['val2']}")
            
            elif file_info['type'] == 'labels':
                result = compare_labels(file1, file2)
                print(f"  Rows: {result['rows1']} vs {result['rows2']}")
                print(f"  Only in dir1: {len(result['only_in_1'])}")
                print(f"  Only in dir2: {len(result['only_in_2'])}")
                print(f"  Common cells: {result['common']}")
                print(f"  Same genes: {result['same_genes']}/{result['common']}")
                print(f"  Different genes: {len(result['diff_genes'])}")
                
                if result['only_in_1'] == [] and result['only_in_2'] == [] and len(result['diff_genes']) == 0:
                    summary['identical'].append(file_info['name'])
                else:
                    summary['different'].append(file_info['name'])
                    if args.verbose:
                        if result['only_in_2']:
                            print(f"  Missing from dir1: {result['only_in_2'][:5]}")
                        for d in result['diff_genes'][:3]:
                            print(f"  {d['cell']}: +{d['only_in_1']}, -{d['only_in_2']}")
            
            elif file_info['type'] == 'recommendations':
                result = compare_recommendations(file1, file2)
                pct_exact = 100 * result['exact_match'] / result['total'] if result['total'] > 0 else 0
                pct_set = 100 * result['same_set'] / result['total'] if result['total'] > 0 else 0
                print(f"  Total cells: {result['total']}")
                print(f"  Exact match: {result['exact_match']} ({pct_exact:.1f}%)")
                print(f"  Same markers (any order): {result['same_set']} ({pct_set:.1f}%)")
                
                if result['exact_match'] == result['total']:
                    summary['identical'].append(file_info['name'])
                else:
                    summary['different'].append(file_info['name'])
                    if args.verbose:
                        for d in result['different'][:3]:
                            print(f"  {d['cell']}:")
                            print(f"    dir1: {d['markers1']}")
                            print(f"    dir2: {d['markers2']}")
        
        except Exception as e:
            print(f"  ERROR: {e}")
            summary['different'].append(f"{file_info['name']} (error)")
    
    # Print summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"\n✓ Identical ({len(summary['identical'])}):")
    for f in summary['identical']:
        print(f"    {f}")
    
    print(f"\n✗ Different ({len(summary['different'])}):")
    for f in summary['different']:
        print(f"    {f}")
    
    if summary['missing']:
        print(f"\n? Missing ({len(summary['missing'])}):")
        for f in summary['missing']:
            print(f"    {f}")
    
    print()
    return 0 if not summary['different'] else 1


if __name__ == "__main__":
    sys.exit(main())

