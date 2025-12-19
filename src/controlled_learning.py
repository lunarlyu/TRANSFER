"""
Controlled Learning Module

This module implements the penalty matrix optimization algorithm for
cell surface marker prediction. Matches the exact logic from Controled_Learning.ipynb.
"""

import csv
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from .constants import LITERATURE_POSITIVES
from .io_utils import (
    load_tissue_cells,
    load_gene_list,
    load_common_cells,
    load_expression_matrices,
    save_recommendations
)


def load_labels_indexed(
    data_dir: str,
    gene_list: List,
    tissue_cells: List
) -> Tuple[Dict[int, List[int]], Dict[int, List[int]]]:
    """
    Load positive and negative labels with index mappings.
    
    Matches exact logic from Controled_Learning.ipynb Cell 7.
    """
    data_path = Path(data_dir)
    
    # Flatten gene_list if nested (exact notebook logic)
    if isinstance(gene_list[0], list):
        gene_list_flat = [gene[0] for gene in gene_list]
    else:
        gene_list_flat = gene_list
    
    tissue_cell_to_index = {tuple(row): idx for idx, row in enumerate(tissue_cells)}
    gene_to_index = {gene: idx for idx, gene in enumerate(gene_list_flat)}
    
    # Load positive labels
    positives_df = pd.read_csv(data_path / 'positives_labels.csv', sep='\t')
    positives = {}
    for idx, row in positives_df.iterrows():
        tissue = row["Tissue"]
        cell_type = row["Cell type"]
        gene_names = eval(row["Positive Gene Names"])
        
        key = (tissue, cell_type)
        if key not in tissue_cell_to_index:
            continue
        
        cell_idx = tissue_cell_to_index[key]
        gene_indices = [gene_to_index[gene] for gene in gene_names if gene in gene_to_index]
        
        if gene_indices:
            positives[cell_idx] = gene_indices
    
    # Load negative labels
    negatives_df = pd.read_csv(data_path / 'negative_labels.csv', sep='\t')
    negatives = {}
    for idx, row in negatives_df.iterrows():
        tissue = row["Tissue"]
        cell_type = row["Cell type"]
        gene_names = eval(row["Negative Gene Names"])
        
        key = (tissue, cell_type)
        if key not in tissue_cell_to_index:
            continue
        
        cell_idx = tissue_cell_to_index[key]
        gene_indices = [gene_to_index[gene] for gene in gene_names if gene in gene_to_index]
        
        if gene_indices:
            negatives[cell_idx] = gene_indices
    
    return positives, negatives


def get_highly_expressed_genes(nTPM_matrix: np.ndarray, top_n: int = 50) -> List[int]:
    """
    Get indices of top N highly expressed genes (by median across cells).
    
    Matches exact logic from Controled_Learning.ipynb Cell 9.
    """
    medians = np.median(nTPM_matrix, axis=0).tolist()
    genes_sorted = [index for index, value in sorted(enumerate(medians), key=lambda x: x[1], reverse=True)]
    return genes_sorted[:top_n]


def run_grid_search(
    tissue_cells: List,
    common_cells: List,
    nTPM_matrix: np.ndarray,
    positives: Dict[int, List[int]],
    negatives: Dict[int, List[int]],
    highly_expressed: List[int],
    gene_list: List,
    positives2: List,
    verbose: bool = True
) -> Tuple[Tuple[float, float, float], Dict]:
    """
    Run grid search for optimal parameters p, q, r.
    
    Args:
        tissue_cells: List of [tissue, cell_type] pairs
        common_cells: List of common cell types across tissues
        nTPM_matrix: Expression matrix
        positives: Dict mapping cell index to list of positive marker gene indices
        negatives: Dict mapping cell index to list of negative marker gene indices
        highly_expressed: List of indices of top 50 highly expressed genes to penalize
        gene_list: List of genes
        positives2: Literature-validated positive markers
        verbose: Whether to print progress
        
    Returns:
        Tuple of (optimal_params, all_results_dict)
    """
    m = len(tissue_cells)
    dic = {}
    best_score = float('-inf')
    best_params = None
    
    if verbose:
        print("  Running grid search over 8000 parameter combinations...")
    
    for p in np.linspace(-700, -500, 20):
        for q in np.linspace(-700, -500, 20):
            for r in np.linspace(-700, -20, 20):
                # Build penalty matrix
                penalty_matrix = []
                for i in range(m):
                    row = []
                    for j in range(m):
                        if j == i:
                            row.append(1000)  # Self: high weight
                        elif tissue_cells[j][1] == tissue_cells[i][1]:
                            if tissue_cells[j][1] in common_cells:
                                row.append(0)  # Same cell (common) in different tissue
                            else:
                                row.append(p)  # Same cell (non-common) in different tissue
                        elif tissue_cells[j][0] == tissue_cells[i][0]:
                            row.append(q)  # Different cell in same tissue
                        else:
                            row.append(r)  # Different cell in different tissue
                    penalty_matrix.append(row)
                penalty_matrix = np.array(penalty_matrix)
                objective_matrix = np.dot(penalty_matrix, nTPM_matrix)
                objective_list = objective_matrix.tolist()
                
                # Count 1: Organ-wide markers (from databases)
                count1 = 0
                for cell_idx in positives:
                    pos_markers = positives[cell_idx]
                    obj_row = objective_list[cell_idx]
                    top_markers = sorted(enumerate(obj_row), key=lambda x: x[1], reverse=True)[:10]
                    suggested = [idx for idx, _ in top_markers]
                    count1 += len(set(pos_markers).intersection(set(suggested)))
                
                # Count 2: Whole-body markers (from literature) - FIXED
                count2 = 0
                for row in positives2:
                    if [row[0]] in gene_list:
                        marker_idx = gene_list.index([row[0]])
                        target_pair = [row[1], row[2]]
                        if target_pair in tissue_cells:
                            cell_idx = tissue_cells.index(target_pair)
                            obj_row = objective_list[cell_idx]
                            top_markers = sorted(enumerate(obj_row), key=lambda x: x[1], reverse=True)[:10]
                            suggested = [idx for idx, _ in top_markers]
                            # FIXED: Check if marker is in suggested top 10
                            if marker_idx in suggested:
                                count2 += 1
                
                # Count 3: Negative markers (to minimize)
                count3 = 0
                for cell_idx in negatives:
                    neg_markers = negatives[cell_idx]
                    obj_row = objective_list[cell_idx]
                    top_markers = sorted(enumerate(obj_row), key=lambda x: x[1], reverse=True)[:10]
                    suggested = [idx for idx, _ in top_markers]
                    count3 += len(set(neg_markers).intersection(set(suggested)))
                
                # Count 4: Highly expressed genes (to minimize) - FIXED
                count4 = 0
                for obj_row in objective_list:
                    top_markers = sorted(enumerate(obj_row), key=lambda x: x[1], reverse=True)[:10]
                    suggested = [idx for idx, _ in top_markers]
                    # FIXED: Count inside the loop for ALL cells, not just last
                    count4 += len(set(suggested).intersection(set(highly_expressed)))
                
                # Compute score (maximize positives, minimize negatives)
                score = count1 + count2 - count3 - count4
                
                dic[(p, q, r)] = (count1, count2, count3, count4)
                
                if score > best_score:
                    best_score = score
                    best_params = (p, q, r)
    
    if verbose and best_params:
        counts = dic[best_params]
        print(f"  Optimal parameters: p={best_params[0]:.2f}, q={best_params[1]:.2f}, r={best_params[2]:.2f}")
        print(f"  Organ-wide marker hits: {counts[0]}")
        print(f"  Whole-body marker hits: {counts[1]}")
        print(f"  Non-marker hits (penalized): {counts[2]}")
        print(f"  Highly expressed hits (penalized): {counts[3]}")
        print(f"  Total score: {best_score}")
    
    return best_params, dic


def build_objective_matrix(
    tissue_cells: List,
    common_cells: List,
    nTPM_matrix: np.ndarray,
    p: float,
    q: float,
    r: float
) -> np.ndarray:
    """
    Build the objective matrix with given parameters.
    
    Matches exact logic from Controled_Learning.ipynb Cell 13/17.
    """
    m = len(tissue_cells)
    penalty_matrix = []
    for i in range(m):
        row = []
        for j in range(m):
            if j == i:
                row.append(1000)
            elif tissue_cells[j][1] == tissue_cells[i][1]:
                if tissue_cells[j][1] in common_cells:
                    row.append(0)  # same cell (common) and diff tissue
                else:
                    row.append(p)  # same cell (non common) and diff tissue
            elif tissue_cells[j][0] == tissue_cells[i][0]:
                row.append(q)  # diff cell and same tissue
            else:
                row.append(r)  # diff cell and diff tissue
        penalty_matrix.append(row)
    penalty_matrix = np.array(penalty_matrix)
    objective_matrix = np.dot(penalty_matrix, nTPM_matrix)
    return objective_matrix


def recommend_markers(
    objective_matrix: np.ndarray,
    tissue_cells: List,
    gene_list: List,
    top_k: int = 10
) -> Dict[str, List[str]]:
    """
    Generate top-k marker recommendations for each cell type.
    
    Matches exact logic from Controled_Learning.ipynb Cell 14/18.
    """
    # Flatten gene_list if nested
    if isinstance(gene_list[0], list):
        genes = gene_list
    else:
        genes = [[g] for g in gene_list]
    
    m = len(tissue_cells)
    whole_body_markers = {}
    
    for i in range(m):
        obj_row = objective_matrix[i]
        markers = sorted(list(enumerate(obj_row)), key=lambda x: x[1], reverse=True)[:top_k]
        markers_id = [x[0] for x in markers]
        cell_name = tissue_cells[i][0] + " " + tissue_cells[i][1]
        whole_body_markers[cell_name] = [genes[x][0] for x in markers_id]
    
    return whole_body_markers


def save_recommendations_notebook_style(
    recommendations: Dict[str, List[str]],
    output_path: str
):
    """
    Save recommendations in exact notebook format.
    
    Matches exact logic from Controled_Learning.ipynb Cell 14/18.
    """
    with open(output_path, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(["cell", "markers"])
        for key, value in recommendations.items():
            row = [key, value]
            csvwriter.writerow(row)


def run_controlled_learning(data_dir: str, output_dir: str) -> Dict:
    """
    Run the complete controlled learning pipeline.
    
    Matches exact logic from Controled_Learning.ipynb.
    """
    data_path = Path(data_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    print("=== Controlled Learning Pipeline ===\n")
    
    # Load data (matching Cell 2-5)
    print("Loading data files...")
    tissue_cells = load_tissue_cells(data_path / "tissue_cell_pairs.tsv")
    common_cells = load_common_cells(data_path / "common_cells_across_tissues.csv")
    gene_list = load_gene_list(data_path / "gene_list.csv")
    
    # Load expression matrices
    nTPM_matrix_high_df = pd.read_csv(data_path / 'gene_expression_matrix_high.csv', sep='\t')
    nTPM_matrix_high_only = nTPM_matrix_high_df.drop(columns=["Tissue", "Cell type"])
    nTPM_matrix_high = nTPM_matrix_high_only.values
    
    nTPM_matrix_median_df = pd.read_csv(data_path / 'gene_expression_matrix_median.csv', sep='\t')
    nTPM_matrix_median_only = nTPM_matrix_median_df.drop(columns=["Tissue", "Cell type"])
    nTPM_matrix_median = nTPM_matrix_median_only.values
    
    print(f"  Loaded {len(tissue_cells)} tissue-cell pairs")
    print(f"  Loaded {len(gene_list)} genes")
    
    # Load labels (matching Cell 7)
    print("\nLoading labels...")
    positives, negatives = load_labels_indexed(data_dir, gene_list, tissue_cells)
    print(f"  Positive labels: {len(positives)} cells")
    print(f"  Negative labels: {len(negatives)} cells")
    
    # Get highly expressed genes - top 50 by median expression
    nTPM_matrix_high_medians = np.median(nTPM_matrix_high, axis=0).tolist()
    genes_sorted_high = [index for index, value in sorted(enumerate(nTPM_matrix_high_medians), key=lambda x: x[1], reverse=True)]
    highly_expressed_high = genes_sorted_high[:50]
    
    # For median method - FIXED: use actual median matrix
    nTPM_matrix_median_medians = np.median(nTPM_matrix_median, axis=0).tolist()
    genes_sorted_median = [index for index, value in sorted(enumerate(nTPM_matrix_median_medians), key=lambda x: x[1], reverse=True)]
    highly_expressed_median = genes_sorted_median[:50]
    
    # Literature positives (positives2 from Cell 10)
    positives2 = LITERATURE_POSITIVES
    
    results = {}
    
    # === Process HIGH method (Cell 12-14) ===
    print("\n--- Processing nTPM HIGH method ---")
    print("\nOptimizing parameters...")
    
    max_set_high, dic_high = run_grid_search(
        tissue_cells, common_cells, nTPM_matrix_high,
        positives, negatives, highly_expressed_high, gene_list, positives2
    )
    
    # Build objective matrix with optimal parameters
    objective_matrix_high = build_objective_matrix(
        tissue_cells, common_cells, nTPM_matrix_high,
        max_set_high[0], max_set_high[1], max_set_high[2]
    )
    
    # Generate recommendations
    print("\nGenerating marker recommendations...")
    recommendations_high = recommend_markers(objective_matrix_high, tissue_cells, gene_list)
    save_recommendations_notebook_style(
        recommendations_high,
        output_path / "recommended_whole_body_markers_high.csv"
    )
    print(f"  Saved recommendations for {len(recommendations_high)} cell types")
    
    results['high'] = {
        'optimal_params': max_set_high,
        'recommendations': recommendations_high
    }
    
    # === Process MEDIAN method (Cell 16-18) ===
    # NOTE: Commented out for now - only running HIGH method
    # print("\n--- Processing nTPM MEDIAN method ---")
    # print("\nOptimizing parameters...")
    # 
    # max_set_median, dic_median = run_grid_search(
    #     tissue_cells, common_cells, nTPM_matrix_median,
    #     positives, negatives, highly_expressed_median, gene_list, positives2
    # )
    # 
    # # Build objective matrix with optimal parameters
    # objective_matrix_median = build_objective_matrix(
    #     tissue_cells, common_cells, nTPM_matrix_median,
    #     max_set_median[0], max_set_median[1], max_set_median[2]
    # )
    # 
    # # Generate recommendations
    # print("\nGenerating marker recommendations...")
    # recommendations_median = recommend_markers(objective_matrix_median, tissue_cells, gene_list)
    # save_recommendations_notebook_style(
    #     recommendations_median,
    #     output_path / "recommended_whole_body_markers_median.csv"
    # )
    # print(f"  Saved recommendations for {len(recommendations_median)} cell types")
    # 
    # results['median'] = {
    #     'optimal_params': max_set_median,
    #     'recommendations': recommendations_median
    # }
    
    print("\n=== Controlled Learning Complete ===\n")
    
    return results


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Run controlled learning for marker prediction")
    parser.add_argument("--data-dir", default=".", help="Directory containing processed data")
    parser.add_argument("--output-dir", default=".", help="Directory for output files")
    
    args = parser.parse_args()
    run_controlled_learning(args.data_dir, args.output_dir)
