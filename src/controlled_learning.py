"""
Controlled Learning Module

This module implements the penalty matrix optimization algorithm for
cell surface marker prediction.

The algorithm optimizes three parameters (p, q, r) that control how expression
from different tissue-cell combinations is weighted when scoring marker candidates.
"""

import csv
from itertools import product
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

from .constants import LITERATURE_POSITIVES
from .io_utils import (
    load_tissue_cells,
    load_gene_list,
    load_common_cells,
)


# =============================================================================
# Penalty Matrix Construction
# =============================================================================

def precompute_penalty_masks(
    tissue_cells: List[List[str]],
    common_cells: List
) -> Dict[str, np.ndarray]:
    """
    Precompute boolean masks for penalty matrix construction.
    
    These masks are independent of p, q, r values and only need to be
    computed once per run.
    
    Args:
        tissue_cells: List of [tissue, cell_type] pairs
        common_cells: List of common cell types across tissues
        
    Returns:
        Dictionary of boolean masks for each penalty category
    """
    m = len(tissue_cells)
    
    # Extract tissue and cell arrays for vectorized comparison
    tissues = np.array([tc[0] for tc in tissue_cells])
    cells = np.array([tc[1] for tc in tissue_cells])
    
    # Basic boolean masks (m x m matrices)
    mask_self = np.eye(m, dtype=bool)
    mask_same_cell = cells[:, None] == cells[None, :]
    mask_same_tissue = tissues[:, None] == tissues[None, :]
    
    # Check if cell[j] is in common_cells
    # Note: common_cells may contain single-element lists like ['macrophages'] 
    # and strings like 'serous glandular cells'. The original notebook code
    # checks `cell in common_cells` directly, which only matches strings.
    # We preserve this exact behavior for compatibility.
    is_common = np.array([cell in common_cells for cell in cells])
    mask_j_common = np.broadcast_to(is_common[None, :], (m, m))
    
    # Combined masks for each penalty case
    return {
        'self': mask_self,
        'same_cell_common': mask_same_cell & mask_j_common & ~mask_self,
        'same_cell_noncommon': mask_same_cell & ~mask_j_common & ~mask_self,
        'diff_cell_same_tissue': ~mask_same_cell & mask_same_tissue,
        'diff_cell_diff_tissue': ~mask_same_cell & ~mask_same_tissue & ~mask_self,
    }


def build_penalty_matrix(
    masks: Dict[str, np.ndarray],
    p: float,
    q: float,
    r: float
) -> np.ndarray:
    """
    Build penalty matrix from precomputed masks and parameter values.
    
    Penalty values:
        - Self (diagonal): 1000
        - Same cell type, common: 0
        - Same cell type, non-common: p
        - Different cell, same tissue: q
        - Different cell, different tissue: r
    
    Args:
        masks: Precomputed boolean masks from precompute_penalty_masks()
        p: Penalty for same cell (non-common) in different tissue
        q: Penalty for different cell in same tissue
        r: Penalty for different cell in different tissue
        
    Returns:
        Penalty matrix of shape (m, m)
    """
    m = masks['self'].shape[0]
    penalty = np.zeros((m, m), dtype=np.float64)
    
    penalty[masks['self']] = 1000.0
    # same_cell_common stays 0 (default)
    penalty[masks['same_cell_noncommon']] = p
    penalty[masks['diff_cell_same_tissue']] = q
    penalty[masks['diff_cell_diff_tissue']] = r
    
    return penalty


# =============================================================================
# Score Computation
# =============================================================================

def compute_top_k_indices(objective_matrix: np.ndarray, k: int = 10) -> np.ndarray:
    """
    Get indices of top-k scoring genes for each cell.
    
    Args:
        objective_matrix: Score matrix of shape (num_cells, num_genes)
        k: Number of top genes to select
        
    Returns:
        Array of shape (num_cells, k) with top-k gene indices per cell
    """
    # argsort gives ascending order, so we take the last k and reverse
    return np.argsort(objective_matrix, axis=1)[:, -k:][:, ::-1]


def compute_scores(
    objective_matrix: np.ndarray,
    top_k_indices: np.ndarray,
    positives: Dict[int, List[int]],
    negatives: Dict[int, List[int]],
    highly_expressed: Set[int],
    literature_lookup: Dict[int, int]
) -> Tuple[int, int, int, int]:
    """
    Compute all four scoring metrics for a given objective matrix.
    
    Args:
        objective_matrix: Score matrix of shape (num_cells, num_genes)
        top_k_indices: Precomputed top-k indices per cell
        positives: Dict mapping cell index to positive marker indices
        negatives: Dict mapping cell index to negative marker indices
        highly_expressed: Set of highly expressed gene indices to penalize
        literature_lookup: Dict mapping cell index to expected literature marker index
        
    Returns:
        Tuple of (count1, count2, count3, count4) scores
    """
    # Count 1: Positive markers (from databases) appearing in top-k
    count1 = 0
    for cell_idx, pos_markers in positives.items():
        top_k_set = set(top_k_indices[cell_idx])
        count1 += len(set(pos_markers) & top_k_set)
    
    # Count 2: Literature-validated markers appearing in top-k
    count2 = 0
    for cell_idx, marker_idx in literature_lookup.items():
        if marker_idx in top_k_indices[cell_idx]:
            count2 += 1
    
    # Count 3: Negative markers appearing in top-k (penalized)
    count3 = 0
    for cell_idx, neg_markers in negatives.items():
        top_k_set = set(top_k_indices[cell_idx])
        count3 += len(set(neg_markers) & top_k_set)
    
    # Count 4: Highly expressed genes in top-k across all cells (penalized)
    highly_expressed_arr = np.array(list(highly_expressed))
    count4 = int(np.sum(np.isin(top_k_indices, highly_expressed_arr)))
    
    return count1, count2, count3, count4


def build_literature_lookup(
    positives2: List,
    gene_list: List,
    tissue_cells: List[List[str]]
) -> Dict[int, int]:
    """
    Build lookup table for literature-validated markers.
    
    Args:
        positives2: Literature marker data
        gene_list: List of genes
        tissue_cells: List of [tissue, cell_type] pairs
        
    Returns:
        Dict mapping cell index to marker gene index
    """
    lookup = {}
    for row in positives2:
        gene_name = row[0]
        target_pair = [row[1], row[2]]
        
        if [gene_name] in gene_list and target_pair in tissue_cells:
            marker_idx = gene_list.index([gene_name])
            cell_idx = tissue_cells.index(target_pair)
            lookup[cell_idx] = marker_idx
    
    return lookup


# =============================================================================
# Grid Search Optimization
# =============================================================================

def run_grid_search(
    tissue_cells: List[List[str]],
    common_cells: List,
    nTPM_matrix: np.ndarray,
    positives: Dict[int, List[int]],
    negatives: Dict[int, List[int]],
    highly_expressed: List[int],
    gene_list: List,
    positives2: List,
    p_range: Tuple[float, float, int] = (-700, -500, 20),
    q_range: Tuple[float, float, int] = (-700, -500, 20),
    r_range: Tuple[float, float, int] = (-700, -20, 20),
    verbose: bool = True
) -> Tuple[Tuple[float, float, float], Dict]:
    """
    Run grid search to find optimal penalty parameters p, q, r.
    
    The objective is to maximize:
        score = count1 + count2 - count3 - count4
    
    Where:
        - count1: Database markers appearing in top-10 (reward)
        - count2: Literature markers appearing in top-10 (reward)
        - count3: Negative markers appearing in top-10 (penalty)
        - count4: Highly expressed genes in top-10 (penalty)
    
    Args:
        tissue_cells: List of [tissue, cell_type] pairs
        common_cells: List of common cell types
        nTPM_matrix: Expression matrix (cells x genes)
        positives: Positive label indices per cell
        negatives: Negative label indices per cell
        highly_expressed: Indices of highly expressed genes to avoid
        gene_list: List of genes
        positives2: Literature-validated markers
        p_range: (start, stop, num_points) for parameter p
        q_range: (start, stop, num_points) for parameter q
        r_range: (start, stop, num_points) for parameter r
        verbose: Print progress information
        
    Returns:
        Tuple of (best_params, all_results_dict)
    """
    # Precompute masks (only depends on tissue/cell structure)
    masks = precompute_penalty_masks(tissue_cells, common_cells)
    
    # Build literature lookup table (only needs to be done once)
    literature_lookup = build_literature_lookup(positives2, gene_list, tissue_cells)
    
    # Convert highly_expressed to set for faster lookup
    highly_expressed_set = set(highly_expressed)
    
    # Generate parameter grid
    p_values = np.linspace(*p_range)
    q_values = np.linspace(*q_range)
    r_values = np.linspace(*r_range)
    
    total_combinations = len(p_values) * len(q_values) * len(r_values)
    if verbose:
        print(f"  Running grid search over {total_combinations} parameter combinations...")
    
    # Grid search
    results = {}
    best_score = float('-inf')
    best_params = None
    
    for p, q, r in product(p_values, q_values, r_values):
        # Build penalty and objective matrices
        penalty_matrix = build_penalty_matrix(masks, p, q, r)
        objective_matrix = penalty_matrix @ nTPM_matrix
        
        # Get top-10 indices for all cells
        top_k_indices = compute_top_k_indices(objective_matrix, k=10)
        
        # Compute scores
        count1, count2, count3, count4 = compute_scores(
            objective_matrix, top_k_indices,
            positives, negatives,
            highly_expressed_set, literature_lookup
        )
        
        score = count1 + count2 - count3 - count4
        results[(p, q, r)] = (count1, count2, count3, count4)
        
        if score > best_score:
            best_score = score
            best_params = (p, q, r)
    
    if verbose and best_params:
        counts = results[best_params]
        print(f"  Optimal parameters: p={best_params[0]:.2f}, q={best_params[1]:.2f}, r={best_params[2]:.2f}")
        print(f"  Organ-wide marker hits: {counts[0]}")
        print(f"  Whole-body marker hits: {counts[1]}")
        print(f"  Non-marker hits (penalized): {counts[2]}")
        print(f"  Highly expressed hits (penalized): {counts[3]}")
        print(f"  Total score: {best_score}")
    
    return best_params, results


# =============================================================================
# Marker Recommendation
# =============================================================================

def recommend_markers(
    objective_matrix: np.ndarray,
    tissue_cells: List[List[str]],
    gene_list: List,
    top_k: int = 10
) -> Dict[str, List[str]]:
    """
    Generate top-k marker recommendations for each cell type.
    
    Args:
        objective_matrix: Score matrix (cells x genes)
        tissue_cells: List of [tissue, cell_type] pairs
        gene_list: List of genes (possibly nested as [[gene1], [gene2], ...])
        top_k: Number of top markers to recommend
        
    Returns:
        Dictionary mapping cell name to list of recommended marker genes
    """
    # Handle nested gene list format
    if gene_list and isinstance(gene_list[0], list):
        genes = gene_list
    else:
        genes = [[g] for g in gene_list]
    
    # Get top-k indices
    top_k_indices = compute_top_k_indices(objective_matrix, k=top_k)
    
    # Build recommendations
    recommendations = {}
    for i, (tissue, cell_type) in enumerate(tissue_cells):
        cell_name = f"{tissue} {cell_type}"
        marker_indices = top_k_indices[i]
        recommendations[cell_name] = [genes[idx][0] for idx in marker_indices]
    
    return recommendations


def save_recommendations_csv(
    recommendations: Dict[str, List[str]],
    output_path: str
):
    """
    Save marker recommendations to CSV file.
    
    Args:
        recommendations: Dictionary mapping cell name to marker list
        output_path: Output file path
    """
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["cell", "markers"])
        for cell_name, markers in recommendations.items():
            writer.writerow([cell_name, markers])


# =============================================================================
# Label Loading
# =============================================================================

def load_labels_indexed(
    data_dir: str,
    gene_list: List,
    tissue_cells: List[List[str]]
) -> Tuple[Dict[int, List[int]], Dict[int, List[int]]]:
    """
    Load positive and negative labels with index mappings.
    
    Args:
        data_dir: Directory containing label CSV files
        gene_list: List of genes
        tissue_cells: List of [tissue, cell_type] pairs
        
    Returns:
        Tuple of (positives, negatives) dictionaries
        Each maps cell index to list of gene indices
    """
    data_path = Path(data_dir)
    
    # Flatten gene_list if nested
    if gene_list and isinstance(gene_list[0], list):
        gene_list_flat = [gene[0] for gene in gene_list]
    else:
        gene_list_flat = gene_list
    
    # Create index mappings
    tissue_cell_to_index = {tuple(row): idx for idx, row in enumerate(tissue_cells)}
    gene_to_index = {gene: idx for idx, gene in enumerate(gene_list_flat)}
    
    # Load positive labels
    positives_df = pd.read_csv(data_path / 'positives_labels.csv', sep='\t')
    positives = {}
    for _, row in positives_df.iterrows():
        key = (row["Tissue"], row["Cell type"])
        if key not in tissue_cell_to_index:
            continue
        
        cell_idx = tissue_cell_to_index[key]
        gene_names = eval(row["Positive Gene Names"])
        gene_indices = [gene_to_index[g] for g in gene_names if g in gene_to_index]
        
        if gene_indices:
            positives[cell_idx] = gene_indices
    
    # Load negative labels
    negatives_df = pd.read_csv(data_path / 'negative_labels.csv', sep='\t')
    negatives = {}
    for _, row in negatives_df.iterrows():
        key = (row["Tissue"], row["Cell type"])
        if key not in tissue_cell_to_index:
            continue
        
        cell_idx = tissue_cell_to_index[key]
        gene_names = eval(row["Negative Gene Names"])
        gene_indices = [gene_to_index[g] for g in gene_names if g in gene_to_index]
        
        if gene_indices:
            negatives[cell_idx] = gene_indices
    
    return positives, negatives


# =============================================================================
# Main Pipeline
# =============================================================================

def run_controlled_learning(
    data_dir: str,
    output_dir: str,
    run_median: bool = False
) -> Dict:
    """
    Run the complete controlled learning pipeline.
    
    Args:
        data_dir: Directory containing input data files
        output_dir: Directory for output files
        run_median: Whether to also run the median expression method
        
    Returns:
        Dictionary with results for each method
    """
    data_path = Path(data_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    print("=== Controlled Learning Pipeline ===\n")
    
    # Load data
    print("Loading data files...")
    tissue_cells = load_tissue_cells(data_path / "tissue_cell_pairs.tsv")
    common_cells = load_common_cells(data_path / "common_cells_across_tissues.csv")
    gene_list = load_gene_list(data_path / "gene_list.csv")
    
    # Load expression matrices
    nTPM_high_df = pd.read_csv(data_path / 'gene_expression_matrix_high.csv', sep='\t')
    nTPM_matrix_high = nTPM_high_df.drop(columns=["Tissue", "Cell type"]).values
    
    print(f"  Loaded {len(tissue_cells)} tissue-cell pairs")
    print(f"  Loaded {len(gene_list)} genes")
    
    # Load labels
    print("\nLoading labels...")
    positives, negatives = load_labels_indexed(data_dir, gene_list, tissue_cells)
    print(f"  Positive labels: {len(positives)} cells")
    print(f"  Negative labels: {len(negatives)} cells")
    
    # Compute highly expressed genes (top 50 by median expression)
    medians_high = np.median(nTPM_matrix_high, axis=0)
    highly_expressed_high = np.argsort(medians_high)[-50:].tolist()
    
    # Literature positives
    positives2 = LITERATURE_POSITIVES
    
    results = {}
    
    # === Process HIGH method ===
    print("\n--- Processing nTPM HIGH method ---")
    print("\nOptimizing parameters...")
    
    # Precompute penalty masks
    masks = precompute_penalty_masks(tissue_cells, common_cells)
    
    best_params_high, _ = run_grid_search(
        tissue_cells, common_cells, nTPM_matrix_high,
        positives, negatives, highly_expressed_high, gene_list, positives2
    )
    
    # Build final objective matrix
    penalty_matrix = build_penalty_matrix(masks, *best_params_high)
    objective_matrix_high = penalty_matrix @ nTPM_matrix_high
    
    # Generate and save recommendations
    print("\nGenerating marker recommendations...")
    recommendations_high = recommend_markers(objective_matrix_high, tissue_cells, gene_list)
    save_recommendations_csv(recommendations_high, output_path / "recommended_whole_body_markers_high.csv")
    print(f"  Saved recommendations for {len(recommendations_high)} cell types")
    
    results['high'] = {
        'optimal_params': best_params_high,
        'recommendations': recommendations_high
    }
    
    # === Process MEDIAN method (optional) ===
    if run_median:
        print("\n--- Processing nTPM MEDIAN method ---")
        
        nTPM_median_df = pd.read_csv(data_path / 'gene_expression_matrix_median.csv', sep='\t')
        nTPM_matrix_median = nTPM_median_df.drop(columns=["Tissue", "Cell type"]).values
        
        medians_median = np.median(nTPM_matrix_median, axis=0)
        highly_expressed_median = np.argsort(medians_median)[-50:].tolist()
        
        print("\nOptimizing parameters...")
        best_params_median, _ = run_grid_search(
            tissue_cells, common_cells, nTPM_matrix_median,
            positives, negatives, highly_expressed_median, gene_list, positives2
        )
        
        penalty_matrix = build_penalty_matrix(masks, *best_params_median)
        objective_matrix_median = penalty_matrix @ nTPM_matrix_median
        
        print("\nGenerating marker recommendations...")
        recommendations_median = recommend_markers(objective_matrix_median, tissue_cells, gene_list)
        save_recommendations_csv(recommendations_median, output_path / "recommended_whole_body_markers_median.csv")
        print(f"  Saved recommendations for {len(recommendations_median)} cell types")
        
        results['median'] = {
            'optimal_params': best_params_median,
            'recommendations': recommendations_median
        }
    
    print("\n=== Controlled Learning Complete ===\n")
    
    return results


# =============================================================================
# CLI Entry Point
# =============================================================================

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Run controlled learning for cell surface marker prediction"
    )
    parser.add_argument(
        "--data-dir", default=".",
        help="Directory containing processed data files"
    )
    parser.add_argument(
        "--output-dir", default=".",
        help="Directory for output files"
    )
    parser.add_argument(
        "--run-median", action="store_true",
        help="Also run the median expression method"
    )
    
    args = parser.parse_args()
    run_controlled_learning(args.data_dir, args.output_dir, args.run_median)
