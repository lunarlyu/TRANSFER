"""
Controlled Learning Module

This module implements the penalty matrix optimization algorithm for
cell surface marker prediction. It includes functions for:
- Building penalty matrices based on tissue-cell relationships
- Grid search optimization for penalty parameters
- Computing objective scores and recommending markers
- Saving results to CSV files
"""

import csv
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd


# Known positive markers from literature (for validation)
LITERATURE_POSITIVES = [
    ["PECAM1", "lung", "endothelial cells", "LNP"],
    ["VCAM1", "vascular", "endothelial cells", "LNP"],
    ["CD4", "pbmc", "t-cells", "LNP"],
    ["CD5", "pbmc", "t-cells", "LNP"],
    ["CD19", "pbmc", "b-cells", "LNP"],
    ["CD3", "pbmc", "t-cells", "LNP"],
    ["NCR1", "pbmc", "nk-cells", "LNP"],
    ["CD14", "pbmc", "macrophages", "LNP"],
    ["MRC1", "pbmc", "macrophages", "LNP"],
    ["ITGAM", "pbmc", "macrophages", ""],
    ["CD28", "pbmc", "t-cells", "EDV"],
    ["CD40", "pbmc", "b-cells", "lenti"],
    ["ENG", "vascular", "endothelial cells", "LNP"],
    ["MRC1", "pbmc", "dendritic cells", "LNP"],
    ["CD8", "pbmc", "t-cells", "LNP"],
    ["PDPN", "skin", "endothelial cells", "LNP"],
    ["PLVAP", "lung", "endothelial cells", "LNP?"],
    ["FCER2", "pbmc", "b-cells", ""],
]


def load_data_files(data_dir: str) -> Dict:
    """
    Load all required data files for controlled learning.
    
    Args:
        data_dir: Directory containing the processed data files
        
    Returns:
        Dictionary containing all loaded data
    """
    data_path = Path(data_dir)
    
    # Load tissue-cell pairs
    with open(data_path / "tissue_cell_pairs.tsv", newline='') as file:
        reader = csv.reader(file, delimiter='\t')
        tissue_cells = list(reader)
    tissue_cells = tissue_cells[1:]  # Remove header
    
    # Load common cells across tissues
    with open(data_path / "Original data" / "common_cells_across_tissues.csv", newline='') as file:
        reader = csv.reader(file, delimiter='\t')
        common_cells = list(reader)
    common_cells.append('serous glandular cells')
    
    # Load gene list
    with open(data_path / "gene_list.csv", newline='') as file:
        reader = csv.reader(file, delimiter='\t')
        gene_list = list(reader)
    
    # Load nTPM matrices
    nTPM_high_df = pd.read_csv(data_path / "gene_expression_matrix_high.csv", sep='\t')
    nTPM_high_only = nTPM_high_df.drop(columns=["Tissue", "Cell type"])
    nTPM_matrix_high = nTPM_high_only.values
    
    nTPM_median_df = pd.read_csv(data_path / "gene_expression_matrix_median.csv", sep='\t')
    nTPM_median_only = nTPM_median_df.drop(columns=["Tissue", "Cell type"])
    nTPM_matrix_median = nTPM_median_only.values
    
    return {
        'tissue_cells': tissue_cells,
        'common_cells': common_cells,
        'gene_list': gene_list,
        'nTPM_matrix_high': nTPM_matrix_high,
        'nTPM_matrix_median': nTPM_matrix_median,
        'nTPM_high_df': nTPM_high_df,
        'nTPM_median_df': nTPM_median_df
    }


def load_labels(data_dir: str, gene_list: List, tissue_cells: List) -> Tuple[Dict, Dict]:
    """
    Load and index positive and negative labels.
    
    Args:
        data_dir: Directory containing label files
        gene_list: List of genes
        tissue_cells: List of tissue-cell pairs
        
    Returns:
        Tuple of (positives dict, negatives dict) with integer indices
    """
    data_path = Path(data_dir)
    
    # Flatten gene list if needed
    if isinstance(gene_list[0], list):
        gene_list = [gene[0] for gene in gene_list]
    
    # Create index mappings
    tissue_cell_to_index = {tuple(row): idx for idx, row in enumerate(tissue_cells)}
    gene_to_index = {gene: idx for idx, gene in enumerate(gene_list)}
    
    # Load positive labels
    positives_df = pd.read_csv(data_path / "positives_labels.csv", sep='\t')
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
    negatives_df = pd.read_csv(data_path / "negative_labels.csv", sep='\t')
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
    
    Args:
        nTPM_matrix: Expression matrix
        top_n: Number of top genes to return
        
    Returns:
        List of gene indices
    """
    medians = np.median(nTPM_matrix, axis=0).tolist()
    sorted_indices = [idx for idx, _ in sorted(enumerate(medians), key=lambda x: x[1], reverse=True)]
    return sorted_indices[:top_n]


def build_penalty_matrix(
    tissue_cells: List,
    common_cells: List,
    p: float,
    q: float,
    r: float
) -> np.ndarray:
    """
    Build the penalty matrix for objective function computation.
    
    The penalty matrix encodes relationships between tissue-cell pairs:
    - Same cell: 1000 (self)
    - Same cell (common) in different tissue: 0
    - Same cell (non-common) in different tissue: p
    - Different cell in same tissue: q
    - Different cell in different tissue: r
    
    Args:
        tissue_cells: List of [tissue, cell] pairs
        common_cells: List of common cell types across tissues
        p: Penalty for same non-common cell in different tissue
        q: Penalty for different cell in same tissue
        r: Penalty for different cell in different tissue
        
    Returns:
        Penalty matrix of shape (n_cells, n_cells)
    """
    m = len(tissue_cells)
    penalty_matrix = []
    
    # Flatten common_cells if nested
    common_cells_flat = []
    for cell in common_cells:
        if isinstance(cell, list):
            common_cells_flat.extend(cell)
        else:
            common_cells_flat.append(cell)
    
    for i in range(m):
        row = []
        for j in range(m):
            if j == i:
                row.append(1000)  # Self
            elif tissue_cells[j][1] == tissue_cells[i][1]:
                # Same cell type
                if tissue_cells[j][1] in common_cells_flat:
                    row.append(0)  # Common cell, different tissue
                else:
                    row.append(p)  # Non-common cell, different tissue
            elif tissue_cells[j][0] == tissue_cells[i][0]:
                row.append(q)  # Different cell, same tissue
            else:
                row.append(r)  # Different cell, different tissue
        penalty_matrix.append(row)
    
    return np.array(penalty_matrix)


def compute_objective_matrix(
    penalty_matrix: np.ndarray,
    nTPM_matrix: np.ndarray
) -> np.ndarray:
    """
    Compute the objective matrix by multiplying penalty and nTPM matrices.
    
    Args:
        penalty_matrix: Penalty matrix
        nTPM_matrix: Expression matrix
        
    Returns:
        Objective matrix
    """
    return np.dot(penalty_matrix, nTPM_matrix)


def evaluate_parameters(
    objective_matrix: np.ndarray,
    positives: Dict[int, List[int]],
    negatives: Dict[int, List[int]],
    highly_expressed: List[int],
    gene_list: List,
    tissue_cells: List,
    literature_positives: List = LITERATURE_POSITIVES,
    top_k: int = 10
) -> Tuple[int, int, int, int]:
    """
    Evaluate parameters by counting hits in various categories.
    
    Args:
        objective_matrix: Computed objective matrix
        positives: Positive labels dict
        negatives: Negative labels dict
        highly_expressed: Indices of highly expressed genes to penalize
        gene_list: List of genes
        tissue_cells: List of tissue-cell pairs
        literature_positives: Known positive markers from literature
        top_k: Number of top markers to consider
        
    Returns:
        Tuple of (count1, count2, count3, count4) representing:
        - count1: Hits in organ-wide markers (from databases)
        - count2: Hits in whole-body markers (from literature)
        - count3: Hits in negative markers (to minimize)
        - count4: Hits in highly expressed genes (to minimize)
    """
    objective_list = objective_matrix.tolist()
    
    # Flatten gene list if needed
    if isinstance(gene_list[0], list):
        gene_list_flat = [gene[0] for gene in gene_list]
    else:
        gene_list_flat = gene_list
    
    # Count 1: Organ-wide markers (from databases)
    count1 = 0
    for cell in positives:
        markers = positives[cell]
        obj_row = objective_list[cell]
        top_markers = sorted(enumerate(obj_row), key=lambda x: x[1], reverse=True)[:top_k]
        suggested = [idx for idx, _ in top_markers]
        count1 += len(set(markers).intersection(set(suggested)))
    
    # Count 2: Whole-body markers (from literature)
    count2 = 0
    for row in literature_positives:
        gene_name = row[0]
        if [gene_name] in gene_list or gene_name in gene_list_flat:
            if [gene_name] in gene_list:
                marker_idx = gene_list.index([gene_name])
            else:
                marker_idx = gene_list_flat.index(gene_name)
            
            target_pair = [row[1], row[2]]
            if target_pair in tissue_cells:
                cell_idx = tissue_cells.index(target_pair)
                obj_row = objective_list[cell_idx]
                top_markers = sorted(enumerate(obj_row), key=lambda x: x[1], reverse=True)[:top_k]
                suggested = [idx for idx, _ in top_markers]
                if marker_idx in suggested:
                    count2 += 1
    
    # Count 3: Negative markers (to minimize)
    count3 = 0
    for cell in negatives:
        non_markers = negatives[cell]
        obj_row = objective_list[cell]
        top_markers = sorted(enumerate(obj_row), key=lambda x: x[1], reverse=True)[:top_k]
        suggested = [idx for idx, _ in top_markers]
        count3 += len(set(non_markers).intersection(set(suggested)))
    
    # Count 4: Highly expressed genes (to minimize)
    count4 = 0
    for obj_row in objective_list:
        top_markers = sorted(enumerate(obj_row), key=lambda x: x[1], reverse=True)[:top_k]
        suggested = [idx for idx, _ in top_markers]
    # Note: The original code only counts the last row's overlap
    for i in suggested:
        if i in highly_expressed:
            count4 += 1
    
    return count1, count2, count3, count4


def grid_search_parameters(
    tissue_cells: List,
    common_cells: List,
    nTPM_matrix: np.ndarray,
    positives: Dict[int, List[int]],
    negatives: Dict[int, List[int]],
    highly_expressed: List[int],
    gene_list: List,
    p_range: Tuple[float, float, int] = (-700, -500, 20),
    q_range: Tuple[float, float, int] = (-700, -500, 20),
    r_range: Tuple[float, float, int] = (-700, -20, 20),
    verbose: bool = True
) -> Tuple[Tuple[float, float, float], Dict]:
    """
    Perform grid search to find optimal penalty parameters.
    
    Args:
        tissue_cells: List of tissue-cell pairs
        common_cells: List of common cells across tissues
        nTPM_matrix: Expression matrix
        positives: Positive labels dict
        negatives: Negative labels dict
        highly_expressed: Indices of highly expressed genes
        gene_list: List of genes
        p_range: (start, end, num_points) for p parameter
        q_range: (start, end, num_points) for q parameter
        r_range: (start, end, num_points) for r parameter
        verbose: Whether to print progress
        
    Returns:
        Tuple of (optimal_params, all_results_dict)
    """
    results = {}
    best_score = 0
    best_params = None
    
    p_values = np.linspace(p_range[0], p_range[1], p_range[2])
    q_values = np.linspace(q_range[0], q_range[1], q_range[2])
    r_values = np.linspace(r_range[0], r_range[1], r_range[2])
    
    total_iterations = len(p_values) * len(q_values) * len(r_values)
    
    if verbose:
        print(f"  Running grid search over {total_iterations} parameter combinations...")
    
    for p in p_values:
        for q in q_values:
            for r in r_values:
                # Build penalty matrix
                penalty_matrix = build_penalty_matrix(tissue_cells, common_cells, p, q, r)
                
                # Compute objective matrix
                objective_matrix = compute_objective_matrix(penalty_matrix, nTPM_matrix)
                
                # Evaluate parameters
                count1, count2, count3, count4 = evaluate_parameters(
                    objective_matrix,
                    positives,
                    negatives,
                    highly_expressed,
                    gene_list,
                    tissue_cells
                )
                
                # Compute score (maximize positives, minimize negatives)
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
    
    return best_params, results


def recommend_markers(
    objective_matrix: np.ndarray,
    tissue_cells: List,
    gene_list: List,
    top_k: int = 10
) -> Dict[str, List[str]]:
    """
    Generate top-k marker recommendations for each cell type.
    
    Args:
        objective_matrix: Computed objective matrix
        tissue_cells: List of tissue-cell pairs
        gene_list: List of genes
        top_k: Number of markers to recommend
        
    Returns:
        Dictionary mapping "tissue cell" to list of recommended marker genes
    """
    # Flatten gene list if needed
    if isinstance(gene_list[0], list):
        genes = [gene[0] for gene in gene_list]
    else:
        genes = gene_list
    
    recommendations = {}
    
    for i, obj_row in enumerate(objective_matrix):
        top_markers = sorted(enumerate(obj_row), key=lambda x: x[1], reverse=True)[:top_k]
        marker_indices = [idx for idx, _ in top_markers]
        marker_names = [genes[idx] for idx in marker_indices]
        
        cell_name = f"{tissue_cells[i][0]} {tissue_cells[i][1]}"
        recommendations[cell_name] = marker_names
    
    return recommendations


def save_recommendations(
    recommendations: Dict[str, List[str]],
    output_path: str
):
    """
    Save marker recommendations to CSV file.
    
    Args:
        recommendations: Dictionary of recommendations
        output_path: Output file path
    """
    with open(output_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["cell", "markers"])
        for cell, markers in recommendations.items():
            writer.writerow([cell, markers])


def run_controlled_learning(data_dir: str, output_dir: str) -> Dict:
    """
    Run the complete controlled learning pipeline.
    
    Args:
        data_dir: Directory containing processed data files
        output_dir: Directory for output files
        
    Returns:
        Dictionary containing results
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    print("=== Controlled Learning Pipeline ===\n")
    
    # Load data
    print("Loading data files...")
    data = load_data_files(data_dir)
    
    tissue_cells = data['tissue_cells']
    common_cells = data['common_cells']
    gene_list = data['gene_list']
    nTPM_matrix_high = data['nTPM_matrix_high']
    nTPM_matrix_median = data['nTPM_matrix_median']
    
    print(f"  Loaded {len(tissue_cells)} tissue-cell pairs")
    print(f"  Loaded {len(gene_list)} genes")
    
    # Load labels
    print("\nLoading labels...")
    positives, negatives = load_labels(data_dir, gene_list, tissue_cells)
    print(f"  Positive labels: {len(positives)} cells")
    print(f"  Negative labels: {len(negatives)} cells")
    
    # Get highly expressed genes
    highly_expressed_high = get_highly_expressed_genes(nTPM_matrix_high)
    highly_expressed_median = get_highly_expressed_genes(nTPM_matrix_median)
    
    results = {}
    
    # === Process HIGH method ===
    print("\n--- Processing nTPM HIGH method ---")
    
    print("\nOptimizing parameters...")
    best_params_high, _ = grid_search_parameters(
        tissue_cells,
        common_cells,
        nTPM_matrix_high,
        positives,
        negatives,
        highly_expressed_high,
        gene_list
    )
    
    # Build final objective matrix with optimal parameters
    penalty_matrix_high = build_penalty_matrix(
        tissue_cells, common_cells,
        best_params_high[0], best_params_high[1], best_params_high[2]
    )
    objective_matrix_high = compute_objective_matrix(penalty_matrix_high, nTPM_matrix_high)
    
    # Generate recommendations
    print("\nGenerating marker recommendations...")
    recommendations_high = recommend_markers(objective_matrix_high, tissue_cells, gene_list)
    save_recommendations(recommendations_high, output_path / "recommended_whole_body_markers_high.csv")
    print(f"  Saved recommendations for {len(recommendations_high)} cell types")
    
    results['high'] = {
        'optimal_params': best_params_high,
        'recommendations': recommendations_high
    }
    
    # === Process MEDIAN method ===
    print("\n--- Processing nTPM MEDIAN method ---")
    
    print("\nOptimizing parameters...")
    best_params_median, _ = grid_search_parameters(
        tissue_cells,
        common_cells,
        nTPM_matrix_median,
        positives,
        negatives,
        highly_expressed_median,
        gene_list
    )
    
    # Build final objective matrix with optimal parameters
    penalty_matrix_median = build_penalty_matrix(
        tissue_cells, common_cells,
        best_params_median[0], best_params_median[1], best_params_median[2]
    )
    objective_matrix_median = compute_objective_matrix(penalty_matrix_median, nTPM_matrix_median)
    
    # Generate recommendations
    print("\nGenerating marker recommendations...")
    recommendations_median = recommend_markers(objective_matrix_median, tissue_cells, gene_list)
    save_recommendations(recommendations_median, output_path / "recommended_whole_body_markers_median.csv")
    print(f"  Saved recommendations for {len(recommendations_median)} cell types")
    
    results['median'] = {
        'optimal_params': best_params_median,
        'recommendations': recommendations_median
    }
    
    print("\n=== Controlled Learning Complete ===\n")
    
    return results


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Run controlled learning for marker prediction")
    parser.add_argument("--data-dir", default=".", help="Directory containing processed data")
    parser.add_argument("--output-dir", default=".", help="Directory for output files")
    
    args = parser.parse_args()
    run_controlled_learning(args.data_dir, args.output_dir)

