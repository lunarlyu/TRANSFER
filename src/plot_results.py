"""
Plot Results Module

This module generates analysis plots and Excel reports from the
marker prediction results. It includes:
- Within-cell analysis: Expression rank of recommended markers
- Across-cell analysis: Target/off-target expression ratios
"""

import ast
from pathlib import Path
from typing import List, Dict

import numpy as np
import pandas as pd

# Use non-interactive backend for headless environments
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from .io_utils import load_tissue_cells, load_gene_list, load_expression_matrices


def within_cell_analysis(
    tissue_cells: List,
    genes: List,
    nTPM_matrix: np.ndarray,
    markers_file: str,
    output_dir: str,
    method_name: str = "high"
) -> List[int]:
    """
    Analyze expression rank of recommended markers within each cell.
    
    For each cell, find the rank of the top recommended marker's expression
    level compared to all other genes in that cell.
    
    Args:
        tissue_cells: List of [tissue, cell] pairs
        genes: List of gene names
        nTPM_matrix: Expression matrix
        markers_file: Path to recommended markers CSV
        output_dir: Directory to save plots
        method_name: Name for output files (e.g., "high" or "median")
        
    Returns:
        List of ranks for each cell
    """
    m = len(tissue_cells)
    cell_strs = [f"{t} {c}" for t, c in tissue_cells]
    
    # Load recommended markers
    rec_df = pd.read_csv(markers_file)
    rec_df['markers'] = rec_df['markers'].apply(lambda x: ast.literal_eval(x))
    
    # Build mapping of cell_id -> marker_ids
    top_markers = {}
    for i in range(min(m, len(rec_df))):
        cell = rec_df.iloc[i, 0]
        if cell not in cell_strs:
            continue
        cell_id = cell_strs.index(cell)
        markers = []
        for j in range(min(10, len(rec_df.iloc[i, 1]))):
            marker = rec_df.iloc[i, 1][j]
            if [marker] in genes:
                marker_id = genes.index([marker])
                markers.append(marker_id)
        if markers:
            top_markers[cell_id] = markers
    
    # Calculate ranks
    ranks = []
    for cell_id in range(m):
        if cell_id not in top_markers or not top_markers[cell_id]:
            ranks.append(0)
            continue
            
        row_data = nTPM_matrix[cell_id, :].copy()
        marker_id = top_markers[cell_id][0]  # Top 1 marker
        sorted_row = sorted(row_data, reverse=True)
        rank = sorted_row.index(row_data[marker_id])
        ranks.append(rank)
    
    # Create histogram plot
    output_path = Path(output_dir)
    bins = np.arange(0, 1000, 10)
    counts, bin_edges = np.histogram(ranks, bins=bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(bin_centers, counts, marker='o', linestyle='-', color='b')
    ax.set_xlabel('Ranks')
    ax.set_ylabel('Counts')
    ax.set_title(f'Rank of Expression Level: Recommended Top Marker ({method_name})')
    
    # Save as PNG (more universally compatible)
    fig.savefig(output_path / f'topmarkers_withinCell_{method_name}.png', 
                format='png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    
    # Print summary statistics
    count_bins = [0, 0, 0, 0, 0]  # <=10, 11-20, 21-30, 31-40, >40
    for r in ranks:
        if r <= 10:
            count_bins[0] += 1
        elif r <= 20:
            count_bins[1] += 1
        elif r <= 30:
            count_bins[2] += 1
        elif r <= 40:
            count_bins[3] += 1
        else:
            count_bins[4] += 1
    
    print(f"  Within-cell analysis ({method_name}):")
    print(f"    Rank distribution: {count_bins}")
    print(f"    Total cells: {len(ranks)}")
    
    return ranks


def across_cell_analysis_top1(
    tissue_cells: List,
    genes: List,
    nTPM_matrix: np.ndarray,
    markers_file: str,
    output_dir: str,
    method_name: str = "high"
) -> None:
    """
    Analyze target/off-target ratio for top 1 marker across cells.
    
    Args:
        tissue_cells: List of [tissue, cell] pairs
        genes: List of gene names
        nTPM_matrix: Expression matrix
        markers_file: Path to recommended markers CSV
        output_dir: Directory to save results
        method_name: Name for output files
    """
    m = len(tissue_cells)
    cell_strs = [f"{t} {c}" for t, c in tissue_cells]
    
    rec_df = pd.read_csv(markers_file)
    rec_df['markers'] = rec_df['markers'].apply(lambda x: ast.literal_eval(x))
    
    # Build mapping
    top_markers = {}
    for i in range(min(m, len(rec_df))):
        cell = rec_df.iloc[i, 0]
        if cell not in cell_strs:
            continue
        cell_id = cell_strs.index(cell)
        marker = rec_df.iloc[i, 1][0]
        if [marker] in genes:
            marker_id = genes.index([marker])
            top_markers[cell_id] = marker_id
    
    # Calculate ratios
    ratios = []
    for i in range(m):
        if i not in top_markers:
            ratios.append(1.0)
            continue
            
        marker_id = top_markers[i]
        column_data = nTPM_matrix[:, marker_id].copy()
        
        # Delete all cells of the same cell type
        cell_type = tissue_cells[i][1]
        delete_idxs = [idx for idx in range(m) if tissue_cells[idx][1] == cell_type]
        column_data = np.delete(column_data, delete_idxs)
        
        off_exp = np.max(column_data) if len(column_data) > 0 else 0
        targ_exp = nTPM_matrix[i, marker_id]
        
        if off_exp == 0 and targ_exp == 0:
            ratios.append(1.0)
        elif off_exp == 0:
            ratios.append(float('inf'))
        else:
            ratios.append(targ_exp / off_exp)
    
    # Sort and save
    sorted_indices = np.argsort(ratios)[::-1]
    x_axis = [tissue_cells[idx][0] + " " + tissue_cells[idx][1] for idx in sorted_indices]
    y_axis = [ratios[idx] for idx in sorted_indices]
    
    df = pd.DataFrame({
        'Organ-cell': x_axis,
        'Target/Off-Target Ratio': y_axis
    })
    
    output_path = Path(output_dir)
    df.to_excel(output_path / f'top1marker_{method_name}.xlsx', index=False)
    print(f"  Saved: top1marker_{method_name}.xlsx")


def across_cell_analysis_top10(
    tissue_cells: List,
    genes: List,
    nTPM_matrix: np.ndarray,
    markers_file: str,
    output_dir: str,
    method_name: str = "high"
) -> None:
    """
    Analyze best target/off-target ratio among top 10 markers.
    
    Args:
        tissue_cells: List of [tissue, cell] pairs
        genes: List of gene names
        nTPM_matrix: Expression matrix
        markers_file: Path to recommended markers CSV
        output_dir: Directory to save results
        method_name: Name for output files
    """
    m = len(tissue_cells)
    cell_strs = [f"{t} {c}" for t, c in tissue_cells]
    
    rec_df = pd.read_csv(markers_file)
    rec_df['markers'] = rec_df['markers'].apply(lambda x: ast.literal_eval(x))
    
    # Build mapping
    top_markers = {}
    for i in range(min(m, len(rec_df))):
        cell = rec_df.iloc[i, 0]
        if cell not in cell_strs:
            continue
        cell_id = cell_strs.index(cell)
        markers = []
        for j in range(min(10, len(rec_df.iloc[i, 1]))):
            marker = rec_df.iloc[i, 1][j]
            if [marker] in genes:
                marker_id = genes.index([marker])
                markers.append(marker_id)
        if markers:
            top_markers[cell_id] = markers
    
    # Calculate best ratio for each cell
    ratios = []
    for i in range(m):
        if i not in top_markers:
            ratios.append(1.0)
            continue
            
        markers = top_markers[i]
        cell_type = tissue_cells[i][1]
        delete_idxs = [idx for idx in range(m) if tissue_cells[idx][1] == cell_type]
        
        ratio_list = []
        for marker_id in markers:
            column_data = nTPM_matrix[:, marker_id].copy()
            column_data = np.delete(column_data, delete_idxs)
            
            off_exp = np.max(column_data) if len(column_data) > 0 else 0
            targ_exp = nTPM_matrix[i, marker_id]
            
            if off_exp == 0 and targ_exp == 0:
                ratio_list.append(1.0)
            elif off_exp == 0:
                ratio_list.append(float('inf'))
            else:
                ratio_list.append(targ_exp / off_exp)
        
        ratios.append(max(ratio_list) if ratio_list else 1.0)
    
    # Sort and save
    sorted_indices = np.argsort(ratios)[::-1]
    x_axis = [tissue_cells[idx][0] + " " + tissue_cells[idx][1] for idx in sorted_indices]
    y_axis = [ratios[idx] for idx in sorted_indices]
    
    df = pd.DataFrame({
        'Organ-cell': x_axis,
        'Highest Target/Off-Target Ratio': y_axis
    })
    
    output_path = Path(output_dir)
    df.to_excel(output_path / f'top10marker_{method_name}.xlsx', index=False)
    print(f"  Saved: top10marker_{method_name}.xlsx")


def across_cell_analysis_top2_combination(
    tissue_cells: List,
    genes: List,
    nTPM_matrix: np.ndarray,
    markers_file: str,
    output_dir: str,
    method_name: str = "high"
) -> None:
    """
    Analyze best 2-marker combination from top 10 markers.
    
    Finds the pair of markers with highest product of target/off-target ratios.
    
    Args:
        tissue_cells: List of [tissue, cell] pairs
        genes: List of gene names
        nTPM_matrix: Expression matrix
        markers_file: Path to recommended markers CSV
        output_dir: Directory to save results
        method_name: Name for output files
    """
    m = len(tissue_cells)
    cell_strs = [f"{t} {c}" for t, c in tissue_cells]
    
    rec_df = pd.read_csv(markers_file)
    rec_df['markers'] = rec_df['markers'].apply(lambda x: ast.literal_eval(x))
    
    # Build mapping
    top_markers = {}
    for i in range(min(m, len(rec_df))):
        cell = rec_df.iloc[i, 0]
        if cell not in cell_strs:
            continue
        cell_id = cell_strs.index(cell)
        markers = []
        for j in range(min(10, len(rec_df.iloc[i, 1]))):
            marker = rec_df.iloc[i, 1][j]
            if [marker] in genes:
                marker_id = genes.index([marker])
                markers.append(marker_id)
        if markers:
            top_markers[cell_id] = markers
    
    # Calculate best 2-marker combination for each cell
    ratios = []
    for i in range(m):
        if i not in top_markers or len(top_markers[i]) < 2:
            ratios.append(1.0)
            continue
            
        markers = top_markers[i]
        cell_type = tissue_cells[i][1]
        delete_idxs = [idx for idx in range(m) if tissue_cells[idx][1] == cell_type]
        
        # Pre-compute ratios for all markers
        marker_ratios = []
        for marker_id in markers:
            column_data = nTPM_matrix[:, marker_id].copy()
            column_data = np.delete(column_data, delete_idxs)
            
            off_exp = np.max(column_data) if len(column_data) > 0 else 0
            targ_exp = nTPM_matrix[i, marker_id]
            
            if off_exp == 0 and targ_exp == 0:
                marker_ratios.append(1.0)
            elif off_exp == 0:
                marker_ratios.append(float('inf'))
            else:
                marker_ratios.append(targ_exp / off_exp)
        
        # Find best pair
        best_product = 0
        for k in range(len(markers)):
            for l in range(k + 1, len(markers)):
                product = marker_ratios[k] * marker_ratios[l]
                if product > best_product:
                    best_product = product
        
        ratios.append(best_product if best_product > 0 else 1.0)
    
    # Sort and save
    sorted_indices = np.argsort(ratios)[::-1]
    x_axis = [tissue_cells[idx][0] + " " + tissue_cells[idx][1] for idx in sorted_indices]
    y_axis = [ratios[idx] for idx in sorted_indices]
    
    df = pd.DataFrame({
        'Organ-cell': x_axis,
        'Highest Target/Off-Target Ratio Product': y_axis
    })
    
    output_path = Path(output_dir)
    df.to_excel(output_path / f'top2_10marker_{method_name}.xlsx', index=False)
    print(f"  Saved: top2_10marker_{method_name}.xlsx")


def run_plotting(
    data_dir: str,
    output_dir: str,
    plots_dir: str,
    verbose: bool = True
) -> None:
    """
    Run all plotting and analysis tasks.
    
    Args:
        data_dir: Directory containing processed data files
        output_dir: Directory containing marker recommendation files
        plots_dir: Directory to save plots and Excel files
        verbose: Whether to print progress
    """
    if verbose:
        print("\n" + "=" * 60)
        print("GENERATING PLOTS AND ANALYSIS")
        print("=" * 60)
    
    # Create plots directory
    plots_path = Path(plots_dir)
    plots_path.mkdir(parents=True, exist_ok=True)
    
    # Load data using shared utilities
    if verbose:
        print("\nLoading data...")
    
    data_path = Path(data_dir)
    tissue_cells = load_tissue_cells(data_path / "tissue_cell_pairs.tsv")
    genes = load_gene_list(data_path / "gene_list.csv")
    nTPM_high, nTPM_median, _, _ = load_expression_matrices(data_dir)
    
    if verbose:
        print(f"  Loaded {len(tissue_cells)} tissue-cell pairs")
        print(f"  Loaded {len(genes)} genes")
    
    output_path = Path(output_dir)
    
    # Run analyses for HIGH method
    if verbose:
        print("\nAnalyzing HIGH method results...")
    
    markers_high = output_path / 'recommended_whole_body_markers_high.csv'
    if markers_high.exists():
        within_cell_analysis(tissue_cells, genes, nTPM_high, str(markers_high), 
                           plots_dir, "high")
        across_cell_analysis_top1(tissue_cells, genes, nTPM_high, str(markers_high),
                                 plots_dir, "high")
        across_cell_analysis_top10(tissue_cells, genes, nTPM_high, str(markers_high),
                                  plots_dir, "high")
        across_cell_analysis_top2_combination(tissue_cells, genes, nTPM_high, 
                                             str(markers_high), plots_dir, "high")
    else:
        print(f"  Warning: {markers_high} not found, skipping HIGH analysis")
    
    # Run analyses for MEDIAN method
    if verbose:
        print("\nAnalyzing MEDIAN method results...")
    
    markers_median = output_path / 'recommended_whole_body_markers_median.csv'
    if markers_median.exists():
        within_cell_analysis(tissue_cells, genes, nTPM_median, str(markers_median),
                           plots_dir, "median")
        across_cell_analysis_top1(tissue_cells, genes, nTPM_median, str(markers_median),
                                 plots_dir, "median")
        across_cell_analysis_top10(tissue_cells, genes, nTPM_median, str(markers_median),
                                  plots_dir, "median")
        across_cell_analysis_top2_combination(tissue_cells, genes, nTPM_median,
                                             str(markers_median), plots_dir, "median")
    else:
        print(f"  Warning: {markers_median} not found, skipping MEDIAN analysis")
    
    if verbose:
        print("\n" + "=" * 60)
        print(f"Plots and analysis saved to: {plots_dir}")
        print("=" * 60)


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Generate plots and analysis from marker prediction results"
    )
    parser.add_argument(
        "--data-dir", 
        default="Results",
        help="Directory containing processed data files (default: Results)"
    )
    parser.add_argument(
        "--output-dir",
        default="Results",
        help="Directory containing marker recommendation files (default: Results)"
    )
    parser.add_argument(
        "--plots-dir",
        default="Plots",
        help="Directory to save plots and Excel files (default: Plots)"
    )
    
    args = parser.parse_args()
    
    run_plotting(
        data_dir=args.data_dir,
        output_dir=args.output_dir,
        plots_dir=args.plots_dir
    )
