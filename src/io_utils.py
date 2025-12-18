"""
I/O Utilities Module

This module provides shared file I/O functions used across the
cell surface marker prediction pipeline.
"""

import csv
from pathlib import Path
from typing import Dict, List, Tuple, Set

import numpy as np
import pandas as pd


def load_surface_proteins(file_path: str) -> Set[str]:
    """
    Load surface protein gene symbols from mass spectrometry data.
    
    Args:
        file_path: Path to the mass_spec_valid_surface_protein.csv file
        
    Returns:
        Set of gene symbols for validated surface proteins
    """
    df = pd.read_csv(file_path, header=0, index_col=0)
    return set(df["ENTREZ gene symbol"])


def load_tissue_cells(file_path: str) -> List[List[str]]:
    """
    Load tissue-cell pairs from TSV file.
    
    Args:
        file_path: Path to tissue_cell_pairs.tsv
        
    Returns:
        List of [tissue, cell_type] pairs
    """
    with open(file_path, newline='') as file:
        reader = csv.reader(file, delimiter='\t')
        tissue_cells = list(reader)
    return tissue_cells[1:]  # Skip header


def load_gene_list(file_path: str) -> List:
    """
    Load gene list from CSV file.
    
    Args:
        file_path: Path to gene_list.csv
        
    Returns:
        List of genes (each as a single-element list for compatibility)
    """
    with open(file_path, newline='') as file:
        reader = csv.reader(file, delimiter='\t')
        genes = list(reader)
    return genes


def load_common_cells(file_path: str) -> List:
    """
    Load common cells across tissues.
    
    Args:
        file_path: Path to common_cells_across_tissues.csv
        
    Returns:
        List of common cell types
    """
    with open(file_path, newline='') as file:
        reader = csv.reader(file, delimiter='\t')
        common_cells = list(reader)
    # Add 'serous glandular cells' as in original notebook
    common_cells.append('serous glandular cells')
    return common_cells


def load_expression_matrices(data_dir: str) -> Tuple[np.ndarray, np.ndarray, pd.DataFrame, pd.DataFrame]:
    """
    Load both HIGH and MEDIAN expression matrices.
    
    Args:
        data_dir: Directory containing the expression matrix CSV files
        
    Returns:
        Tuple of (nTPM_matrix_high, nTPM_matrix_median, nTPM_high_df, nTPM_median_df)
    """
    data_path = Path(data_dir)
    
    # Load HIGH matrix
    nTPM_high_df = pd.read_csv(data_path / "gene_expression_matrix_high.csv", sep='\t')
    nTPM_high_only = nTPM_high_df.drop(columns=["Tissue", "Cell type"])
    nTPM_matrix_high = nTPM_high_only.values
    
    # Load MEDIAN matrix
    nTPM_median_df = pd.read_csv(data_path / "gene_expression_matrix_median.csv", sep='\t')
    nTPM_median_only = nTPM_median_df.drop(columns=["Tissue", "Cell type"])
    nTPM_matrix_median = nTPM_median_only.values
    
    return nTPM_matrix_high, nTPM_matrix_median, nTPM_high_df, nTPM_median_df


def load_labels(
    data_dir: str,
    gene_list: List,
    tissue_cells: List
) -> Tuple[Dict[int, List[int]], Dict[int, List[int]]]:
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
        gene_list_flat = [gene[0] for gene in gene_list]
    else:
        gene_list_flat = gene_list
    
    # Create index mappings
    tissue_cell_to_index = {tuple(row): idx for idx, row in enumerate(tissue_cells)}
    gene_to_index = {gene: idx for idx, gene in enumerate(gene_list_flat)}
    
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


def save_labels_to_csv(
    labels: Dict[Tuple[str, str], List[str]],
    output_path: str,
    label_type: str = "Positive"
):
    """
    Save labels dictionary to CSV file.
    
    Args:
        labels: Labels dictionary mapping (tissue, cell) to gene list
        output_path: Output file path
        label_type: Type of labels ("Positive" or "Negative")
    """
    records = []
    for (tissue, cell_type), gene_names in labels.items():
        records.append({
            "Tissue": tissue,
            "Cell type": cell_type,
            f"{label_type} Gene Names": gene_names
        })
    
    df = pd.DataFrame(records)
    df.to_csv(output_path, sep='\t', index=False)


def save_recommendations(
    recommendations: Dict[str, List[str]],
    output_path: str
):
    """
    Save marker recommendations to CSV file.
    
    Args:
        recommendations: Dictionary mapping cell name to list of markers
        output_path: Output file path
    """
    with open(output_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["cell", "markers"])
        # Sort by cell name for deterministic output order
        for cell in sorted(recommendations.keys()):
            markers = recommendations[cell]
            writer.writerow([cell, markers])


def flatten_gene_list(gene_list: List) -> List[str]:
    """
    Flatten gene list if it contains nested lists.
    
    Args:
        gene_list: List of genes, possibly nested
        
    Returns:
        Flat list of gene names
    """
    if gene_list and isinstance(gene_list[0], list):
        return [gene[0] for gene in gene_list]
    return gene_list

