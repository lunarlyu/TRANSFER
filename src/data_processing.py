"""
Data Processing Module

This module handles processing and cleaning of gene expression data for
cell surface marker prediction. It includes functions for:
- Loading and cleaning gene expression data
- Filtering for membrane/surface proteins
- Building tissue-cell expression matrices
- Generating positive and negative labels from marker databases
"""

import csv
import re
import ast
import shutil
from pathlib import Path
from typing import Dict, List, Set, Tuple

import numpy as np
import pandas as pd

from constants import IMMUNE_SUBTYPES, NON_MEMBRANE_GENES
from io_utils import load_surface_proteins, save_labels_to_csv


def clean_expression_data(
    df: pd.DataFrame,
    surface_genes: Set[str],
    non_membrane_genes: Set[str] = NON_MEMBRANE_GENES
) -> pd.DataFrame:
    """
    Clean and filter the gene expression data.
    
    Applies the following filters:
    1. Remove "mixed" cell types
    2. Keep only surface genes
    3. Remove non-membrane genes
    4. Remove rows with NaN in Cell type or Gene name
    5. Remove PRF1 gene
    6. Set nTPM < 0.1 to 0
    7. Group immune cells into "immune" tissue
    
    Args:
        df: Raw gene expression DataFrame
        surface_genes: Set of validated surface protein genes
        non_membrane_genes: Set of genes to exclude as non-membrane
        
    Returns:
        Cleaned DataFrame
    """
    # 1. Remove "mixed" cell types
    df_clean = df[~df["Cell type"].str.contains("mixed", case=False, na=False)]
    
    # 2. Keep only surface genes
    df_clean = df_clean[df_clean['Gene name'].isin(surface_genes)]
    
    # 3. Remove non-membrane genes
    df_clean = df_clean[~df_clean['Gene name'].isin(non_membrane_genes)]
    
    # 4. Remove NaN in Cell type or Gene name
    df_clean = df_clean.dropna(subset=["Cell type", "Gene name"])
    
    # 5. Remove PRF1 gene
    df_clean = df_clean[df_clean["Gene name"] != "PRF1"]
    
    # 6. Set nTPM < 0.1 to 0
    df_clean = df_clean.copy()
    df_clean.loc[df_clean["nTPM"] < 0.1, "nTPM"] = 0
    
    # 7. Group immune cells into "immune" tissue
    def reassign_immune_subtype(row):
        for subtype in IMMUNE_SUBTYPES:
            if subtype == row["Cell type"]:
                return pd.Series(["immune", row["Cell type"]])
        return pd.Series([row["Tissue"], row["Cell type"]])
    
    df_clean[["Tissue", "Cell type"]] = df_clean.apply(reassign_immune_subtype, axis=1)
    
    return df_clean


def build_expression_matrix_high(df_clean: pd.DataFrame) -> Tuple[pd.DataFrame, List[str], List[List[str]]]:
    """
    Build expression matrix using highest nTPM from clusters.
    
    Args:
        df_clean: Cleaned expression DataFrame
        
    Returns:
        Tuple of (nTPM matrix DataFrame, gene list, tissue-cell pairs)
    """
    # Combine clusters by highest nTPM
    df_high = (
        df_clean.sort_values(
            ["Tissue", "Cell type", "Gene name", "nTPM"],
            ascending=[True, True, True, False]
        )
        .drop_duplicates(subset=["Tissue", "Cell type", "Gene name"], keep="first")
    )
    
    # Get sorted gene list
    gene_list = sorted(df_high["Gene name"].unique())
    
    # Create tissue-cell pairs
    tissue_cells = df_high[["Tissue", "Cell type"]].drop_duplicates().values.tolist()
    
    # Build nTPM matrix
    nTPM_matrix_df = (
        df_high
        .pivot_table(index=["Tissue", "Cell type"], columns="Gene name", values="nTPM", fill_value=0)
        .reindex(columns=gene_list)
        .reset_index()
    )
    
    return nTPM_matrix_df, gene_list, tissue_cells


def build_expression_matrix_median(df_clean: pd.DataFrame) -> Tuple[pd.DataFrame, List[str], List[List[str]]]:
    """
    Build expression matrix using median nTPM from clusters.
    
    Args:
        df_clean: Cleaned expression DataFrame
        
    Returns:
        Tuple of (nTPM matrix DataFrame, gene list, tissue-cell pairs)
    """
    # Combine clusters by median nTPM
    df_median = (
        df_clean
        .groupby(["Tissue", "Cell type", "Gene name"], as_index=False)["nTPM"]
        .median()
    )
    
    # Get sorted gene list
    gene_list = sorted(df_median["Gene name"].unique())
    
    # Create tissue-cell pairs
    tissue_cells = df_median[["Tissue", "Cell type"]].drop_duplicates().values.tolist()
    
    # Build nTPM matrix
    nTPM_matrix_df = (
        df_median
        .pivot_table(index=["Tissue", "Cell type"], columns="Gene name", values="nTPM", fill_value=0)
        .reindex(columns=gene_list)
        .reset_index()
    )
    
    return nTPM_matrix_df, gene_list, tissue_cells


def load_marker_matches(file_path: str) -> Tuple[Dict, Set]:
    """
    Load cell marker name matches from CSV file.
    
    The file format is a comma-separated file with two columns:
    - Column 1: Protein Atlas cell (tissue, cell type) as string representation of tuple
    - Column 2: Database cell (tissue, cell type) as string representation of tuple, or None
    
    The original notebook reads with tab delimiter to get each line as single element,
    then parses as Python literal.
    
    Args:
        file_path: Path to the name match CSV file
        
    Returns:
        Tuple of (matched_cells dict, name_no_match set)
    """
    matched_cells = {}
    name_no_match = set()
    
    with open(file_path, 'r', encoding='utf-8') as file:
        reader = csv.reader(file, delimiter='\t')
        next(reader, None)  # Skip header
        for row in reader:
            if not row:
                continue
            
            # Replace curly quotes with straight quotes and remove backslashes
            # Handle multiple types of curly quotes using unicode escape sequences only
            cleaned = row[0]
            # Left single quotation marks (U+2018, U+201B)
            cleaned = cleaned.replace('\u2018', "'")
            cleaned = cleaned.replace('\u201B', "'")
            # Right single quotation marks (U+2019)
            cleaned = cleaned.replace('\u2019', "'")
            # Remove backslashes
            cleaned = cleaned.replace('\\', '')
            
            try:
                formatted = ast.literal_eval(cleaned)
            except (SyntaxError, ValueError) as e:
                # Skip malformed rows
                continue
            
            # formatted should be a tuple of two string representations
            if not isinstance(formatted, tuple) or len(formatted) < 2:
                if isinstance(formatted, tuple) and len(formatted) == 1:
                    name_no_match.add(formatted[0])
                continue
                
            if formatted[1] is None or formatted[1] == 'None':
                name_no_match.add(formatted[0])
            else:
                try:
                    # Parse the inner tuple string
                    inner_tuple = ast.literal_eval(formatted[1])
                    matched_cells[formatted[0]] = inner_tuple
                except (SyntaxError, ValueError):
                    name_no_match.add(formatted[0])
    
    return matched_cells, name_no_match


def get_cellmarker_genes(
    matched_cells: Dict,
    cellmarker_df: pd.DataFrame
) -> Dict[Tuple[str, str], List[str]]:
    """
    Get marker genes from CellMarker database.
    
    Args:
        matched_cells: Dictionary mapping cell names to (tissue, cell) tuples
        cellmarker_df: CellMarker DataFrame
        
    Returns:
        Dictionary mapping (tissue, cell) to list of marker genes
    """
    result = {}
    for key in matched_cells:
        tissue, cell = matched_cells[key]
        genes = cellmarker_df[
            (cellmarker_df['species'] == 'Human') &
            (cellmarker_df['tissue_type'] == tissue) &
            (cellmarker_df['cell_name'] == cell)
        ]['Symbol'].tolist()
        result[key] = list(set(genes))
    return result


def get_panglao_genes(
    matched_cells: Dict,
    panglao_df: pd.DataFrame
) -> Dict[Tuple[str, str], List[str]]:
    """
    Get marker genes from PanglaoDB database.
    
    Args:
        matched_cells: Dictionary mapping cell names to (tissue, cell) tuples
        panglao_df: PanglaoDB DataFrame
        
    Returns:
        Dictionary mapping (tissue, cell) to list of marker genes
    """
    result = {}
    for key in matched_cells:
        tissue, cell = matched_cells[key]
        genes = panglao_df[
            (panglao_df['species'] != 'Mm') &
            (panglao_df['organ'] == tissue) &
            (panglao_df['cell type'] == cell)
        ]['official gene symbol'].tolist()
        result[key] = list(set(genes))
    return result


def combine_marker_databases(
    cellmarker_genes: Dict,
    panglao_genes: Dict
) -> Dict[Tuple[str, str], List[str]]:
    """
    Combine markers from CellMarker and PanglaoDB databases.
    
    Args:
        cellmarker_genes: Markers from CellMarker
        panglao_genes: Markers from PanglaoDB
        
    Returns:
        Combined dictionary of markers
    """
    markers = {}
    
    # First add all Panglao markers
    for key, genes in panglao_genes.items():
        if key in cellmarker_genes:
            combined = list(set(cellmarker_genes[key] + genes))
        else:
            combined = genes
        markers[key] = combined
    
    # Then add CellMarker markers not in Panglao
    for key, genes in cellmarker_genes.items():
        if key not in markers:
            markers[key] = genes
    
    return markers


def build_positive_labels(
    markers: Dict[Tuple[str, str], List[str]],
    gene_list: List[str],
    tissue_cells: List[List[str]]
) -> Dict[Tuple[str, str], List[str]]:
    """
    Build positive labels dictionary from marker databases.
    
    Args:
        markers: Combined marker dictionary
        gene_list: List of valid genes in the dataset
        tissue_cells: List of valid tissue-cell pairs
        
    Returns:
        Dictionary mapping (tissue, cell) to list of positive marker genes
    """
    positives = {}
    gene_set = set(gene_list)
    
    for key, genes in markers.items():
        tissue, cell = key
        
        # Filter to genes that exist in our gene list
        positive_genes = [gene for gene in genes if gene in gene_set]
        
        if not positive_genes:
            continue
        
        # Handle immune cell subtypes
        if cell in IMMUNE_SUBTYPES:
            positives[("immune", cell)] = positive_genes
        elif [tissue, cell] in tissue_cells:
            positives[(tissue, cell)] = positive_genes
    
    return positives


def build_negative_labels(
    positives: Dict[Tuple[str, str], List[str]],
    tissue_cells: List[List[str]]
) -> Dict[Tuple[str, str], List[str]]:
    """
    Build negative labels dictionary.
    
    Markers that are positive for one cell type become negative labels
    for other cell types in the same tissue.
    
    Args:
        positives: Positive labels dictionary
        tissue_cells: List of tissue-cell pairs
        
    Returns:
        Dictionary mapping (tissue, cell) to list of negative marker genes
    """
    negatives = {}
    
    # Build tissue to cells mapping
    tissue_to_cells = {}
    for tissue, cell in tissue_cells:
        if tissue not in tissue_to_cells:
            tissue_to_cells[tissue] = set()
        tissue_to_cells[tissue].add(cell)
    
    # Build negatives
    for (tissue, cell_type), positive_genes in positives.items():
        for other_cell_type in tissue_to_cells.get(tissue, []):
            if other_cell_type != cell_type:
                key = (tissue, other_cell_type)
                if key in negatives:
                    negatives[key] = list(set(negatives[key] + positive_genes))
                else:
                    negatives[key] = positive_genes.copy()
    
    return negatives


def run_data_processing(data_dir: str, output_dir: str) -> Dict:
    """
    Run the complete data processing pipeline.
    
    Args:
        data_dir: Directory containing input data files
        output_dir: Directory for output files
        
    Returns:
        Dictionary containing processed data for downstream analysis
    """
    data_path = Path(data_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    print("=== Data Processing Pipeline ===\n")
    
    # Load surface proteins
    print("Loading surface proteins...")
    surface_genes = load_surface_proteins(data_path / "mass_spec_valid_surface_protein.csv")
    print(f"  Loaded {len(surface_genes)} surface protein genes")
    
    # Load and clean expression data
    print("\nLoading and cleaning expression data...")
    df = pd.read_csv(data_path / "rna_single_cell_type_tissue.tsv", sep='\t')
    df_clean = clean_expression_data(df, surface_genes)
    print(f"  Cleaned data: {len(df_clean)} rows")
    
    # Build expression matrices
    print("\nBuilding expression matrices...")
    nTPM_high_df, gene_list_high, tissue_cells_high = build_expression_matrix_high(df_clean)
    nTPM_median_df, gene_list_median, tissue_cells_median = build_expression_matrix_median(df_clean)
    
    print(f"  High method matrix: {nTPM_high_df.shape}")
    print(f"  Median method matrix: {nTPM_median_df.shape}")
    
    # Use median method results as primary
    gene_list = gene_list_median
    tissue_cells = tissue_cells_median
    
    # Save expression matrices
    nTPM_high_df.to_csv(output_path / "gene_expression_matrix_high.csv", sep='\t', index=False)
    nTPM_median_df.to_csv(output_path / "gene_expression_matrix_median.csv", sep='\t', index=False)
    
    # Save tissue-cell pairs
    tissue_cell_df = pd.DataFrame(tissue_cells, columns=["Tissue", "Cell type"])
    tissue_cell_df.to_csv(output_path / "tissue_cell_pairs.tsv", sep='\t', index=False)
    
    # Save gene list
    gene_list_df = pd.DataFrame(gene_list, columns=["Gene name"])
    gene_list_df.to_csv(output_path / "gene_list.csv", index=False, header=False)
    
    print(f"\n  Saved {len(tissue_cells)} tissue-cell pairs")
    print(f"  Saved {len(gene_list)} genes")
    
    # Load marker databases
    print("\nLoading marker databases...")
    cellmarker_df = pd.read_excel(data_path / "Cell_marker_Human.xlsx")
    panglao_df = pd.read_csv(data_path / "PanglaoDB_markers_27_Mar_2020.tsv", delimiter='\t')
    
    # Load marker matches
    cellmarker_matches, _ = load_marker_matches(data_path / "CellMarker_name_match.csv")
    panglao_matches, _ = load_marker_matches(data_path / "PanglaoDB_name_match.csv")
    
    # Get markers from databases
    cellmarker_genes = get_cellmarker_genes(cellmarker_matches, cellmarker_df)
    panglao_genes = get_panglao_genes(panglao_matches, panglao_df)
    
    # Combine markers
    markers = combine_marker_databases(cellmarker_genes, panglao_genes)
    
    # Build labels
    print("\nBuilding positive and negative labels...")
    positives = build_positive_labels(markers, gene_list, tissue_cells)
    negatives = build_negative_labels(positives, tissue_cells)
    
    positive_count = sum(len(genes) for genes in positives.values())
    negative_count = sum(len(genes) for genes in negatives.values())
    
    print(f"  Positive labels: {len(positives)} cells, {positive_count} total labels")
    print(f"  Negative labels: {len(negatives)} cells, {negative_count} total labels")
    
    # Save labels
    save_labels_to_csv(positives, output_path / "positives_labels.csv", "Positive")
    save_labels_to_csv(negatives, output_path / "negative_labels.csv", "Negative")
    
    # Copy common cells file to output for controlled learning
    common_cells_src = data_path / "common_cells_across_tissues.csv"
    common_cells_dst = output_path / "common_cells_across_tissues.csv"
    shutil.copy(common_cells_src, common_cells_dst)
    print(f"  Copied common_cells_across_tissues.csv to output directory")
    
    print("\n=== Data Processing Complete ===\n")
    
    return {
        'nTPM_matrix_high_df': nTPM_high_df,
        'nTPM_matrix_median_df': nTPM_median_df,
        'gene_list': gene_list,
        'tissue_cells': tissue_cells,
        'positives': positives,
        'negatives': negatives
    }


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Process gene expression data")
    parser.add_argument("--data-dir", default=".", help="Directory containing input data")
    parser.add_argument("--output-dir", default=".", help="Directory for output files")
    
    args = parser.parse_args()
    run_data_processing(args.data_dir, args.output_dir)
