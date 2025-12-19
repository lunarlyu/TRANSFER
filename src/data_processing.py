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
import ast
import shutil
from pathlib import Path
from typing import Dict, List, Set, Tuple

import pandas as pd

from .constants import IMMUNE_SUBTYPES, NON_MEMBRANE_GENES
from .io_utils import load_surface_proteins, save_labels_to_csv


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
        .pivot_table(index=["Tissue", "Cell type"], columns="Gene name", values="nTPM", fill_value=0.0)
        .reindex(columns=gene_list)
        .reset_index()
    )
    
    # Ensure all values are floats (not integers) for consistent output format
    for col in gene_list:
        nTPM_matrix_df[col] = nTPM_matrix_df[col].astype(float)
    
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
        .pivot_table(index=["Tissue", "Cell type"], columns="Gene name", values="nTPM", fill_value=0.0)
        .reindex(columns=gene_list)
        .reset_index()
    )
    
    # Ensure all values are floats (not integers) for consistent output format
    for col in gene_list:
        nTPM_matrix_df[col] = nTPM_matrix_df[col].astype(float)
    
    return nTPM_matrix_df, gene_list, tissue_cells


def load_and_get_cellmarker_genes(
    match_file_path: str,
    cellmarker_df: pd.DataFrame
) -> Tuple[Dict[Tuple[str, str], List[str]], int, int]:
    """
    Load CellMarker name matches and get marker genes.
        
    Args:
        match_file_path: Path to CellMarker_name_match.csv
        cellmarker_df: CellMarker DataFrame
        
    Returns:
        Tuple of (matched_cells dict with genes, match count, marker count)
    """
    matched_cells = {}
    name_no_match = set()
    
    with open(match_file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        next(reader, None)
        for row in reader:
            cleaned = row[0].replace('\u2018', "'").replace('\u2019', "'").replace('\\', '')
            formatted = ast.literal_eval(cleaned)
            if len(formatted) > 1:
                if formatted[1] == None:
                    name_no_match.add(formatted[0])
                else:
                    matched_cells[formatted[0]] = ast.literal_eval(formatted[1])
            else:
                name_no_match.add(formatted[0])
    
    # Get genes for each matched cell (overwrites matched_cells with gene lists)
    marker_count = 0
    for key in matched_cells:
        tissue, cell = matched_cells[key]
        genes = cellmarker_df[
            (cellmarker_df['species'] == 'Human') &
            (cellmarker_df['tissue_type'] == tissue) &
            (cellmarker_df['cell_name'] == cell)
        ]['Symbol'].tolist()
        genes = list(set(genes))
        marker_count += len(genes)
        matched_cells[key] = genes
    
    return matched_cells, len(matched_cells), marker_count


def load_and_get_panglao_genes(
    match_file_path: str,
    panglao_df: pd.DataFrame
) -> Tuple[Dict[Tuple[str, str], List[str]], int, int]:
    """
    Load PanglaoDB name matches and get marker genes.
    
    Note: CSV has leading spaces in some entries that cause ast.literal_eval to fail,
    so we strip whitespace before parsing.
    
    Args:
        match_file_path: Path to PanglaoDB_name_match.csv
        panglao_df: PanglaoDB DataFrame
        
    Returns:
        Tuple of (matched_cells dict with genes, match count, marker count)
    """
    matched_cells = {}
    name_no_match = set()
    
    with open(match_file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        next(reader, None)
        for row in reader:
            cleaned = row[0].replace('\u2018', "'").replace('\u2019', "'").replace('\\', '')
            formatted = ast.literal_eval(cleaned)
            if formatted[1] == None:
                name_no_match.add(formatted[0])
            else:
                # Strip whitespace to handle leading spaces in the CSV
                inner_str = formatted[1].strip() if isinstance(formatted[1], str) else formatted[1]
                matched_cells[formatted[0]] = ast.literal_eval(inner_str)
    
    # Get genes for each matched cell (overwrites matched_cells with gene lists)
    marker_count = 0
    for key in matched_cells:
        tissue, cell = matched_cells[key]
        genes = panglao_df[
            (panglao_df['species'] != 'Mm') &
            (panglao_df['organ'] == tissue) &
            (panglao_df['cell type'] == cell)
        ]['official gene symbol'].tolist()
        genes = list(set(genes))
        marker_count += len(genes)
        matched_cells[key] = genes
    
    return matched_cells, len(matched_cells), marker_count


def combine_marker_databases(
    cellmarker_genes: Dict[Tuple[str, str], List[str]],
    panglao_genes: Dict[Tuple[str, str], List[str]]
) -> Tuple[Dict[Tuple[str, str], List[str]], int]:
    """
    Combine markers from CellMarker and PanglaoDB databases.
        
    Args:
        cellmarker_genes: Markers from CellMarker (matched_cells1)
        panglao_genes: Markers from PanglaoDB (matched_cells2)
        
    Returns:
        Tuple of (combined markers dict, total marker count)
    """
    total_count = 0
    markers = {}
    
    # First iterate through PanglaoDB (matched_cells2)
    for key in panglao_genes:
        if key in cellmarker_genes:
            genes = cellmarker_genes[key] + panglao_genes[key]
            genes = list(set(genes))
        else:
            genes = panglao_genes[key]
        markers[key] = genes
        total_count += len(genes)
    
    # Then add CellMarker entries not in PanglaoDB
    for key in cellmarker_genes:
        if key not in markers:
            genes = cellmarker_genes[key]
            markers[key] = genes
            total_count += len(genes)
    
    return markers, total_count


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
    
    for key, genes in markers.items():
        tissue, cell = key
        
        # Filter to genes that exist in our gene list
        positive_genes = [gene for gene in genes if gene in gene_list]
        
        if not positive_genes:
            continue  # Skip if no matching genes
        
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
    
    cellmarker_genes, cm_match_count, cm_marker_count = load_and_get_cellmarker_genes(
        data_path / "CellMarker_name_match.csv", cellmarker_df
    )
    panglao_genes, pg_match_count, pg_marker_count = load_and_get_panglao_genes(
        data_path / "PanglaoDB_name_match.csv", panglao_df
    )
    
    print(f"  CellMarker: {cm_match_count} cells matched, {cm_marker_count} markers")
    print(f"  PanglaoDB: {pg_match_count} cells matched, {pg_marker_count} markers")
    
    markers, total_markers = combine_marker_databases(cellmarker_genes, panglao_genes)
    print(f"  Combined: {len(markers)} cells, {total_markers} markers")
    
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
