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
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional

import numpy as np
import pandas as pd


# Define immune cell subtypes that should be grouped together
IMMUNE_SUBTYPES = [
    "t-cells", "b-cells", "nk-cells", "monocytes", "macrophages",
    "granulocytes", "plasma cells", "dendritic cells"
]

# Non-membrane genes to exclude (extracted from UniProt data)
NON_MEMBRANE_GENES = {
    'YBX1', 'PPT1', 'C6orf120', 'RPL35', 'AEBP1', 'ANGPT1', 'PRG2', 'LAMB2',
    'LMAN2L', 'NAGPA', 'PLD3', 'IFI30', 'MAN2A2', 'NCLN', 'TOR1AIP2', 'SEMA4C',
    'IL17D', 'LAMB1', 'STS', 'RPL28', 'MAN2B1', 'ATP6AP1', 'CRELD2', 'BRI3BP',
    'ADAM15', 'ORM1', 'GBA', 'LIPA', 'SCARA3', 'SLC5A3', 'LTBP1', 'PXDN',
    'ST8SIA4', 'B3GNT2', 'HEPACAM2', 'NUP210', 'CTSC', 'PFN1', 'GALNS',
    'TUBB2C', 'KNG1', 'MOXD1', 'WNT5B', 'FKTN', 'WNT11', 'ELOVL4', 'POSTN',
    'GDF15', 'SNRPD2', 'WNT5A', 'SLC43A3', 'ABI3BP', 'GRN', 'EPX', 'SERPINI1',
    'PLD4', 'APOB', 'LMNA', 'FAM20B', 'CRB2', 'RPL35A', 'UBA1', 'LZTFL1',
    'LAMA5', 'RPL9', 'LTF', 'EMILIN1', 'COL14A1', 'PPT2', 'MIA3', 'ATRX',
    'HGSNAT', 'COL11A1', 'HSP90B1', 'APOH', 'BTD', 'TFPI2', 'VWF', 'NXPH4',
    'PLAU', 'FBN2', 'PCBP2', 'CFH', 'FKBP9', 'METTL9', 'HNRNPU', 'ATF6B',
    'NPTX1', 'C4BPB', 'SV2A', 'TUBB', 'PPIB', 'ITPR3', 'ANGPTL4', 'CHGA',
    'CD46', 'OGFOD1', 'DCN', 'GDPD5', 'LAMC2', 'F5', 'PRELP', 'ASAH1',
    'QSOX1', 'ISLR', 'GMPS', 'TGFB1', 'PLOD3', 'GMPPB', 'TMEM9', 'GLB1',
    'APOD', 'ITIH2', 'GAPDHS', 'LGMN', 'IGLON5', 'CASD1', 'SERPINB2',
    'SERPINF1', 'NXF2', 'OGN', 'DMXL2', 'FN1', 'BDNF', 'PRSS23', 'PLBD1',
    'RCN3', 'RPL18', 'CCDC80', 'ELANE', 'HEXA', 'ABCC5', 'PON1', 'HP',
    'IGF2R', 'SULF2', 'RPL21', 'TIMP1', 'YWHAQ', 'CCDC126', 'TOR1B',
    'PTGFRN', 'TPP1', 'COL6A3', 'HSPA13', 'ABCC4', 'ACTN4', 'RNF13', 'FGL2',
    'COL6A1', 'PCBP3', 'DKK3', 'FGB', 'BCHE', 'RSL1D1', 'NAGA', 'LEMD2',
    'CHSY3', 'SERPINA1', 'GALNT5', 'ADAMTS15', 'POFUT2', 'HSD11B1', 'COL3A1',
    'TOR3A', 'HMCN1', 'C4BPA', 'CALR', 'NOP56', 'GNPTAB', 'SERPINH1',
    'CASP7', 'THSD4', 'PLAT', 'VCAN', 'EDEM3', 'CGN', 'RPL17', 'PSAP',
    'MAN2B2', 'LMF2', 'TLL2', 'TLL1', 'GNS', 'CTSD', 'ERAP1', 'TYR', 'HYOU1',
    'PTPRQ', 'RPL3', 'HEXB', 'F2', 'ARSF', 'GNPTG', 'TF', 'COL12A1',
    'SLCO3A1', 'STT3B', 'OLFML2B', 'COL1A1', 'EXTL2', 'GZMH', 'FKBP14',
    'PCYOX1', 'M6PR', 'SEMA3C', 'TMX4', 'ARSB', 'COL6A6', 'LMAN2', 'PPIAL4A',
    'LRPAP1', 'DNASE1L1', 'HSPD1', 'SPANXB1', 'INTS12', 'APLF', 'TUBA4A',
    'THBS1', 'TNFRSF11B', 'B4GALT3', 'EFEMP1', 'CDSN', 'COL5A2', 'RPL32',
    'CHST3', 'LTBP3', 'WDR13', 'LAMA1', 'HNRNPK', 'FKBP7', 'MPO', 'MMP13',
    'ERGIC2', 'ANGPTL1', 'ACP2', 'SLC29A2', 'GLCE', 'CTSA', 'ADAMTSL4',
    'OLFML3', 'UST', 'TRAP1', 'CERCAM', 'LGALS3BP', 'SERPING1', 'NPTX2',
    'FGA', 'CTSL1', 'INHBE', 'IGFBP3', 'UGGT1', 'PBXIP1', 'B4GALT5', 'SEL1L',
    'GLA', 'HSD17B2', 'MAN2A1', 'ANGPTL2', 'SULF1', 'FLNA', 'LAMC1', 'CES1',
    'BIRC6', 'SUMF1', 'ACTB', 'FLG2', 'SLIT2', 'SLC22A4', 'RPL37A', 'SLIT3',
    'RPS23', 'LOXL4', 'PIGT', 'GPI', 'ENTPD7', 'SERINC1', 'LUM', 'OSTM1',
    'OLFML2A', 'LIPG', 'A2M', 'COL2A1', 'FBLN2', 'FUT11', 'TPI1', 'HPX',
    'C4A', 'SEMA3A', 'SFRS3', 'CCT4', 'SERPINE2', 'CHST7', 'TXNDC15', 'RALY',
    'LAMP3', 'PI16', 'FICD', 'SLCO4C1', 'UQCRC2', 'CREG1', 'EEF1D', 'RPL27A',
    'FSTL1', 'CRTAP', 'DPP7', 'NAV2', 'ZNF844', 'VEGFC', 'RCN1', 'TFPI',
    'DMTF1', 'ASPN', 'P4HA2', 'MPHOSPH9', 'XRCC6', 'NELL1', 'ASPH', 'NTN3',
    'TNC', 'RPS15A', 'AGPS', 'MYH9', 'PRSS35', 'COL5A1', 'ABCC2', 'SERPINA3',
    'VTN', 'TUBB6', 'DBH', 'SMPD1', 'CPZ', 'TMED9', 'CP', 'MXRA5', 'AKAP11',
    'MMRN1', 'ERAP2', 'COL7A1', 'NAGLU', 'PLTP', 'MFAP4', 'MGAT4B', 'PRKCSH',
    'C1QTNF2', 'ZC3HAV1', 'P4HA1', 'TNXB', 'EMILIN2', 'OLFM1', 'SIL1',
    'CRLF1', 'SCARB2', 'GREM1', 'PA2G4', 'RPL19', 'UGGT2', 'FKBP10',
    'COL5A3', 'IDH2', 'EPDR1', 'PLOD2', 'BGN', 'RERG', 'PLBD2', 'FBN1',
    'PTBP1', 'AHSG', 'ACAA1', 'PLOD1', 'GUSB', 'BMP1', 'LOX', 'MFAP5',
    'LAMB3', 'C3', 'ST6GALNAC3', 'COL1A2', 'F13A1'
}


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
    
    Args:
        file_path: Path to the name match CSV file
        
    Returns:
        Tuple of (matched_cells dict, name_no_match set)
    """
    matched_cells = {}
    name_no_match = set()
    
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        next(reader, None)  # Skip header
        for row in reader:
            formatted = ast.literal_eval(
                row[0].replace(''', "'").replace(''', "'").replace('\\', '')
            )
            if len(formatted) > 1:
                if formatted[1] is None:
                    name_no_match.add(formatted[0])
                else:
                    matched_cells[formatted[0]] = ast.literal_eval(formatted[1])
            else:
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


def save_labels_to_csv(
    labels: Dict[Tuple[str, str], List[str]],
    output_path: str,
    label_type: str = "Positive"
):
    """
    Save labels dictionary to CSV file.
    
    Args:
        labels: Labels dictionary
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
    surface_genes = load_surface_proteins(data_path / "Original data" / "mass_spec_valid_surface_protein.csv")
    print(f"  Loaded {len(surface_genes)} surface protein genes")
    
    # Load and clean expression data
    print("\nLoading and cleaning expression data...")
    df = pd.read_csv(data_path / "Original data" / "rna_single_cell_type_tissue.tsv", sep='\t')
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
    cellmarker_df = pd.read_excel(data_path / "Original data" / "Cell_marker_Human.xlsx")
    panglao_df = pd.read_csv(data_path / "Original data" / "PanglaoDB_markers_27_Mar_2020.tsv", delimiter='\t')
    
    # Load marker matches
    cellmarker_matches, _ = load_marker_matches(data_path / "Original data" / "CellMarker_name_match.csv")
    panglao_matches, _ = load_marker_matches(data_path / "Original data" / "PanglaoDB_name_match.csv")
    
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

