# Cell Surface Marker Prediction Pipeline

A computational pipeline for predicting cell-specific surface markers for targeted delivery applications. This tool analyzes single-cell gene expression data to identify optimal cell surface markers that can distinguish specific cell types across tissues.

## 1. System Requirements

### Software Dependencies

- **Python**: >= 3.8
- **numpy**: >= 1.21.0
- **pandas**: >= 1.3.0
- **openpyxl**: >= 3.0.0 (for reading/writing Excel files)
- **matplotlib**: >= 3.4.0 (for generating plots)

### Operating Systems

The software has been tested on:

- macOS 14.x (Sonoma)
- macOS 13.x (Ventura)
- Ubuntu 20.04 LTS
- Ubuntu 22.04 LTS
- Windows 10/11 (with Python installed)

### Hardware Requirements

No non-standard hardware is required. The software runs on any standard desktop computer.

**Recommended specifications:**

- RAM: 8 GB minimum (16 GB recommended for larger datasets)
- CPU: Any modern multi-core processor
- Disk space: 500 MB for software and demo data

## 2. Installation Guide

### Instructions

1. **Clone or download the repository:**

   ```bash
   git clone <repository-url>
   ```

2. **Create a virtual environment (recommended):**

   ```bash
   python3 -m venv venv
   source venv/bin/activate  # On macOS/Linux
   # or
   venv\Scripts\activate  # On Windows
   ```

3. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

### Typical Install Time

Installation typically takes **1-2 minutes** on a standard desktop computer with a stable internet connection.

## 3. Demo

### Demo Dataset

The repository includes a complete demo dataset in the `Data/` directory containing:

| File                                  | Description                                              |
| ------------------------------------- | -------------------------------------------------------- |
| `rna_single_cell_type_tissue.tsv`     | Single-cell RNA expression data from Human Protein Atlas |
| `mass_spec_valid_surface_protein.csv` | Mass spectrometry-validated surface proteins             |
| `Cell_marker_Human.xlsx`              | CellMarker database annotations                          |
| `PanglaoDB_markers_27_Mar_2020.tsv`   | PanglaoDB marker annotations                             |
| `CellMarker_name_match.csv`           | Cell type name mappings for CellMarker                   |
| `PanglaoDB_name_match.csv`            | Cell type name mappings for PanglaoDB                    |
| `common_cells_across_tissues.csv`     | Cell types found across multiple tissues                 |

### Instructions to Run Demo

Run the complete pipeline with default settings:

```bash
python main.py
```

Or with explicit parameters:

```bash
python main.py --data-dir Data --output-dir Results --plots-dir Plots
```

### Expected Output

The pipeline generates the following output files in the `Results/` directory:

| File                                        | Description                            |
| ------------------------------------------- | -------------------------------------- |
| `gene_expression_matrix_high.csv`           | Expression matrix using HIGH method    |
| `gene_expression_matrix_median.csv`         | Expression matrix using MEDIAN method  |
| `tissue_cell_pairs.tsv`                     | List of tissue-cell type combinations  |
| `gene_list.csv`                             | List of surface protein genes analyzed |
| `positives_labels.csv`                      | Known positive marker labels           |
| `negative_labels.csv`                       | Negative marker labels                 |
| `recommended_whole_body_markers_high.csv`   | **Predicted markers (HIGH method)**    |
| `recommended_whole_body_markers_median.csv` | **Predicted markers (MEDIAN method)**  |

The pipeline also generates analysis plots in the `Plots/` directory:

| File                          | Description                             |
| ----------------------------- | --------------------------------------- |
| `topmarkers_withinCell_*.svg` | Histogram of marker expression ranks    |
| `top1marker_*.xlsx`           | Target/off-target ratios for top marker |
| `top10marker_*.xlsx`          | Best ratio among top 10 markers         |
| `top2_10marker_*.xlsx`        | Best 2-marker combination ratios        |

### Expected Console Output

```
============================================================
  Cell Surface Marker Prediction Pipeline
============================================================

Data directory: /path/to/Data
Output directory: /path/to/Results
Plots directory: /path/to/plots

============================================================
  STEP 1: Data Processing
============================================================

=== Data Processing Pipeline ===

Loading surface proteins...
  Loaded 1492 surface protein genes

Loading and cleaning expression data...
  Cleaned data: 8247 rows

Building expression matrices...
  High method matrix: (155, 1101)
  Median method matrix: (155, 1101)

  Saved 155 tissue-cell pairs
  Saved 1099 genes

Loading marker databases...

Building positive and negative labels...
  Positive labels: 155 cells, 1847 total labels
  Negative labels: 155 cells, 8174 total labels

=== Data Processing Complete ===

============================================================
  STEP 2: Controlled Learning
============================================================

Running for HIGH method...
  Running grid search over 8000 parameter combinations...
  Optimal parameters found: p=-700.00, q=-700.00, r=-20.00

Running for MEDIAN method...
  Running grid search over 8000 parameter combinations...
  Optimal parameters found: p=-700.00, q=-700.00, r=-20.00

============================================================
  Pipeline Complete
============================================================

Pipeline execution completed successfully!
```

### Expected Run Time

| Step                              | Time              |
| --------------------------------- | ----------------- |
| Data Processing                   | ~30 seconds       |
| Controlled Learning (grid search) | ~5-10 minutes     |
| Plot Generation                   | ~30 seconds       |
| **Total**                         | **~6-11 minutes** |

_Times measured on a standard desktop (Intel i7, 16 GB RAM)_

## 4. Instructions for Use

### Running on Your Own Data

To run the pipeline on your own data, prepare the following input files in your data directory:

#### Required Input Files

1. **`rna_single_cell_type_tissue.tsv`** - Gene expression data with columns:

   - `Gene name`: Gene symbol
   - `Tissue`: Tissue type
   - `Cell type`: Cell type name
   - `nTPM`: Normalized TPM expression value
   - `pTPM`: Percent TPM (optional)

2. **`mass_spec_valid_surface_protein.csv`** - Surface protein list with column:

   - `Gene name`: Gene symbols of validated surface proteins

3. **`Cell_marker_Human.xlsx`** - CellMarker database file (download from [CellMarker](http://xteam.xbio.top/CellMarker/))

4. **`PanglaoDB_markers_27_Mar_2020.tsv`** - PanglaoDB marker file (download from [PanglaoDB](https://panglaodb.se/))

5. **`CellMarker_name_match.csv`** - Mapping between your cell type names and CellMarker names

6. **`PanglaoDB_name_match.csv`** - Mapping between your cell type names and PanglaoDB names

7. **`common_cells_across_tissues.csv`** - List of cell types present in multiple tissues

### Command-Line Options

```bash
python main.py [OPTIONS]

Options:
  --data-dir PATH       Directory containing input data (default: Data)
  --output-dir PATH     Directory for output files (default: Results)
  --plots-dir PATH      Directory for plot outputs (default: plots)
  --skip-processing     Skip data processing, use existing processed files
  --skip-plots          Skip plot generation step
```

### Examples

```bash
# Run with custom directories
python main.py --data-dir /path/to/mydata --output-dir /path/to/results

# Skip data processing (reuse existing processed files)
python main.py --skip-processing

# Run without generating plots
python main.py --skip-plots

# Run individual modules
python -m src.data_processing --data-dir Data --output-dir Results
python -m src.controlled_learning --data-dir Results --output-dir Results
python -m src.plot_results --data-dir Results --plots-dir plots
```

### Output Interpretation

The main output files are `recommended_whole_body_markers_*.csv`, which contain:

- **Column 1**: Cell type (format: "Tissue CellType")
- **Column 2**: Ranked list of recommended surface markers

Higher-ranked markers have better specificity for distinguishing the target cell type from other cell types across all tissues.

## License

[Add your license information here]

## Citation

[Add citation information here]

## Contact

[Add contact information here]
