#!/usr/bin/env python3
"""
Cell Surface Marker Prediction Pipeline

This script runs the complete pipeline for predicting cell-specific
surface markers for targeted delivery applications.

The pipeline consists of two main stages:
1. Data Processing: Clean and prepare gene expression data
2. Controlled Learning: Optimize penalty parameters and predict markers

Usage:
    python main.py [--data-dir DATA_DIR] [--output-dir OUTPUT_DIR] [--skip-processing] [--skip-plots]

Arguments:
    --data-dir       Directory containing input data (default: Data)
    --output-dir     Directory for output files (default: Results)
    --plots-dir      Directory for plot outputs (default: plots)
    --skip-processing  Skip data processing step (use existing processed files)
    --skip-plots     Skip plot generation step
"""

import argparse
import sys
from pathlib import Path

# Add src directory to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from data_processing import run_data_processing
from controlled_learning import run_controlled_learning
from plot_results import run_plotting


def main():
    """Run the complete cell surface marker prediction pipeline."""
    
    parser = argparse.ArgumentParser(
        description="Cell Surface Marker Prediction Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Run full pipeline with default directories
    python main.py

    # Run with custom data and output directories
    python main.py --data-dir /path/to/data --output-dir /path/to/output

    # Skip data processing (use existing processed files)
    python main.py --skip-processing
        """
    )
    
    parser.add_argument(
        "--data-dir",
        default="Data",
        help="Directory containing input data files (default: Data)"
    )
    
    parser.add_argument(
        "--output-dir",
        default="Results",
        help="Directory for output files (default: Results)"
    )
    
    parser.add_argument(
        "--skip-processing",
        action="store_true",
        help="Skip data processing step and use existing processed files"
    )
    
    parser.add_argument(
        "--plots-dir",
        default="plots",
        help="Directory for plot outputs (default: plots)"
    )
    
    parser.add_argument(
        "--skip-plots",
        action="store_true",
        help="Skip plot generation step"
    )
    
    args = parser.parse_args()
    
    # Resolve paths
    data_dir = Path(args.data_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    plots_dir = Path(args.plots_dir).resolve()
    
    print("=" * 60)
    print("  Cell Surface Marker Prediction Pipeline")
    print("=" * 60)
    print(f"\nData directory: {data_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Plots directory: {plots_dir}")
    print()
    
    # Validate data directory
    required_files = [
        "mass_spec_valid_surface_protein.csv",
        "rna_single_cell_type_tissue.tsv",
        "Cell_marker_Human.xlsx",
        "PanglaoDB_markers_27_Mar_2020.tsv",
        "CellMarker_name_match.csv",
        "PanglaoDB_name_match.csv",
        "common_cells_across_tissues.csv"
    ]
    
    if not args.skip_processing:
        missing_files = []
        for f in required_files:
            if not (data_dir / f).exists():
                missing_files.append(f)
        
        if missing_files:
            print("ERROR: Missing required input files:")
            for f in missing_files:
                print(f"  - {f}")
            print("\nPlease ensure all required files are in the data directory.")
            sys.exit(1)
    
    # Step 1: Data Processing
    if not args.skip_processing:
        print("\n" + "=" * 60)
        print("  STEP 1: Data Processing")
        print("=" * 60 + "\n")
        
        try:
            processing_results = run_data_processing(str(data_dir), str(output_dir))
            print(f"Data processing completed successfully.")
        except Exception as e:
            print(f"ERROR in data processing: {e}")
            sys.exit(1)
    else:
        print("\n[Skipping data processing - using existing files]\n")
        
        # Validate that processed files exist
        processed_files = [
            "gene_expression_matrix_high.csv",
            "gene_expression_matrix_median.csv",
            "tissue_cell_pairs.tsv",
            "gene_list.csv",
            "positives_labels.csv",
            "negative_labels.csv"
        ]
        
        missing = []
        for f in processed_files:
            if not (output_dir / f).exists():
                missing.append(f)
        
        if missing:
            print("ERROR: Missing processed files (run without --skip-processing first):")
            for f in missing:
                print(f"  - {f}")
            sys.exit(1)
    
    # Step 2: Controlled Learning
    print("\n" + "=" * 60)
    print("  STEP 2: Controlled Learning")
    print("=" * 60 + "\n")
    
    try:
        learning_results = run_controlled_learning(str(output_dir), str(output_dir))
        print(f"Controlled learning completed successfully.")
    except Exception as e:
        print(f"ERROR in controlled learning: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # Step 3: Generate Plots
    if not args.skip_plots:
        print("\n" + "=" * 60)
        print("  STEP 3: Generate Plots and Analysis")
        print("=" * 60 + "\n")
        
        try:
            run_plotting(
                data_dir=str(output_dir),
                output_dir=str(output_dir),
                plots_dir=str(plots_dir)
            )
            print(f"Plot generation completed successfully.")
        except Exception as e:
            print(f"ERROR in plot generation: {e}")
            import traceback
            traceback.print_exc()
            # Don't exit - plots are optional
    else:
        print("\n[Skipping plot generation]\n")
    
    # Summary
    print("\n" + "=" * 60)
    print("  Pipeline Complete")
    print("=" * 60)
    print("\nOutput files generated:")
    
    output_files = [
        "gene_expression_matrix_high.csv",
        "gene_expression_matrix_median.csv",
        "tissue_cell_pairs.tsv",
        "gene_list.csv",
        "positives_labels.csv",
        "negative_labels.csv",
        "recommended_whole_body_markers_high.csv",
        "recommended_whole_body_markers_median.csv"
    ]
    
    for f in output_files:
        path = output_dir / f
        if path.exists():
            size = path.stat().st_size / 1024
            print(f"  ✓ {f} ({size:.1f} KB)")
        else:
            print(f"  ✗ {f} (not created)")
    
    # List plot files if generated
    if not args.skip_plots and plots_dir.exists():
        print(f"\nPlot files (in {plots_dir.name}/):")
        plot_files = list(plots_dir.glob("*"))
        for f in sorted(plot_files):
            if f.is_file():
                size = f.stat().st_size / 1024
                print(f"  ✓ {f.name} ({size:.1f} KB)")
    
    print("\n" + "=" * 60)
    print("  Results Summary")
    print("=" * 60)
    
    if 'high' in learning_results:
        params = learning_results['high']['optimal_params']
        print(f"\nHIGH method optimal parameters:")
        print(f"  p = {params[0]:.2f}")
        print(f"  q = {params[1]:.2f}")
        print(f"  r = {params[2]:.2f}")
    
    if 'median' in learning_results:
        params = learning_results['median']['optimal_params']
        print(f"\nMEDIAN method optimal parameters:")
        print(f"  p = {params[0]:.2f}")
        print(f"  q = {params[1]:.2f}")
        print(f"  r = {params[2]:.2f}")
    
    print("\nPipeline execution completed successfully!")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

