"""
Cell Surface Marker Prediction Pipeline

This package contains modules for predicting cell-specific surface markers
for targeted delivery applications.

Modules:
--------
constants
    Configuration constants (immune subtypes, non-membrane genes, literature positives)

io_utils
    Shared I/O utilities for loading and saving data files

data_processing
    Processing and cleaning gene expression data, building expression matrices,
    and generating positive/negative labels from marker databases

controlled_learning
    Penalty matrix optimization algorithm for marker prediction,
    including grid search and marker recommendations

plot_results
    Analysis plots and Excel reports (within-cell ranks, target/off-target ratios)

Usage:
------
Run the complete pipeline:
    $ python main.py

Or run individual modules:
    $ python -m src.data_processing --data-dir Data --output-dir Results
    $ python -m src.controlled_learning --data-dir Results --output-dir Results
    $ python -m src.plot_results --data-dir Results --plots-dir Plots
"""

__version__ = "1.0.0"

# Expose main functions for programmatic use
from .data_processing import run_data_processing
from .controlled_learning import run_controlled_learning
from .plot_results import run_plotting
