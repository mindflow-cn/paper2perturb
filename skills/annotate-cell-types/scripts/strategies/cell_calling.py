"""Lightweight cell calling for raw 10x MEX matrices.

DO NOT annotate millions of empty barcodes. Filter to likely real cells
using conservative UMI/gene-count knee detection before annotation.
DropletUtils/CellBender are deep-mode only.
"""

import json
import logging
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp

logger = logging.getLogger(__name__)


def _knee_threshold(
    umi_counts: np.ndarray,
    min_cells: int = 1000,
    max_cells: int = 200_000,
) -> int:
    """Conservative knee detection on sorted UMI counts.

    Returns the number of barcodes to keep (those above the knee).
    """
    sorted_counts = np.sort(umi_counts)[::-1]
    if len(sorted_counts) < 2:
        return len(sorted_counts)

    # Look for the inflection point: biggest drop in log-log slope
    log_rank = np.log10(np.arange(1, len(sorted_counts) + 1))
    log_counts = np.log10(np.maximum(sorted_counts, 1))

    # Smooth second derivative approach: find where slope changes most
    # Use a windowed approach
    window = max(10, len(sorted_counts) // 100)
    slopes = np.zeros(len(sorted_counts) - window)
    for i in range(len(slopes)):
        dy = log_counts[i + window] - log_counts[i]
        dx = log_rank[i + window] - log_rank[i]
        slopes[i] = dy / dx if dx > 0 else 0

    # Find where slope steepens most (second derivative min)
    slope_diff = np.diff(slopes)
    if len(slope_diff) > 0:
        knee_idx = np.argmin(slope_diff) + window // 2
    else:
        knee_idx = len(sorted_counts) // 2

    knee_idx = max(min_cells, min(knee_idx, max_cells))
    return int(knee_idx)


def call_cells(
    adata: "ad.AnnData",
    min_umis: int = 100,
    min_genes: int = 50,
    max_cells: int = 200_000,
) -> "ad.AnnData":
    """Filter raw MEX to likely real cells.

    Strategy:
    1. Compute UMI and gene counts per barcode
    2. Knee filter on UMI counts
    3. Conservative min UMI/gene thresholds
    """
    import anndata as ad

    n_barcodes = adata.n_obs
    logger.info(f"Cell calling: {n_barcodes} barcodes -> filtering")

    if n_barcodes <= max_cells:
        # Already reasonable — apply only min thresholds
        umi_sums = np.array(adata.X.sum(axis=1)).ravel()
        gene_counts = np.array((adata.X > 0).sum(axis=1)).ravel()
        keep = (umi_sums >= min_umis) & (gene_counts >= min_genes)
        logger.info(
            f"Min-threshold filter: {keep.sum()} / {n_barcodes} barcodes kept"
        )
        return adata[keep].copy()

    # Large raw matrix — knee filter first
    umi_sums = np.array(adata.X.sum(axis=1)).ravel()
    n_keep = _knee_threshold(umi_sums, max_cells=max_cells)

    # Get threshold UMI at knee
    sorted_umis = np.sort(umi_sums)[::-1]
    umi_threshold = sorted_umis[n_keep - 1] if n_keep <= len(sorted_umis) else 0

    gene_counts = np.array((adata.X > 0).sum(axis=1)).ravel()
    keep = (umi_sums >= max(umi_threshold, min_umis)) & (gene_counts >= min_genes)

    logger.info(
        f"Knee filter: {keep.sum()} / {n_barcodes} barcodes kept "
        f"(UMI threshold >= {max(umi_threshold, min_umis)})"
    )
    return adata[keep].copy()


def main():
    """CLI entry point for standalone cell calling.

    Usage:
        python3 cell_calling.py input.h5ad output.h5ad [--summary summary.json]
    """
    import anndata as ad

    if len(sys.argv) < 3:
        print(f"Usage: python3 {__file__} input.h5ad output.h5ad [--summary summary.json]")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]
    summary_path = None
    if len(sys.argv) > 3 and sys.argv[3] == "--summary" and len(sys.argv) > 4:
        summary_path = sys.argv[4]

    adata = ad.read_h5ad(input_path)
    n_before = adata.n_obs
    adata_filtered = call_cells(adata)
    n_after = adata_filtered.n_obs

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    adata_filtered.write(output_path, compression="gzip")

    summary = {
        "method": "knee_filter",
        "n_barcodes_before": int(n_before),
        "n_cells_after": int(n_after),
        "fraction_kept": float(n_after / n_before) if n_before > 0 else 0.0,
    }
    print(json.dumps(summary, indent=2))

    if summary_path:
        Path(summary_path).parent.mkdir(parents=True, exist_ok=True)
        with open(summary_path, "w") as f:
            json.dump(summary, f, indent=2)


if __name__ == "__main__":
    main()
