"""Metadata-based labeling for cell-line and in-vitro datasets.

For homogeneous datasets (cell lines, organoids, in-vitro differentiation),
label all cells with the matching target cell type from the Case-level table.
"""

import json
import logging
import sys
from pathlib import Path
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)


def label_from_metadata(
    adata: "ad.AnnData",
    target_cell_types: list[str],
    strategy_category: str = "cell_line_or_in_vitro",
) -> dict:
    """Label all cells with the best-matching target cell type.

    For cell_line_or_in_vitro, uses target_cell_types[0] as the homogeneous label.
    For other categories, sets 'Unknown' with annotation_method='Metadata'.

    Returns:
        dict with summary for manifest.
    """
    if "cell_type" not in adata.obs.columns:
        adata.obs["cell_type"] = None
    if "annotation_method" not in adata.obs.columns:
        adata.obs["annotation_method"] = "Unresolved"

    # Only label cells that are still missing
    mask = adata.obs["cell_type"].isna() | (
        adata.obs["cell_type"].astype(str).str.strip().str.lower().isin(
            ["", "nan", "none", "null"]
        )
    )

    if mask.sum() == 0:
        return {
            "method": "metadata_label",
            "cells_labeled": 0,
            "label_used": None,
        }

    label = _pick_best_label(target_cell_types, adata)

    adata.obs.loc[mask, "cell_type"] = label
    adata.obs.loc[mask, "annotation_method"] = "Metadata"

    logger.info(
        f"Metadata labeling: {mask.sum()} cells -> '{label}'"
    )

    return {
        "method": "metadata_label",
        "cells_labeled": int(mask.sum()),
        "label_used": label,
    }


def _pick_best_label(
    target_cell_types: list[str],
    adata: "ad.AnnData",
) -> str:
    """Pick the most specific target cell type label."""
    if not target_cell_types:
        return "Unknown"

    # Favour longer (more specific) labels
    sorted_types = sorted(target_cell_types, key=len, reverse=True)

    # If adata.obs has a sample or tissue column, try to refine
    for col in ["sample_system", "tissue", "source_type"]:
        if col in adata.obs.columns:
            # Check if any target matches the metadata
            meta_vals = adata.obs[col].dropna().astype(str).str.lower().unique()
            for ct in sorted_types:
                for mv in meta_vals:
                    if ct.lower() in mv or mv in ct.lower():
                        return ct

    return sorted_types[0]


def main():
    """CLI for standalone metadata labeling.

    Usage:
        python3 metadata_label.py input.h5ad output.h5ad \
            --target-cell-types '["CD8 T cells"]' \
            [--summary summary.json]
    """
    import anndata as ad

    if len(sys.argv) < 3:
        print(f"Usage: python3 {__file__} input.h5ad output.h5ad ...")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]

    target_cell_types = []
    summary_path = None
    i = 3
    while i < len(sys.argv):
        if sys.argv[i] == "--target-cell-types" and i + 1 < len(sys.argv):
            target_cell_types = json.loads(sys.argv[i + 1])
            i += 2
        elif sys.argv[i] == "--summary" and i + 1 < len(sys.argv):
            summary_path = sys.argv[i + 1]
            i += 2
        else:
            i += 1

    adata = ad.read_h5ad(input_path)
    summary = label_from_metadata(adata, target_cell_types)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    adata.write(output_path, compression="gzip")

    print(json.dumps(summary, indent=2))
    if summary_path:
        Path(summary_path).parent.mkdir(parents=True, exist_ok=True)
        with open(summary_path, "w") as f:
            json.dump(summary, f, indent=2)


if __name__ == "__main__":
    main()
