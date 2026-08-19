"""Detect and standardize original cell annotations.

Candidate annotation columns in adata.obs:
  - cell_type, celltype, cell.types
  - cell_type_fine, cell_type_major
  - annotation, cell_annotation
  - cell_ontology_class
  - seurat_annotations
  - predicted.celltype
  - major_celltype, minor_celltype
  - subtype, labels

Never treat as cell type:
  - cluster, cluster_id, leiden, seurat_clusters
  - sample, sample_id
  - condition, treatment, dose, time
"""

import logging
from typing import Optional

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)

CANDIDATE_COLUMNS = [
    "cell_type",
    "celltype",
    "cell.types",
    "cell_type_fine",
    "cell_type_major",
    "annotation",
    "cell_annotation",
    "cell_ontology_class",
    "seurat_annotations",
    "predicted.celltype",
    "major_celltype",
    "minor_celltype",
    "subtype",
    "labels",
]

BLACKLIST_COLUMNS = [
    "cluster",
    "cluster_id",
    "leiden",
    "seurat_clusters",
    "sample",
    "sample_id",
    "condition",
    "treatment",
    "dose",
    "time",
]


def _is_biological_labels(series: pd.Series) -> bool:
    """Check if a series contains biological-looking labels.

    Returns True if labels look like cell type names (not cluster numbers, etc.).
    """
    vals = series.dropna().astype(str)
    if len(vals) == 0:
        return False

    unique = vals.unique()

    # Reject purely numeric clusters
    if all(v.replace(".", "").replace("-", "").isdigit() for v in unique):
        return False

    # Reject short numeric-prefix clusters like "0", "1", "2" ...
    if len(unique) <= 3 and all(
        v.strip().isdigit() or v.strip() == "" for v in unique
    ):
        return False

    # Favour strings with spaces or biological words
    has_spaces = sum(1 for v in unique if " " in v)
    has_biological = sum(
        1 for v in unique
        if any(
            kw in v.lower()
            for kw in [
                "cell", "cyte", "blast", "phil", "phage",
                "epithelial", "endothelial", "fibroblast",
                "neuron", "keratinocyte", "immune", "t cell",
                "b cell", "nk ", "monocyte", "macrophage",
                "dendritic", "plasma", "stem", "progenitor",
            ]
        )
    )

    if has_spaces > 0 or has_biological > 0:
        return True

    # Accept short non-numeric codes (e.g. "ME", "BR", "HN" from GEO metadata)
    # as valid labels when all values are alphanumeric abbreviations
    has_abbrev = all(
        v.replace("_", "").replace("-", "").replace(".", "").isalnum()
        for v in unique
    )
    return has_abbrev and len(unique) >= 1


def detect_original_annotation(adata: "ad.AnnData") -> Optional[str]:
    """Find the best original annotation column in adata.obs.

    Selection rules:
    1. Choose candidate with biological-looking labels
    2. Prefer high non-null coverage
    3. Prefer specific column names over generic ones
    4. Prefer higher label diversity (not just 2-3 categories)

    Returns:
        Column name if found, else None.
    """
    obs_cols = set(c.lower() for c in adata.obs.columns)

    # Check each candidate, in priority order
    candidates_found = []
    for col in CANDIDATE_COLUMNS:
        # Match case-insensitively
        match = None
        for obs_col in adata.obs.columns:
            if obs_col.lower() == col.lower():
                match = obs_col
                break
        if match is None:
            continue

        # Must not be in blacklist
        if match.lower() in BLACKLIST_COLUMNS:
            continue

        series = adata.obs[match]
        coverage = series.notna().mean()
        n_unique = series.dropna().nunique()

        if not _is_biological_labels(series):
            logger.debug(f"Skipping {match}: doesn't look like biological labels")
            continue

        # Score: higher is better
        # Coverage matters most, then diversity, then name specificity
        name_score = 2 if match.lower() in ("cell_type", "celltype", "cell.types") else 0
        name_score += 1 if match.lower() in ("cell_type_fine", "cell_type_major") else 0

        candidates_found.append({
            "column": match,
            "coverage": coverage,
            "n_unique": n_unique,
            "name_score": name_score,
        })

    if not candidates_found:
        return None

    # Sort by coverage desc, then name_score desc, then n_unique desc
    candidates_found.sort(
        key=lambda x: (x["coverage"], x["name_score"], x["n_unique"]),
        reverse=True,
    )

    best = candidates_found[0]
    logger.info(
        f"Detected original annotation: {best['column']} "
        f"(coverage={best['coverage']:.1%}, n_unique={best['n_unique']})"
    )
    return best["column"]


def standardize_original_annotation(
    adata: "ad.AnnData",
    source_col: str,
) -> "ad.AnnData":
    """Standardize original annotation to obs['cell_type'] and obs['annotation_method'].

    - Copies values from source_col to obs['cell_type']
    - Sets obs['annotation_method'] = 'Original' for cells with non-null labels
    - Does NOT overwrite existing cell_type if already set
    """
    if "cell_type" not in adata.obs.columns:
        adata.obs["cell_type"] = adata.obs[source_col].astype(str)
        adata.obs.loc[
            adata.obs[source_col].isna(), "cell_type"
        ] = None
    else:
        # Fill only missing
        mask = adata.obs["cell_type"].isna() | (
            adata.obs["cell_type"].astype(str).str.strip().isin(["", "nan", "None", "none"])
        )
        adata.obs.loc[mask, "cell_type"] = adata.obs.loc[
            mask, source_col
        ].astype(str)

    if "annotation_method" not in adata.obs.columns:
        adata.obs["annotation_method"] = "Unresolved"
    # Set Original for cells that have a non-null, non-empty label from source
    has_label = adata.obs["cell_type"].notna() & (
        adata.obs["cell_type"].astype(str).str.strip().str.lower().isin(
            ["", "nan", "none", "null"]
        )
        == False  # noqa: E712
    )
    adata.obs.loc[has_label, "annotation_method"] = "Original"

    logger.info(
        f"Standardized original annotation: "
        f"{(adata.obs['annotation_method'] == 'Original').sum()} cells "
        f"marked as Original"
    )

    return adata
