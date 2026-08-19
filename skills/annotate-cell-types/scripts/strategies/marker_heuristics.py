"""Marker-based heuristics for fast coarse cell-type annotation.

Two modes:
  1. Fill: annotate cells missing cell_type labels.
  2. Refine: cluster-level refinement of broad CellTypist labels into
     specific target cell types using target-specific markers.

Refinement operates on CLUSTERS (CellTypist label groups), not per-cell,
to avoid O(n_cells) scanning.
"""

import json
import logging
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp

logger = logging.getLogger(__name__)

# Detailed marker sets organized by lineage
MARKER_SETS = {
    # --- Immune ---
    "T cells": ["CD3D", "CD3E", "CD3G"],
    "CD8 T cells": ["CD8A", "CD8B"],
    "CD4 T cells": ["CD4"],
    "TEMRA T cells": ["CD8A", "GZMB", "NKG7", "PRF1", "KLRG1", "B3GAT1"],
    "NK cells": ["NKG7", "GNLY", "KLRD1", "KLRF1", "PRF1"],
    "B cells": ["CD19", "CD79A", "CD79B", "MS4A1"],
    "Plasma cells": ["MZB1", "SDC1", "IGHG1", "JCHAIN", "XBP1"],
    "Monocytes / Macrophages": ["CD14", "CD68", "FCGR3A", "CSF1R", "LYZ"],
    "Dendritic cells": ["FCER1A", "CLEC10A", "CLEC9A", "CD1C", "XCR1"],
    "DC3 / myeloid DC": ["LST1", "FCER1A", "CD1C", "CLEC10A", "LAMP3", "CCR7"],
    "Neutrophils": ["CSF3R", "FCGR3B", "CXCR2", "ELANE", "MPO"],

    # --- Skin / Epithelial ---
    "Keratinocytes": ["KRT5", "KRT14", "KRT1", "KRT10"],
    "Basal keratinocytes": ["KRT5", "KRT14", "COL17A1", "ITGA6", "ITGB4"],
    "Spinous keratinocytes": ["KRT1", "KRT10", "DSG1", "DSC1", "IVL"],
    "Supra-spinous / suprabasal keratinocytes": ["KRT1", "KRT10", "FLG", "LOR", "DSG1"],
    "Proliferating cells": ["MKI67", "TOP2A", "PCNA", "CCNB1", "CDK1"],
    "Proliferating keratinocytes": [
        "KRT5", "KRT14", "MKI67", "TOP2A", "PCNA",
    ],
    "Melanocytes": ["MLANA", "PMEL", "TYR", "DCT", "TYRP1"],

    # --- Endothelial / Vascular ---
    "Endothelial cells": ["PECAM1", "CDH5", "VWF", "CLDN5", "ENG"],
    "Capillary endothelial": ["PLVAP", "RGCC", "CA4", "AQP1"],
    "PLVAP+ capillaries": ["PLVAP", "RGCC", "CA4", "PECAM1"],
    "Lymphatic endothelial": ["PROX1", "LYVE1", "PDPN", "CCL21"],

    # --- Mesenchymal ---
    "Fibroblasts": ["COL1A1", "COL1A2", "DCN", "LUM", "FAP", "PDGFRA"],
    "WNT5A+ / IL24+ fibroblasts": ["WNT5A", "IL24", "COL1A1", "DCN"],
    "Myofibroblasts": ["ACTA2", "MYH11", "TAGLN"],
    "Pericytes": ["RGS5", "PDGFRB", "CSPG4", "MCAM", "CD248", "ACTA2"],

    # --- Epithelial / Organ ---
    "Epithelial cells": ["EPCAM", "KRT8", "KRT18", "CDH1"],
    "Hepatocytes": ["ALB", "CYP3A4", "ASGR1", "HNF4A"],
    "Neurons": ["RBFOX3", "TUBB3", "MAP2", "SYN1", "SYP"],
    "Glial cells": ["GFAP", "OLIG2", "MBP", "SOX10", "AQP4"],
}

# Map target keywords to relevant marker set names so we prioritise them
TARGET_TO_MARKER_PRIORITY = {
    "keratinocyte": [
        "Basal keratinocytes", "Spinous keratinocytes",
        "Supra-spinous / suprabasal keratinocytes",
        "Proliferating keratinocytes", "Keratinocytes",
    ],
    "t cell": [
        "T cells", "CD8 T cells", "CD4 T cells", "TEMRA T cells",
    ],
    "endothelial": [
        "Endothelial cells", "Capillary endothelial",
        "PLVAP+ capillaries", "Lymphatic endothelial",
    ],
    "capillary": [
        "Capillary endothelial", "PLVAP+ capillaries", "Endothelial cells",
    ],
    "fibroblast": [
        "Fibroblasts", "WNT5A+ / IL24+ fibroblasts", "Myofibroblasts",
    ],
    "pericyte": ["Pericytes"],
    "dendritic": [
        "Dendritic cells", "DC3 / myeloid DC",
    ],
    "b cell": ["B cells", "Plasma cells"],
    "macrophage": ["Monocytes / Macrophages"],
    "monocyte": ["Monocytes / Macrophages"],
    "neutrophil": ["Neutrophils"],
    "nk": ["NK cells"],
    "melanocyte": ["Melanocytes"],
    "epithelial": ["Epithelial cells"],
    "neuron": ["Neurons"],
    "glial": ["Glial cells"],
    "hepatocyte": ["Hepatocytes"],
}


def _normalize(s: str) -> str:
    return str(s).strip().lower().replace("_", " ").replace("-", " ")


def _normalize_gene(name: str) -> str:
    return str(name).strip().upper()


def _get_prioritized_marker_sets(target_cell_types: list[str]) -> list[tuple[str, list[str]]]:
    """Return marker sets ordered by relevance to target cell types."""
    if not target_cell_types:
        return list(MARKER_SETS.items())

    all_target_text = " ".join(t.lower() for t in target_cell_types)

    priority_sets = []
    for kw, set_names in TARGET_TO_MARKER_PRIORITY.items():
        if kw in all_target_text:
            for sn in set_names:
                if sn in MARKER_SETS and sn not in priority_sets:
                    priority_sets.append(sn)

    # Add remaining marker sets
    remaining = [(k, v) for k, v in MARKER_SETS.items() if k not in priority_sets]
    result = [(sn, MARKER_SETS[sn]) for sn in priority_sets] + remaining
    return result


def annotate_markers(
    adata: "ad.AnnData",
    target_cell_types: Optional[list[str]] = None,
    min_score: float = 0.1,
    refine: bool = False,
    refine_min_score: float = 0.3,
    max_refine_cells: int = 200000,
) -> dict:
    """Assign cell types using marker gene expression (sparse-safe).

    Two modes:
      1. Fill mode (refine=False): only annotate cells missing cell_type.
      2. Refine mode (refine=True): re-evaluate cells whose current label
         does not match any target AND whose current confidence is low.
         Original labels never overwritten. Refinement only applied when
         marker score >= refine_min_score (default 0.3).

    Score = fraction of markers detected > 0.
    """
    import re as _re

    if "cell_type" not in adata.obs.columns:
        adata.obs["cell_type"] = None
    if "annotation_method" not in adata.obs.columns:
        adata.obs["annotation_method"] = "Unresolved"

    ct_values = adata.obs["cell_type"].astype(str)

    # Fill mask: empty/missing
    fill_mask = ct_values.isna() | ct_values.str.strip().str.lower().isin(
        ["", "nan", "none", "null"]
    )

    # Refine mask: only cells that:
    #   - Are NOT Original
    #   - Have a label that doesn't match any target
    #   - Have low confidence (or no confidence column)
    #   - Are not Unknown
    refine_mask = pd.Series(False, index=adata.obs.index)
    if refine and target_cell_types and len(target_cell_types) > 0:
        is_original = adata.obs["annotation_method"] == "Original"
        target_norm = {_normalize(t) for t in target_cell_types}
        current_norm = ct_values.apply(_normalize)
        in_target = current_norm.apply(lambda x: any(
            t in x or x in t for t in target_norm
        ))
        # Only refine if CellTypist confidence is low (or unavailable)
        low_conf = pd.Series(True, index=adata.obs.index)
        if "cell_type_confidence" in adata.obs.columns:
            conf = adata.obs["cell_type_confidence"]
            low_conf = conf.isna() | (conf < 0.5)

        refine_mask = (
            ~is_original
            & ~fill_mask
            & ~in_target
            & low_conf
            & (ct_values != "Unknown")
        )

    # Cap refine candidates if too many
    n_refine_raw = refine_mask.sum()
    if n_refine_raw > max_refine_cells:
        logger.info(
            "Limiting refinement to %d / %d candidates (max_refine_cells=%d)",
            max_refine_cells, n_refine_raw, max_refine_cells,
        )
        refine_indices = refine_mask.index[refine_mask]
        sampled = np.random.choice(
            refine_indices, size=max_refine_cells, replace=False,
        )
        refine_mask_limited = pd.Series(False, index=adata.obs.index)
        refine_mask_limited.loc[sampled] = True
        refine_mask = refine_mask_limited

    mask = fill_mask | refine_mask
    n_missing = fill_mask.sum()
    n_refine = refine_mask.sum()

    if mask.sum() == 0:
        return {
            "method": "marker_heuristics",
            "cells_assessed": 0,
            "cells_labeled": 0,
            "cells_refined": 0,
        }

    # Build gene index
    gene_index = {_normalize_gene(g): i for i, g in enumerate(adata.var_names)}
    marker_sets = _get_prioritized_marker_sets(target_cell_types or [])

    # Subset to candidate cells
    X_subset = adata[mask].X
    n_cells_subset = X_subset.shape[0]

    best_type = np.full(n_cells_subset, "Unknown", dtype=object)
    best_score = np.zeros(n_cells_subset)

    for ct, markers in marker_sets:
        present_idx = [
            idx for m in markers
            if (idx := gene_index.get(m.upper())) is not None
        ]
        if not present_idx:
            continue

        marker_subset = X_subset[:, present_idx]
        if sp.issparse(marker_subset):
            scores = np.array((marker_subset > 0).sum(axis=1)).ravel() / len(present_idx)
        else:
            scores = np.mean(marker_subset > 0, axis=1)

        improve = scores > best_score
        best_type[improve] = ct
        best_score[improve] = scores[improve]

    # Apply threshold — stricter for refined cells
    mask_indices = mask.index[mask]
    def _get_mask_value(mask_series, idx):
        if idx not in mask_series.index:
            return False
        val = mask_series.loc[idx]
        if isinstance(val, pd.Series):
            return bool(val.iloc[0])
        return bool(val)

    is_refined_cell = np.array([
        _get_mask_value(refine_mask, idx)
        for idx in mask_indices
    ], dtype=bool)
    threshold_arr = np.where(is_refined_cell, refine_min_score, min_score)
    assignable = best_score >= threshold_arr
    # Refined cells that fail threshold: keep original label
    for i in range(n_cells_subset):
        if is_refined_cell[i] and not assignable[i]:
            orig_idx = mask_indices[i]
            best_type[i] = ct_values.loc[orig_idx]
    best_type[~assignable & ~is_refined_cell] = "Unknown"

    # Write back — track whether label actually changed
    old_labels = ct_values.loc[mask_indices].values
    labels = pd.Series(best_type, index=mask_indices)
    label_changed = labels.values != old_labels

    # cell_type
    adata.obs.loc[mask, "cell_type"] = labels.values

    # annotation_method: only write MarkerHeuristics if label actually changed
    # AND the cell wasn't Original. Cells keeping their previous label
    # retain their previous annotation_method.
    override_idx = mask_indices[labels != "Unknown"]
    if len(override_idx) > 0:
        changed = pd.Series(label_changed, index=mask_indices).loc[override_idx]
        is_orig = adata.obs.loc[override_idx, "annotation_method"] == "Original"
        # Write MarkerHeuristics only for non-Original cells whose label changed
        write_mh = override_idx[changed & ~is_orig]
        if len(write_mh) > 0:
            adata.obs.loc[write_mh, "annotation_method"] = "MarkerHeuristics"

    # Track: which refined cells actually changed vs kept original
    refined_changed = 0
    refined_kept = 0
    missing_filled = 0
    if refine:
        refined_indices = mask_indices[is_refined_cell]
        for i, idx in enumerate(refined_indices):
            if label_changed[list(mask_indices).index(idx)]:
                refined_changed += 1
            else:
                refined_kept += 1
        fill_indices = mask_indices[~np.array(is_refined_cell)]
        for idx in fill_indices:
            i = list(mask_indices).index(idx)
            if labels.iloc[i] != "Unknown":
                missing_filled += 1

    # confidence
    if "cell_type_confidence" not in adata.obs.columns:
        adata.obs["cell_type_confidence"] = np.nan
    conf_series = pd.Series(best_score, index=mask_indices)
    adata.obs.loc[mask, "cell_type_confidence"] = conf_series.values

    n_labeled = int((labels != "Unknown").sum())
    logger.info(
        "Marker heuristics: %d missing filled, %d/%d refinement candidates -> "
        "%d overridden, %d kept original",
        missing_filled, n_refine if refine else 0, n_refine if refine else 0,
        refined_changed, refined_kept,
    )

    return {
        "method": "marker_heuristics",
        "cells_assessed": int(mask.sum()),
        "cells_labeled": n_labeled,
        "missing_filled": missing_filled,
        "refinement_candidates": int(refine_mask.sum()) if refine else 0,
        "refinement_successful_overrides": refined_changed,
        "refinement_kept_original": refined_kept,
        "cell_types_assigned": sorted(set(str(t) for t in best_type[assignable])),
    }


# ── Cluster-level target refinement ─────────────────────────────
# Operates on CellTypist label groups (clusters), not individual cells.
# For each unmatched target, finds its broad-candidate CellTypist labels,
# computes target-specific marker scores at cluster level, and labels
# high-scoring cells within those clusters.

# Per-target marker sets for refinement (compact, specific)
TARGET_REFINEMENT_MARKERS = {
    # Keratinocyte subtypes
    "basal keratinocytes": ["KRT5", "KRT14", "COL17A1", "ITGA6"],
    "spinous keratinocytes": ["KRT1", "KRT10", "DSG1", "DSC1"],
    "supra-spinous keratinocytes": ["KRT1", "KRT10", "FLG", "LOR"],
    "proliferating keratinocytes": ["MKI67", "TOP2A", "PCNA", "KRT5"],
    # Immune
    "temra t cells": ["CD3D", "CD8A", "GZMB", "KLRG1", "NKG7"],
    "dc3-1": ["LAMP3", "CCR7", "FCER1A", "CLEC10A"],
    "dc3-2": ["LAMP3", "CCR7", "FCER1A", "CLEC10A"],
    # Vascular
    "plvap+ capillaries": ["PLVAP", "RGCC", "PECAM1", "VWF"],
    # Fibroblast subtypes
    "wnt5a+ fibroblasts": ["WNT5A", "IL24", "COL1A1", "DCN"],
    # Kidney tubule subtypes (keyword-matched in refine_by_target_markers)
    "_kidney_pt": ["SLC34A1", "LRP2", "CUBN", "SLC5A2", "SLC22A6", "MIOX",
                    "ALDOB", "FBP1", "PCK1", "GGT1", "AQP1", "PDZK1"],
    "_kidney_tal": ["SLC12A1", "UMOD", "KCNJ1", "CLDN16", "CLDN19", "BSND",
                     "CLCNKB", "CASR", "EGF"],
}


def refine_by_target_markers(
    adata: "ad.AnnData",
    target_cell_types: list[str],
    broad_candidate_map: dict,
    min_cluster_score: float = 0.15,
    min_cell_score: float = 0.2,
) -> dict:
    """Cluster-level target refinement.

    For each unmatched target, finds CellTypist label clusters that are
    broad candidates (from BROAD_CANDIDATE_MAP), computes target-specific
    marker scores, and labels qualifying cells.

    Does NOT overwrite Original labels. Only labels cells whose current
    label is a broad candidate.

    Returns summary dict.
    """
    if "cell_type" not in adata.obs.columns:
        return {"refined_targets": 0, "cells_labeled": 0}

    ct_values = adata.obs["cell_type"].astype(str)
    gene_index = {g.upper(): i for i, g in enumerate(adata.var_names)}
    refined_targets = []
    total_labeled = 0

    for target in target_cell_types:
        target_n = _normalize(target)
        markers = TARGET_REFINEMENT_MARKERS.get(target_n)
        candidate_labels = broad_candidate_map.get(target_n, [])

        # Keyword-based fallback for kidney tubule and other epithelial subtypes
        if not markers or not candidate_labels:
            kw_markers = None
            kw_candidates = []
            if any(kw in target_n for kw in ["pt ", "proximal tubul", "pt cell",
                                               "proximal convoluted"]):
                kw_markers = TARGET_REFINEMENT_MARKERS.get("_kidney_pt")
                kw_candidates = ["epithelial cells", "epithelial cell"]
            elif any(kw in target_n for kw in ["tal ", "thick ascend", "loop of henle",
                                                 "distal tubul", "macula densa"]):
                kw_markers = TARGET_REFINEMENT_MARKERS.get("_kidney_tal")
                kw_candidates = ["epithelial cells", "epithelial cell"]
            if kw_markers:
                markers = kw_markers
            if kw_candidates:
                candidate_labels = kw_candidates

        if not markers:
            continue
        # Find markers present in data
        present_idx = [gene_index[m] for m in markers if m in gene_index]
        if len(present_idx) < 2:
            continue

        # Find broad-candidate CellTypist labels for this target
        candidate_labels = broad_candidate_map.get(target_n, [])
        if not candidate_labels:
            continue

        candidate_labels_norm = {_normalize(cl) for cl in candidate_labels}
        # Find cells in these candidate clusters
        in_cluster = ct_values.apply(
            lambda x: _normalize(x) in candidate_labels_norm
        )
        if in_cluster.sum() == 0:
            continue

        # Compute marker scores for candidate cells only
        X_subset = adata[in_cluster].X
        if X_subset.shape[0] == 0:
            continue

        marker_expr = X_subset[:, present_idx]
        if sp.issparse(marker_expr):
            scores = np.array((marker_expr > 0).sum(axis=1)).ravel() / len(present_idx)
        else:
            scores = np.mean(marker_expr > 0, axis=1)

        # Label cells above threshold
        good = scores >= min_cell_score
        if good.sum() == 0:
            continue

        good_indices = in_cluster.index[in_cluster][good]
        is_original = adata.obs.loc[good_indices, "annotation_method"] == "Original"
        write_idx = good_indices[~is_original]

        if len(write_idx) == 0:
            continue

        # Write refined label
        adata.obs.loc[write_idx, "cell_type"] = target  # use original casing
        adata.obs.loc[write_idx, "annotation_method"] = "CellTypist+MarkerHeuristics"
        if "cell_type_confidence" in adata.obs.columns:
            adata.obs.loc[write_idx, "cell_type_confidence"] = scores[good][~is_original.values]

        total_labeled += len(write_idx)
        refined_targets.append(target)
        logger.info(
            "Refined '%s': %d / %d candidate cells labeled (score>=%.2f)",
            target, len(write_idx), in_cluster.sum(), min_cell_score,
        )

    return {
        "method": "cluster_refinement",
        "refined_targets": refined_targets,
        "cells_labeled": total_labeled,
    }


def main():
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
    summary = annotate_markers(adata)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    adata.write(output_path, compression="gzip")
    print(json.dumps(summary, indent=2))
    if summary_path:
        Path(summary_path).parent.mkdir(parents=True, exist_ok=True)
        with open(summary_path, "w") as f:
            json.dump(summary, f, indent=2)


if __name__ == "__main__":
    main()
