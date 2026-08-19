"""Select cells matching Case-level target cell types.

Matching accumulates across match types per target:
  1. Exact normalized match
  2. Alias match (from ALIAS_MAP)
  3. Direct model-label match (DIRECT_MODEL_MAP — reliable 1:1)
  4. Broad candidate match (BROAD_CANDIDATE_MAP — ambiguous, needs refinement)
  5. Parent-child match
  6. Cell-line containment match

Narrow targets (CD8 T cells) must NOT match broad parents (all T cells).
Broad targets (T cells) MAY match narrow children (CD8 T cells, CD4 T cells).

Quality tracking: records whether each target was matched directly, via
broad candidate, or not at all.
"""

import logging
import re

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# ── Direct model-label map: reliable 1:1 mappings ──────────────
# Only map when the model label is a clear synonym for the target.
# DO NOT add broad labels like "Differentiated_KC -> spinous KC" here.
DIRECT_MODEL_MAP = {
    "pericytes": [
        "pericyte", "pericyte_1", "pericyte_2",
    ],
    "t cells": [
        "t cell", "t cells",
    ],
    "cd8 t cells": [
        "cd8 t cell", "cd8 t cells", "cd8+ t cell",
    ],
    "cd4 t cells": [
        "cd4 t cell", "cd4 t cells", "cd4+ t cell",
    ],
    "b cells": [
        "b cell", "b cells",
    ],
    "nk cells": [
        "nk", "nk cell", "nk cells",
    ],
    "plasma cells": [
        "plasma cell", "plasma cells", "plasmablast",
    ],
    "neutrophils": [
        "neutrophil", "neutrophils",
    ],
    "monocytes / macrophages": [
        "monocyte", "macrophage",
    ],
    "endothelial cells": [
        "endothelial cell", "endothelial cells",
    ],
    "keratinocytes": [
        "keratinocyte", "keratinocytes",
    ],
    "melanocytes": [
        "melanocyte", "melanocytes",
    ],
    "human cervical cancer tumor cells": [
        "tumor cells", "tumor cell",
        "cancer cells", "malignant cells",
    ],
    "fibroblasts": [
        "fibroblast", "fibroblasts",
    ],
}

# ── Broad candidate map: model label -> target cluster candidate ─
# These model labels are BROAD categories that COULD include the target.
# They should only be used to narrow the candidate pool for marker
# refinement, NOT treated as confirmed matches.
BROAD_CANDIDATE_MAP = {
    "basal keratinocytes": [
        "undifferentiated kc", "basal cell", "basal resting",
    ],
    "spinous keratinocytes": [
        "differentiated kc", "spinous", "spinous layer",
    ],
    "supra-spinous keratinocytes": [
        "suprabasal", "suprabasal keratinocyte", "granular keratinocyte",
    ],
    "proliferating keratinocytes": [
        "proliferating keratinocyte", "proliferative kc", "cycling kc",
    ],
    "temra t cells": [
        "temra", "temra t cell", "cd8 temra",
    ],
    "dc3-1": [
        "dc1", "dc2", "dc3", "dc3 dendritic", "dc3 / myeloid dc",
    ],
    "dc3-2": [
        "dc1", "dc2", "dc3", "dc3 dendritic", "dc3 / myeloid dc",
    ],
    "dendritic cells / dc3": [
        "dc1", "dc2", "dc3", "dc3 dendritic", "dc3 / myeloid dc",
    ],
    "plvap+ capillaries": [
        "capillary endothelial", "ec capillary", "ec general capillary",
    ],
    "wnt5a+ fibroblasts": [
        "wnt5a fibroblast", "wnt5a+ fibroblast", "inflammatory fibroblast",
    ],
    "proliferating cells": [
        "proliferating", "cycling", "dividing",
    ],
}

# ── Alias map (unchanged) ──────────────────────────────────────

ALIAS_MAP = {
    "cd8 t cells": [
        "cd8 t cell", "cd8-positive t cell", "cytotoxic t cell",
        "cd8+ t cell", "cd8+ t cells",
    ],
    "cd4 t cells": [
        "cd4 t cell", "helper t cell", "cd4+ t cell", "cd4+ t cells",
    ],
    "t cells": [
        "t cell", "t lymphocyte", "t-cells", "t-lymphocytes",
    ],
    "b cells": ["b cell", "b lymphocyte", "b-cells", "b-lymphocytes", "cll cells", "cll cell"],
    "nk cells": ["nk cell", "natural killer cell", "natural killer"],
    "macrophages": ["macrophage"],
    "monocytes": ["monocyte"],
    "dendritic cells": ["dendritic cell", "dc", "dcs"],
    "endothelial cells": [
        "endothelial cell", "endothelium", "capillaries", "capillary",
        "endothelial",
    ],
    "fibroblasts": ["fibroblast"],
    "epithelial cells": ["epithelial cell", "epithelium"],
    "neutrophils": ["neutrophil"],
    "plasma cells": ["plasma cell", "plasmablast", "plasmablasts"],
    "pericytes": ["pericyte", "perivascular"],
    "keratinocytes": ["keratinocyte"],
    "melanocytes": ["melanocyte"],
}

PARENT_CHILD_MAP = {
    "t cells": [
        "cd4 t cells", "cd8 t cells", "treg", "th17", "th1", "th2",
        "tcm", "tem", "temra", "gamma delta t cells", "mait",
        "cytotoxic t cells", "helper t cells", "temra t cells",
    ],
    "b cells": ["plasma cells", "memory b cells", "naive b cells"],
    "keratinocytes": [
        "basal keratinocytes", "spinous keratinocytes",
        "supra-spinous keratinocytes", "proliferating keratinocytes",
        "differentiated keratinocytes",
    ],
    "dendritic cells": [
        "dendritic cells / dc3", "dc1", "dc2", "dc3",
        "plasmacytoid dendritic cells",
    ],
    "endothelial cells": [
        "plvap+ capillaries", "capillary endothelial",
        "lymphatic endothelial", "vascular endothelial",
    ],
    "fibroblasts": [
        "wnt5a+ fibroblasts", "inflammatory fibroblasts",
        "myofibroblasts", "dermal fibroblasts",
    ],
}


def _normalize_label(label: str) -> str:
    s = str(label).strip().lower()
    s = s.replace("_", " ").replace("-", " ")
    s = re.sub(r"\s+", " ", s)
    return s


def _exact_match(cell_labels: pd.Series, target: str) -> pd.Series:
    target_n = _normalize_label(target)
    return cell_labels.apply(_normalize_label) == target_n


def _alias_match(cell_labels: pd.Series, target: str) -> pd.Series:
    target_n = _normalize_label(target)
    aliases = set(ALIAS_MAP.get(target_n, []))
    # Reverse lookup: if target is an alias of a canonical type, add that
    # type's full alias set (including the canonical name itself).
    for canonical, alias_list in ALIAS_MAP.items():
        if target_n in [_normalize_label(a) for a in alias_list]:
            aliases.add(canonical)
            aliases.update(_normalize_label(a) for a in alias_list)
    aliases.add(target_n)
    if target_n.endswith("s"):
        aliases.add(target_n[:-1])
    else:
        aliases.add(target_n + "s")
    normalized_labels = cell_labels.apply(_normalize_label)
    mask = pd.Series(False, index=cell_labels.index)
    for alias in aliases:
        mask |= normalized_labels == alias
    return mask


def _direct_model_match(cell_labels: pd.Series, target: str) -> pd.Series:
    """Direct 1:1 model-label match (reliable)."""
    target_n = _normalize_label(target)
    model_labels = DIRECT_MODEL_MAP.get(target_n, [])
    if not model_labels:
        return pd.Series(False, index=cell_labels.index)
    normalized_labels = cell_labels.apply(_normalize_label)
    mask = pd.Series(False, index=cell_labels.index)
    for ml in model_labels:
        mask |= normalized_labels == _normalize_label(ml)
    return mask


def _broad_candidate_match(cell_labels: pd.Series, target: str) -> pd.Series:
    """Broad candidate match (ambiguous, needs refinement)."""
    target_n = _normalize_label(target)
    model_labels = BROAD_CANDIDATE_MAP.get(target_n, [])
    if not model_labels:
        return pd.Series(False, index=cell_labels.index)
    normalized_labels = cell_labels.apply(_normalize_label)
    mask = pd.Series(False, index=cell_labels.index)
    for ml in model_labels:
        mask |= normalized_labels == _normalize_label(ml)
    return mask


def _parent_child_match(cell_labels: pd.Series, target: str) -> pd.Series:
    target_n = _normalize_label(target)
    children = PARENT_CHILD_MAP.get(target_n, [])
    if not children:
        return pd.Series(False, index=cell_labels.index)
    normalized_labels = cell_labels.apply(_normalize_label)
    mask = pd.Series(False, index=cell_labels.index)
    for child in children:
        mask |= normalized_labels == _normalize_label(child)
        for alias in ALIAS_MAP.get(_normalize_label(child), []):
            mask |= normalized_labels == alias
    return mask


def _containment_match(
    cell_labels: pd.Series, target: str, strategy_category: str,
) -> pd.Series:
    if strategy_category != "cell_line_or_in_vitro":
        return pd.Series(False, index=cell_labels.index)
    target_n = _normalize_label(target)
    normalized_labels = cell_labels.apply(_normalize_label)
    return normalized_labels.str.contains(target_n, regex=False, na=False)


def select_cells(
    adata: "ad.AnnData",
    target_cell_types: list[str],
    strategy_category: str = "unknown",
    allow_empty: bool = False,
) -> tuple["ad.AnnData", dict]:
    """Select cells matching target cell types.

    For each target, accumulates: exact, alias, direct-model, broad-candidate,
    parent-child, and containment matches.

    Tracks match quality: which targets were directly matched vs only via
    broad candidates vs unmatched.

    Returns (selected_adata, selection_summary).
    """
    if "cell_type" not in adata.obs.columns:
        raise ValueError("adata.obs lacks 'cell_type' column — annotate first")

    n_before = adata.n_obs
    if not target_cell_types:
        if allow_empty:
            empty = adata[0:0].copy()
            return empty, {
                "n_cells_before_selection": n_before,
                "n_cells_after_selection": 0,
                "matched_cell_type_labels": [],
                "unmatched_target_cell_types": [],
                "available_cell_types": [],
            }
        raise ValueError("No target cell types specified and --allow-empty-selection not set")

    cell_labels = adata.obs["cell_type"].astype(str)
    available = sorted(set(_normalize_label(l) for l in cell_labels.unique()))

    match_mask = pd.Series(False, index=adata.obs.index)
    directly_matched_targets = []
    broad_matched_targets = []
    unmatched_targets = []

    for target in target_cell_types:
        target_n = _normalize_label(target)
        target_match = pd.Series(False, index=adata.obs.index)
        direct_hit = False
        broad_hit = False

        # 1. Exact
        m = _exact_match(cell_labels, target)
        if m.any():
            target_match |= m
            direct_hit = True

        # 2. Alias
        m = _alias_match(cell_labels, target)
        if m.any():
            target_match |= m
            direct_hit = True

        # 3. Direct model-label (reliable)
        m = _direct_model_match(cell_labels, target)
        if m.any():
            target_match |= m
            direct_hit = True

        # 4. Broad candidate (ambiguous — flag for refinement)
        m = _broad_candidate_match(cell_labels, target)
        if m.any():
            target_match |= m
            broad_hit = True

        # 5. Parent-child
        m = _parent_child_match(cell_labels, target)
        if m.any():
            target_match |= m
            direct_hit = True

        # 6. Containment
        m = _containment_match(cell_labels, target, strategy_category)
        if m.any():
            target_match |= m
            direct_hit = True

        if target_match.any():
            match_mask |= target_match
            if direct_hit:
                directly_matched_targets.append(target_n)
            elif broad_hit:
                broad_matched_targets.append(target_n)
        else:
            unmatched_targets.append(target_n)

    n_after = int(match_mask.sum())

    present_labels = sorted(set(
        _normalize_label(l) for l in cell_labels[match_mask].unique()
    ))

    result = {
        "n_cells_before_selection": n_before,
        "n_cells_after_selection": n_after,
        "matched_cell_type_labels": present_labels,
        "unmatched_target_cell_types": unmatched_targets,
        "available_cell_types": available,
        "directly_matched_targets": directly_matched_targets,
        "broad_candidate_targets": broad_matched_targets,
    }

    if n_after == 0 and not allow_empty:
        raise ValueError(
            f"No cells matched target cell types {target_cell_types}. "
            f"Available cell types in data: {available}. "
            f"Use --allow-empty-selection to permit empty output."
        )

    selected = adata[match_mask].copy()
    logger.info(
        "Cell selection: %s / %s cells matched (direct=%d, broad=%d, unmatched=%d)",
        n_after, n_before, len(directly_matched_targets),
        len(broad_matched_targets), len(unmatched_targets),
    )

    return selected, result
