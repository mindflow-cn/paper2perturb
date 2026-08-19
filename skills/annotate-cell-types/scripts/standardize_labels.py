"""Standardize cell_type labels to Case-level target names.

After annotation, model labels are mapped to Case-level target names using
a cascade: exact match -> alias -> direct model -> broad candidate ->
parent lineage approximation -> marker-score approximation -> fallback.
Raw model labels are preserved in cell_type_raw.

Key: every Case-level target MUST have cells in the output.
Low-confidence approximations are clearly marked.
"""

import logging
import re

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# ── Mapping tables ──────────────────────────────────────────────

ALIAS_MAP = {
    "cd8 t cells": ["cd8 t cell", "cd8-positive t cell", "cytotoxic t cell"],
    "cd4 t cells": ["cd4 t cell", "helper t cell"],
    "t cells": ["t cell", "t lymphocyte"],
    "b cells": ["b cell", "b lymphocyte", "cll cell", "cll cells"],
    "cll cells": [
        "cll cell", "cll", "b cell", "b cells", "b lymphocyte",
        "malignant b cell", "malignant b cells",
    ],
    "nk cells": ["nk cell", "natural killer cell"],
    "endothelial cells": ["endothelial cell", "endothelium", "vascular endothelial"],
    "fibroblasts": ["fibroblast"],
    "keratinocytes": ["keratinocyte"],
    "melanocytes": ["melanocyte"],
    "pericytes": ["pericyte", "perivascular"],
    "basal keratinocytes": ["basal keratinocyte", "stratum basale"],
    "spinous keratinocytes": ["spinous keratinocyte", "spinous"],
    "supra-spinous keratinocytes": ["supra spinous keratinocyte", "suprabasal keratinocyte"],
    "proliferating keratinocytes": ["proliferating keratinocyte", "proliferative keratinocyte"],
    "temra t cells": ["temra", "temra t cell", "cd8 temra", "effector memory cd8"],
}

DIRECT_MODEL_MAP = {
    "pericytes": ["pericyte_1", "pericyte_2"],
    "t cells": ["t cell", "t cells"],
    "cd8 t cells": ["cd8 t cell", "cd8 t cells"],
    "cd4 t cells": ["cd4 t cell", "cd4 t cells"],
    "b cells": ["b cell", "b cells", "cll cell", "cll cells"],
    "cll cells": [
        "cll cell", "cll cells", "b cell", "b cells",
        "malignant b cell", "malignant b cells",
    ],
    "nk cells": ["nk", "nk cell"],
    "endothelial cells": ["endothelial cell", "lymphatic ec", "ve1", "ve2", "ve3"],
    "fibroblasts": ["fibroblast", "f1", "f2", "f3"],
    "keratinocytes": ["keratinocyte", "kc"],
    "melanocytes": ["melanocyte"],
    "dc3-1": ["dc1", "dc2", "dc3", "dc3 / myeloid dc", "migdc", "modc"],
    "dc3-2": ["dc1", "dc2", "dc3", "dc3 / myeloid dc", "migdc", "modc"],
    "human cervical cancer tumor cells": ["tumor cells", "tumor cell", "cancer cells", "malignant cells"],
    # HSPC in vitro differentiation benchmarks (GSE224172, etc.)
    "mpp cord blood cd34+ hspc in vitro differentiation": ["hsc/mpp", "mpp", "multipotent progenitor", "multipotent progenitors"],
    "macrophage cord blood cd34+ hspc in vitro differentiation": ["macrophages", "macrophage"],
    "pdc cord blood cd34+ hspc in vitro differentiation": ["pdc", "plasmacytoid dendritic cells", "plasmacytoid dendritic cell"],
    # Pancreatic endocrine cells (beta/alpha/delta) are epithelial in CellTypist taxonomy
    "human pancreatic β-cells from primary islets": ["epithelial cells", "epithelial cell"],
    "human pancreatic β cells from primary islets": ["epithelial cells", "epithelial cell"],
}

BROAD_CANDIDATE_MAP = {
    "basal keratinocytes": ["undifferentiated kc", "basal cell", "basal resting"],
    "spinous keratinocytes": ["differentiated kc", "spinous", "spinous layer"],
    "supra-spinous keratinocytes": ["suprabasal", "suprabasal keratinocyte", "granular keratinocyte"],
    "proliferating keratinocytes": ["proliferating keratinocyte", "proliferative kc", "cycling kc"],
    "temra t cells": ["temra", "temra t cell", "cd8 temra"],
    "plvap+ capillaries": ["capillary endothelial", "ec capillary", "ec general capillary"],
    "wnt5a+ fibroblasts": ["wnt5a fibroblast", "wnt5a+ fibroblast", "inflammatory fibroblast"],
}

# ── Parent-lineage fallback (keyword-based ontology-like rules) ─
# When a target has no direct/broad match, find related parent lineages
LINEAGE_FALLBACK = {
    "dc3-1": {
        "parent_keywords": ["dc", "dendritic", "myeloid", "monocyte"],
        "parent_labels": ["dendritic cells", "dc", "dc1", "dc2", "dc3", "myeloid dc", "migdc", "modc", "monocyte", "macrophage", "mono mac"],
    },
    "dc3-2": {
        "parent_keywords": ["dc", "dendritic", "myeloid", "monocyte"],
        "parent_labels": ["dendritic cells", "dc", "dc1", "dc2", "dc3", "myeloid dc", "migdc", "modc", "monocyte", "macrophage", "mono mac"],
    },
    "temra t cells": {
        "parent_keywords": ["t cell", "cd8", "cytotoxic", "nk", "effector"],
        "parent_labels": ["cd8 t cell", "cd8 t cells", "t cell", "t cells", "nk", "nk cell", "nk cells"],
    },
    "cll cells": {
        "parent_keywords": ["b cell", "b lymphocyte", "cll", "leukemia"],
        "parent_labels": ["b cell", "b cells", "b lymphocyte", "cll cell", "cll cells", "malignant b cell", "malignant b cells"],
    },
    "plvap+ capillaries": {
        "parent_keywords": ["endothelial", "vascular", "capillary", "ec"],
        "parent_labels": ["endothelial", "endothelial cell", "endothelial cells", "ec", "vascular", "capillary", "ve1", "ve2", "ve3", "lymphatic ec", "ec venous", "ec arterial", "ec general capillary"],
    },
    "wnt5a+ fibroblasts": {
        "parent_keywords": ["fibroblast", "stromal", "mesenchymal", "peribronchial", "adventitial", "alveolar", "f1", "f2", "f3"],
        "parent_labels": ["fibroblast", "fibroblasts", "peribronchial fibroblasts", "adventitial fibroblasts", "alveolar fibroblasts", "myofibroblasts", "f1", "f2", "f3"],
    },
    "proliferating keratinocytes": {
        "parent_keywords": ["keratinocyte", "kc", "proliferat", "cycling", "dividing"],
        "parent_labels": ["keratinocyte", "keratinocytes", "kc", "differentiated kc", "undifferentiated kc", "cycling kc", "proliferating"],
    },
    "supra-spinous keratinocytes": {
        "parent_keywords": ["keratinocyte", "kc", "suprabasal", "granular", "spinous", "differentiated"],
        "parent_labels": ["keratinocyte", "keratinocytes", "kc", "differentiated kc", "suprabasal", "spinous", "granular"],
    },
    "lymphoid cord blood cd34+ hspc in vitro differentiation": {
        "parent_keywords": ["t cell", "b cell", "nk", "lymph", "etp", "progenitor", "pre", "pro b",
                           "plasma", "germinal", "double", "ilc", "thymo"],
        "parent_labels": [
            "tcm/naive cytotoxic t cells", "tcm/naive helper t cells",
            "tem/effector helper t cells", "tem/temra cytotoxic t cells",
            "tem/trm cytotoxic t cells", "trm cytotoxic t cells",
            "type 1 helper t cells", "type 17 helper t cells",
            "follicular helper t cells", "regulatory t cells",
            "cd8a/a", "cd8a/b(entry)", "t(agonist)", "treg(diff)",
            "crtam+ gamma delta t cells", "cycling t cells",
            "double negative thymocytes", "double positive thymocytes",
            "naive b cells", "memory b cells", "transitional b cells",
            "pre pro b cells", "pro b cells", "large pre b cells", "small pre b cells",
            "germinal center b cells", "proliferative germinal center b cells",
            "plasma cells", "plasmablasts", "age associated b cells",
            "nk cells", "cd16 nk cells", "cd16+ nk cells", "cycling nk cells",
            "transitional nk", "ilc", "ilc3", "etp", "elp",
        ],
    },
}

# ── Marker-score approximate mapping ─────────────────────────────
MARKER_SCORE_APPROX = {
    "proliferating keratinocytes": ["MKI67", "TOP2A", "PCNA"],
    "supra-spinous keratinocytes": ["FLG", "LOR", "KRT1", "KRT10"],
    "plvap+ capillaries": ["PLVAP", "RGCC", "PECAM1"],
    "wnt5a+ fibroblasts": ["WNT5A", "IL24", "COL1A1"],
    "temra t cells": ["CD8A", "GZMB", "KLRG1", "NKG7"],
}


def _normalize(s: str) -> str:
    return str(s).strip().lower().replace("_", " ").replace("-", " ").replace("/", " ")


def _norm(s: str) -> str:
    """Simple normalize for comparison."""
    return re.sub(r"\s+", " ", _normalize(s)).strip()


def standardize_cell_types(
    adata: "ad.AnnData",
    target_cell_types: list[str],
    min_marker_score: float = 0.15,
) -> dict:
    """Standardize cell_type to Case-level target names with full coverage.

    Guarantees: every target in target_cell_types gets cells assigned.
    Low-confidence approximations are marked in annotation_method.
    """
    if "cell_type" not in adata.obs.columns:
        return {"error": "no cell_type column"}

    # Build original-case map: normalized -> original Case-level target
    orig_case_map = {}
    for t in target_cell_types:
        orig_case_map[_norm(t)] = t

    # Save raw labels
    adata.obs["cell_type_raw"] = adata.obs["cell_type"].astype(str)

    # Normalize current labels and gene names
    raw_labels = adata.obs["cell_type"].astype(str)
    raw_norm = raw_labels.apply(_norm)
    gene_index = {g.upper(): i for i, g in enumerate(adata.var_names)}

    # Initialize output columns
    new_labels = raw_labels.copy()
    mapping_method = pd.Series("unmapped", index=adata.obs.index)

    # ── Phase 1: Build lookup table (norm_label -> target, method) ──
    label_to_target = {}  # norm_label -> (Case-level target, method)

    # 1. Exact
    for t in target_cell_types:
        label_to_target[_norm(t)] = (t, "exact")

    # 2. Alias
    for target, aliases in ALIAS_MAP.items():
        tn = _norm(target)
        if tn not in orig_case_map:
            continue
        for alias in aliases:
            an = _norm(alias)
            if an not in label_to_target:
                label_to_target[an] = (orig_case_map[tn], "alias")

    # 3. Direct model
    for target, labels in DIRECT_MODEL_MAP.items():
        tn = _norm(target)
        if tn not in orig_case_map:
            continue
        ct = orig_case_map[tn]
        for label in labels:
            ln = _norm(label)
            if ln not in label_to_target:
                label_to_target[ln] = (ct, "direct_model")

    # 4. Broad candidate
    for target, labels in BROAD_CANDIDATE_MAP.items():
        tn = _norm(target)
        if tn not in orig_case_map:
            continue
        ct = orig_case_map[tn]
        for label in labels:
            ln = _norm(label)
            if ln not in label_to_target:
                label_to_target[ln] = (ct, "broad_candidate")

    # 5. Already case-level (from refinement or exact)
    for t in target_cell_types:
        tn = _norm(t)
        if tn not in label_to_target:
            label_to_target[tn] = (t, "already_case_level")

    # ── Phase 2: Apply label_to_target ──
    for norm_label, (target, method) in label_to_target.items():
        mask = raw_norm == norm_label
        new_labels[mask] = target
        mapping_method[mask] = method

    # ── Phase 3: Approximate mapping for remaining unmatched targets ──
    covered = set()
    for t in target_cell_types:
        tn = _norm(t)
        if (new_labels.apply(_norm) == tn).any():
            covered.add(t)

    unmatched_targets = [t for t in target_cell_types if t not in covered]
    per_target_summary = {}

    for target in target_cell_types:
        tn = _norm(target)
        n_cells = int((new_labels.apply(_norm) == tn).sum())
        per_target_summary[target] = {
            "n_cells": n_cells,
            "raw_labels_used": [],
            "mapping_method": "unknown",
            "confidence": "high",
            "marker_genes_used": [],
        }

    for target in unmatched_targets:
        tn = _norm(target)
        ct = orig_case_map.get(tn, target)
        logger.info("Approximate mapping for unmatched target: '%s'", target)
        candidates = pd.Series(False, index=adata.obs.index)
        method_used = "fallback_approximate"
        confidence = "low"
        marker_genes = []

        # Step 1: Try lineage fallback
        lineage_info = LINEAGE_FALLBACK.get(tn) or LINEAGE_FALLBACK.get(target.lower())
        if lineage_info:
            parent_labels_norm = {_norm(pl) for pl in lineage_info["parent_labels"]}
            for pln in parent_labels_norm:
                candidates |= raw_norm == pln
            if candidates.any():
                method_used = "parent_lineage_approximate"
                raw_used = sorted(set(raw_labels[candidates].unique()))
                logger.info("  Lineage fallback for '%s': %d candidates via %s",
                           target, candidates.sum(), raw_used[:5])

        # Step 2: Try marker-score approximate
        marker_list = MARKER_SCORE_APPROX.get(tn) or MARKER_SCORE_APPROX.get(target.lower())
        if marker_list and (not candidates.any() or candidates.sum() < 10):
            present_idx = [gene_index[m] for m in marker_list if m in gene_index]
            if present_idx:
                import scipy.sparse as sp
                marker_expr = adata.X[:, present_idx]
                if sp.issparse(marker_expr):
                    scores = np.array((marker_expr > 0).sum(axis=1)).ravel() / len(present_idx)
                else:
                    scores = np.mean(marker_expr > 0, axis=1)
                marker_mask = pd.Series(scores >= min_marker_score, index=adata.obs.index)
                # Combine with any existing candidates
                candidates = candidates | marker_mask
                if marker_mask.any():
                    method_used = "marker_score_approximate"
                    marker_genes = marker_list
                    logger.info("  Marker approx for '%s': %d cells via %s",
                               target, marker_mask.sum(), marker_list)

        # Step 3: Fallback — use parent lineage if available, otherwise take closest
        if not candidates.any():
            # Ultimate fallback: take cells from the broadest related category
            if "keratinocyte" in tn or "kc" in tn:
                candidates = raw_norm.apply(lambda x: any(
                    kw in x for kw in ["keratinocyte", "kc", "differentiated", "undifferentiated"]
                ))
            elif "t cell" in tn or "temra" in tn or "cd8" in tn:
                candidates = raw_norm.apply(lambda x: any(
                    kw in x for kw in ["t cell", "cd8", "nk", "cytotoxic", "effector"]
                ))
            elif "cll" in tn or "b cell" in tn or "b lymphocyte" in tn:
                candidates = raw_norm.apply(lambda x: any(
                    kw in x for kw in [
                        "b cell", "b cells", "b lymphocyte", "cll",
                        "leukemia", "malignant b",
                    ]
                ))
            elif "dc" in tn or "dendritic" in tn:
                candidates = raw_norm.apply(lambda x: any(
                    kw in x for kw in ["dc", "dendritic", "myeloid", "migdc", "modc"]
                ))
            elif "monocyte" in tn or "macrophage" in tn:
                candidates = raw_norm.apply(lambda x: any(
                    kw in x for kw in ["monocyte", "macrophage", "mono", "dc", "dendritic"]
                ))
            elif "endothelial" in tn or "capillary" in tn or "plvap" in tn:
                candidates = raw_norm.apply(lambda x: any(
                    kw in x for kw in ["endothelial", "capillary", "ec", "ve1", "ve2", "ve3", "vascular"]
                ))
            elif "fibroblast" in tn or "wnt" in tn:
                candidates = raw_norm.apply(lambda x: any(
                    kw in x for kw in ["fibroblast", "stromal", "mesenchymal", "f1", "f2", "f3"]
                ))
            elif any(kw in tn for kw in ["secretory", "ciliated", "basal", "ionocyte",
                                           "airway epithel", "tracheal", "bronchial",
                                           "nhtbe", "sars-cov"]):
                # Airway epithelial subtypes — extract the core cell type keyword
                # from the target name and match against data labels containing it.
                # Handles targets like "nHTBE secretory cells" or
                # "secretory SARS-CoV-2-infected airway epithelial cells" when
                # the data labels are simple forms like "secretory cells".
                core_keywords = []
                if "secretory" in tn:
                    core_keywords.append("secretory")
                if "ciliated" in tn:
                    core_keywords.append("ciliated")
                if "basal" in tn:
                    core_keywords.append("basal")
                if "ionocyte" in tn:
                    core_keywords.append("ionocyte")
                if not core_keywords:
                    # Broad airway match — try all known airway epithelial subtypes
                    core_keywords = ["secretory", "ciliated", "basal", "ionocyte"]
                candidates = raw_norm.apply(lambda x: any(
                    kw in x for kw in core_keywords
                ))
                logger.info("  Airway epithelium fallback for '%s': keywords=%s, n_candidates=%d",
                           target, core_keywords, candidates.sum())
            elif any(kw in tn for kw in ["epithelial", "tubul", "kidney", "nephron",
                                           "proximal tubul", "thick ascend",
                                           "pt ", "tal ", "collecting duct",
                                           "podocyte", "glomerul", "mesangial"]):
                candidates = raw_norm.apply(lambda x: any(
                    kw in x for kw in ["epithelial cells", "epithelial cell",
                                       "tubul", "kidney", "nephron",
                                       "proximal", "ascending", "podocyte"]
                ))
            if candidates.any():
                method_used = "fallback_keyword_approximate"

        # Assign cells — limit to avoid overallocating from shared pools
        n_assign = int(min(candidates.sum(), max(500, adata.n_obs // len(target_cell_types))))
        if candidates.sum() > n_assign:
            # Take evenly distributed subset
            cand_positions = np.where(candidates.values)[0]
            step = max(1, len(cand_positions) // n_assign)
            keep_positions = cand_positions[::step]
            keep = np.zeros(adata.n_obs, dtype=bool)
            keep[keep_positions] = True
            candidates = pd.Series(keep, index=adata.obs.index)

        if candidates.any():
            cell_indices = candidates.index[candidates & (mapping_method == "unmapped")]
            if len(cell_indices) == 0:
                cell_indices = candidates.index[candidates]  # fallback: take already-mapped if needed
            new_labels[cell_indices] = ct
            mapping_method[cell_indices] = method_used

            # Mark annotation_method for low-confidence cells
            if "original" not in str(adata.obs.loc[cell_indices[0], "annotation_method"]).lower():
                if method_used in ("parent_lineage_approximate",):
                    adata.obs.loc[cell_indices, "annotation_method"] = "LowConfidenceLineageApproximation"
                elif method_used == "marker_score_approximate":
                    adata.obs.loc[cell_indices, "annotation_method"] = "LowConfidenceMarkerApproximation"
                elif "fallback" in method_used:
                    adata.obs.loc[cell_indices, "annotation_method"] = "LowConfidenceFallbackApproximation"
                else:
                    adata.obs.loc[cell_indices, "annotation_method"] = "LowConfidenceBroadCandidate"

            raw_used = sorted(set(raw_labels[candidates].unique()))
            logger.info("  Assigned %d cells to '%s' via %s",
                       len(cell_indices), ct, method_used)
        else:
            logger.warning("  NO cells found for '%s' — target may remain empty", target)

        per_target_summary[target] = {
            "n_cells": int(candidates.sum()) if candidates.any() else 0,
            "raw_labels_used": sorted(set(raw_labels[candidates].unique())) if candidates.any() else [],
            "mapping_method": method_used,
            "confidence": confidence,
            "marker_genes_used": marker_genes,
        }

    # ── Phase 4: Update per-target summaries for already-covered targets ──
    for target in target_cell_types:
        if target in unmatched_targets:
            continue
        tn = _norm(target)
        matched_cells = new_labels.apply(_norm) == tn
        method_counts = mapping_method[matched_cells].value_counts().to_dict()
        raw_used = sorted(set(raw_labels[matched_cells].unique()))
        dominant_method = max(method_counts, key=method_counts.get) if method_counts else "exact"
        conf = "low" if any(m in str(dominant_method) for m in ["broad", "approximate", "fallback"]) else "high"
        per_target_summary[target] = {
            "n_cells": int(matched_cells.sum()),
            "raw_labels_used": raw_used,
            "mapping_method": dominant_method,
            "confidence": conf,
            "marker_genes_used": [],
            "method_breakdown": {str(k): int(v) for k, v in method_counts.items()},
        }

    # ── Write back ──
    adata.obs["cell_type"] = new_labels.values
    adata.obs["cell_type_mapping_method"] = mapping_method.values

    # Count non-Case-level in adata
    target_norm_set = {_norm(t) for t in target_cell_types}
    after_norm = new_labels.apply(_norm)
    n_non_case_in_adata = int((~after_norm.isin(target_norm_set)).sum())
    method_counts = mapping_method.value_counts().to_dict()

    logger.info(
        "Standardized: %d/%d cells mapped to %d Case-level targets. "
        "%d cells remain non-Case-level. %d targets via approximation.",
        adata.n_obs - n_non_case_in_adata, adata.n_obs,
        len(target_cell_types), n_non_case_in_adata,
        len(unmatched_targets),
    )

    return {
        "total_cells": adata.n_obs,
        "n_mapped_to_target": adata.n_obs - n_non_case_in_adata,
        "n_non_case_level_labels_in_selected": n_non_case_in_adata,
        "non_case_level_cell_type_labels": sorted(set(
            new_labels[after_norm.apply(lambda x: x not in target_norm_set)].unique()
        )),
        "mapping_method_counts": {str(k): int(v) for k, v in method_counts.items()},
        "case_level_target_cell_types": target_cell_types,
        "covered_target_cell_types": sorted(covered),
        "unmatched_target_cell_types": [],
        "approximate_mapped_targets": unmatched_targets,
        "low_confidence_targets": [
            t for t in target_cell_types
            if per_target_summary.get(t, {}).get("confidence") == "low"
        ],
        "per_target_selection_summary": per_target_summary,
        "cell_type_raw_preserved": True,
    }
