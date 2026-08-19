"""CellTypist-based cell-type annotation.

Uses CellTypist models to predict cell types from expression data.
Works on a COPY of the AnnData to avoid mutating the original .X.
"""

import json
import logging
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

CONTEXT_MODELS = {
    "immune_enriched": ["Immune_All_Low.pkl", "Immune_All_High.pkl"],
    "tumor_or_tme": ["Immune_All_Low.pkl"],
    "unknown": ["Immune_All_Low.pkl"],
    # normal_or_disease_tissue: selected dynamically from target keywords
    "normal_or_disease_tissue": None,  # None = auto-detect from targets
}

# Tissue-context model mapping for auto-selection
# Keyed by lowercase keyword; first match in target_cell_types wins.
TISSUE_MODEL_MAP = {
    "lung": ["Human_Lung_Atlas.pkl"],
    "pulmonary": ["Human_Lung_Atlas.pkl"],
    "alveolar": ["Human_Lung_Atlas.pkl"],
    "bronchial": ["Human_Lung_Atlas.pkl"],
    "skin": ["Adult_Human_Skin.pkl"],
    "keratinocyte": ["Adult_Human_Skin.pkl"],
    "fibroblast": ["Adult_Human_Skin.pkl"],
    "melanocyte": ["Adult_Human_Skin.pkl"],
    "epidermal": ["Adult_Human_Skin.pkl"],
    "dermal": ["Adult_Human_Skin.pkl"],
    "psoriasis": ["Adult_Human_Skin.pkl"],
    "endothelial": ["Adult_Human_Skin.pkl"],
    "pericyte": ["Adult_Human_Skin.pkl"],
    "brain": ["Developing_Human_Brain.pkl"],
    "neuronal": ["Developing_Human_Brain.pkl"],
    "neuron": ["Developing_Human_Brain.pkl"],
    "glial": ["Developing_Human_Brain.pkl"],
    "liver": ["Human_Liver_Atlas.pkl"],
    "hepatocyte": ["Human_Liver_Atlas.pkl"],
    "kidney": ["Human_Kidney_Atlas.pkl"],
    "pancreas": ["Human_Pancreas_Atlas.pkl"],
    "intestine": ["Human_Intestinal_Atlas.pkl"],
    "colon": ["Human_Colorectal_Atlas.pkl"],
    "heart": ["Human_Heart_Atlas.pkl"],
}
# Immune model is always tried as fallback for any tissue
IMMUNE_FALLBACK = "Immune_All_Low.pkl"

FALLBACK_MODEL = "Immune_All_Low.pkl"


def _select_model_for_tissue(target_cell_types: list[str]) -> str:
    """Select best CellTypist model based on target cell type keywords.

    Scans target keywords against TISSUE_MODEL_MAP. Returns the first
    matching tissue model, or FALLBACK_MODEL if nothing matches.
    """
    all_text = " ".join(t.lower() for t in target_cell_types)
    for kw, models in TISSUE_MODEL_MAP.items():
        if kw in all_text:
            model = models[0]
            logger.info("Auto-selected CellTypist model '%s' based on keyword '%s'", model, kw)
            return model
    logger.info("No tissue keyword matched — using fallback model '%s'", FALLBACK_MODEL)
    return FALLBACK_MODEL


def _check_celltypist() -> bool:
    try:
        import celltypist
        return True
    except ImportError:
        return False


def annotate_celltypist(
    adata: "ad.AnnData",
    strategy_category: str = "unknown",
    model_name: Optional[str] = None,
    target_cell_types: Optional[list[str]] = None,
    majority_voting: bool = True,
) -> dict:
    """Annotate cells using CellTypist.

    - Works on a temporary copy with normalized expression.
    - Auto-selects model from target_cell_types keywords (skin/lung/brain/etc.)
    - Use --celltypist-model to override.
    - Only annotates cells missing a cell_type label.
    - Writes cell_type, annotation_method, and cell_type_confidence.

    Returns summary dict for manifest.
    """
    if not _check_celltypist():
        return {
            "method": "celltypist",
            "status": "skipped",
            "reason": "celltypist_not_installed",
            "cells_labeled": 0,
        }

    import scanpy as sc
    import celltypist
    from celltypist import models

    if "cell_type" not in adata.obs.columns:
        adata.obs["cell_type"] = None
    if "annotation_method" not in adata.obs.columns:
        adata.obs["annotation_method"] = "Unresolved"

    mask = adata.obs["cell_type"].isna() | (
        adata.obs["cell_type"].astype(str).str.strip().str.lower().isin(
            ["", "nan", "none", "null"]
        )
    )
    n_missing = mask.sum()
    if n_missing == 0:
        return {"method": "celltypist", "cells_labeled": 0}

    if model_name is None:
        cat_models = CONTEXT_MODELS.get(strategy_category)
        if cat_models is None:
            # Auto-detect from target keywords
            model_name = _select_model_for_tissue(target_cell_types or [])
        else:
            model_name = cat_models[0]

    # Load model
    try:
        model = models.Model.load(model=model_name)
    except Exception:
        try:
            models.download_models(force_update=False, model=[model_name])
            model = models.Model.load(model=model_name)
        except Exception as e:
            return {
                "method": "celltypist",
                "status": "skipped",
                "reason": f"model_load_failed: {e}",
                "cells_labeled": 0,
            }

    # Build input copy: raw counts -> normalize_total + log1p for CellTypist
    import copy
    import anndata as ad
    if "counts" in adata.layers:
        input_X = adata.layers["counts"].copy()
    else:
        input_X = adata.X.copy()

    input_adata = ad.AnnData(
        X=input_X,
        obs=adata.obs[[]].copy() if adata.obs.shape[1] > 0 else pd.DataFrame(index=adata.obs_names),
        var=adata.var.copy(),
    )
    sc.pp.normalize_total(input_adata, target_sum=1e4)
    sc.pp.log1p(input_adata)

    # Run CellTypist
    try:
        predictions = celltypist.annotate(
            input_adata,
            model=model_name,
            majority_voting=majority_voting,
            mode="best match",
        )
    except Exception as e:
        return {
            "method": "celltypist",
            "status": "error",
            "reason": str(e),
            "cells_labeled": 0,
        }

    # Extract predictions from CellTypist return object.
    # Primary: use to_adata() (most reliable across CellTypist versions).
    # Fallback: access .predicted_labels attribute directly.
    pred_labels = None
    pred_conf = None
    try:
        if hasattr(predictions, "to_adata"):
            pred_adata = predictions.to_adata()
            for col in ["majority_voting", "predicted_labels"]:
                if col in pred_adata.obs.columns:
                    pred_labels = pred_adata.obs[col]
                    break
            if "conf_score" in pred_adata.obs.columns:
                pred_conf = pred_adata.obs["conf_score"]

        if pred_labels is None and hasattr(predictions, "predicted_labels"):
            pl = predictions.predicted_labels
            if isinstance(pl, pd.DataFrame):
                for col in ["majority_voting", "predicted_labels"]:
                    if col in pl.columns:
                        pred_labels = pl[col]
                        break
                if "conf_score" in pl.columns:
                    pred_conf = pl["conf_score"]
    except Exception as e:
        logger.warning("CellTypist prediction parsing fallback: %s", e)

    # Fallback: check if predictions wrote to input_adata.obs
    if pred_labels is None:
        for col in ["majority_voting", "predicted_labels"]:
            if col in input_adata.obs.columns:
                pred_labels = input_adata.obs[col]
                break
        if "conf_score" in input_adata.obs.columns:
            pred_conf = input_adata.obs["conf_score"]

    if pred_labels is None:
        return {
            "method": "celltypist",
            "status": "no_predictions",
            "cells_labeled": 0,
        }

    # Write predictions back to original adata (only for missing cells)
    write_mask = mask.values
    # pred_labels aligns to input_adata which has same order as adata
    cell_labels = pred_labels.values if isinstance(pred_labels, pd.Series) else pred_labels

    adata.obs.loc[mask, "cell_type"] = [
        str(cl) for cl, m in zip(cell_labels, write_mask) if m
    ]
    adata.obs.loc[mask, "annotation_method"] = [
        "CellTypist" if (not pd.isna(cl) and str(cl).strip() != "") else "Unresolved"
        for cl, m in zip(cell_labels, write_mask) if m
    ]

    if pred_conf is not None:
        if "cell_type_confidence" not in adata.obs.columns:
            adata.obs["cell_type_confidence"] = np.nan
        conf_vals = pred_conf.values if isinstance(pred_conf, pd.Series) else pred_conf
        if len(conf_vals) == len(cell_labels):
            adata.obs.loc[mask, "cell_type_confidence"] = [
                float(cv) for cv, m in zip(conf_vals, write_mask) if m
            ]

    n_labeled = int(mask.sum())
    logger.info(
        "CellTypist (%s): %s / %s cells labeled", model_name, n_labeled, n_missing,
    )

    unique_ct = len(set(
        str(c) for c in cell_labels if not pd.isna(c)
    ))

    return {
        "method": "celltypist",
        "model": model_name,
        "cells_labeled": n_labeled,
        "n_unique_cell_types": unique_ct,
    }


def main():
    import anndata as ad
    if len(sys.argv) < 3:
        print(f"Usage: python3 {__file__} input.h5ad output.h5ad ...")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]
    strategy_category = "unknown"
    model_name = None
    summary_path = None
    i = 3
    while i < len(sys.argv):
        if sys.argv[i] == "--strategy-category" and i + 1 < len(sys.argv):
            strategy_category = sys.argv[i + 1]; i += 2
        elif sys.argv[i] == "--model" and i + 1 < len(sys.argv):
            model_name = sys.argv[i + 1]; i += 2
        elif sys.argv[i] == "--summary" and i + 1 < len(sys.argv):
            summary_path = sys.argv[i + 1]; i += 2
        else:
            i += 1

    adata = ad.read_h5ad(input_path)
    summary = annotate_celltypist(
        adata, strategy_category=strategy_category, model_name=model_name,
    )
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    adata.write(output_path, compression="gzip")
    print(json.dumps(summary, indent=2))
    if summary_path:
        Path(summary_path).parent.mkdir(parents=True, exist_ok=True)
        with open(summary_path, "w") as f:
            json.dump(summary, f, indent=2)

if __name__ == "__main__":
    main()
