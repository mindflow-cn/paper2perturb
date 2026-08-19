"""Generate annotation manifest JSON for a benchmark.

Complete manifest with input detection, strategy, methods used/skipped/failed,
selection counts, warnings, and dependency status.
"""

import json
import logging
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


def build_manifest(
    dataset_id: str,
    benchmark_id: str,
    input_type: str,
    annotated_h5ad: str,
    selected_h5ad: str,
    target_cell_types: list[str],
    strategy_category: str,
    annotation_methods_used: list[str],
    selection_summary: dict,
    has_original_annotation: bool = False,
    needs_manual_review: bool = False,
    notes: Optional[list[str]] = None,
    input_info: Optional[dict] = None,
    case_table: Optional[str] = None,
    original_annotation_column: Optional[str] = None,
    original_annotation_coverage: Optional[float] = None,
    strategy_reason: str = "",
    skipped_methods: Optional[list[dict]] = None,
    failed_methods: Optional[list[dict]] = None,
    dependency_status: Optional[dict] = None,
    preflight_status: Optional[dict] = None,
    cell_calling_summary: Optional[dict] = None,
    normalization_status: str = "none",
    status: str = "annotated",
    standardization_summary: Optional[dict] = None,
) -> dict:
    """Build complete annotation manifest dict."""
    manifest = {
        "dataset_id": dataset_id,
        "benchmark_id": benchmark_id,
        "status": status,
        "case_table": case_table or "",
        "input_type": input_type,
        "input_detection": input_info or {},
        "prepared_h5ad": annotated_h5ad,
        "selected_h5ad": selected_h5ad,
        "target_cell_types": target_cell_types,
        "target_cell_type_source": "case_table.cell_type",
        "strategy_category": strategy_category,
        "strategy_reason": strategy_reason,
        "annotation_methods_used": annotation_methods_used,
        "skipped_methods": skipped_methods or [],
        "failed_methods": failed_methods or [],
        "has_original_annotation": has_original_annotation,
        "original_annotation_column": original_annotation_column or "",
        "original_annotation_coverage": original_annotation_coverage,
        "cell_calling_summary": cell_calling_summary,
        "normalization_status": normalization_status,
        "selection": selection_summary,
        "needs_manual_review": needs_manual_review,
        "case_level_target_cell_types": target_cell_types,
        "cell_type_standardization_summary": standardization_summary or {},
        "dependency_status": dependency_status,
        "preflight_status": preflight_status,
        "notes": notes or [],
    }
    return manifest


def build_failure_manifest(
    dataset_id: str,
    benchmark_id: str,
    input_type: str = "unknown",
    status: str = "failed_annotation",
    error_message: str = "",
    input_info: Optional[dict] = None,
    target_cell_types: Optional[list[str]] = None,
    preflight_status: Optional[dict] = None,
) -> dict:
    """Build a manifest for failure cases to aid debugging."""
    return {
        "dataset_id": dataset_id,
        "benchmark_id": benchmark_id,
        "status": status,
        "input_type": input_type,
        "input_detection": input_info or {},
        "target_cell_types": target_cell_types or [],
        "preflight_status": preflight_status,
        "error": error_message,
        "prepared_h5ad": "",
        "selected_h5ad": "",
        "annotation_methods_used": [],
        "notes": [error_message] if error_message else [],
    }


def save_manifest(manifest: dict, output_path: str):
    """Write manifest to JSON file."""
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(manifest, f, indent=2, default=str)
    logger.info("Manifest saved: %s", output_path)


def check_needs_manual_review(
    adata: "ad.AnnData",
    selection_summary: dict = None,
    total_annotated: int = None,
    original_target_count: int = None,
    strategy_category: str = "",
) -> tuple[bool, list[str]]:
    """Check if any cells have unresolved annotations or selection quality issues.

    Args:
        original_target_count: number of unique targets from Case-level table.
            Used for accurate unmatched ratio calculation.
        strategy_category: e.g. 'cell_line_or_in_vitro', 'normal_or_disease_tissue'.
    """
    notes = []
    needs_review = False

    if "cell_type" in adata.obs.columns:
        ct = adata.obs["cell_type"].astype(str)
        unknown_mask = ct.str.lower().isin([
            "unknown", "nan", "none", "", "null",
        ])
        n_unknown = int(unknown_mask.sum())
        if n_unknown > 0:
            notes.append(
                f"{n_unknown} cells have cell_type='Unknown' "
                f"({n_unknown / adata.n_obs:.1%})"
            )
            if n_unknown / adata.n_obs > 0.3:
                needs_review = True
                notes.append(">30% cells unresolved — manual review recommended")

        if "annotation_method" in adata.obs.columns:
            n_unresolved = int(
                (adata.obs["annotation_method"] == "Unresolved").sum()
            )
            if n_unresolved > 0:
                notes.append(
                    f"{n_unresolved} cells have annotation_method='Unresolved'"
                )

    # Selection quality checks
    if selection_summary:
        n_before = selection_summary.get("n_cells_before_selection", 0)
        n_after = selection_summary.get("n_cells_after_selection", 0)
        unmatched = selection_summary.get("unmatched_target_cell_types", [])
        matched = selection_summary.get("matched_cell_type_labels", [])
        approximate = selection_summary.get("approximate_mapped_targets", [])
        low_conf = selection_summary.get("low_confidence_targets", [])

        n_targets = original_target_count or (len(matched) + len(unmatched))

        # Determine if 100% selection is expected (e.g. cell_line with metadata)
        is_cell_line_strategy = "cell_line" in (strategy_category or "").lower()
        all_exact_metadata = (
            is_cell_line_strategy
            and "metadata" in str(
                adata.obs["annotation_method"].iloc[0]
            ).lower() if "annotation_method" in adata.obs.columns and adata.n_obs > 0
            else False
        )
        has_approximation = bool(approximate) or bool(low_conf)

        if n_before > 0:
            fraction = n_after / n_before
            if fraction < 0.01:
                needs_review = True
                notes.append(
                    f"Very low selection rate ({fraction:.2%})")
            elif fraction > 0.5:
                if all_exact_metadata and not unmatched and not has_approximation:
                    notes.append(
                        f"Full selection ({fraction:.0%}) — expected for "
                        f"{strategy_category} with metadata labeling"
                    )
                else:
                    needs_review = True
                    notes.append(
                        f"High selection rate ({fraction:.2%}) — targets may be too broad")

        if unmatched:
            notes.append(
                f"{len(unmatched)}/{n_targets} targets STILL unmatched (should be 0): {unmatched}"
            )
            needs_review = True
        else:
            notes.append(f"All {n_targets}/{n_targets} targets covered in selected.h5ad")

        if approximate:
            notes.append(
                f"{len(approximate)} targets mapped via approximation: {approximate}"
            )
            needs_review = True
        if low_conf:
            notes.append(
                f"{len(low_conf)} targets at low confidence: {low_conf}"
            )

    return needs_review, notes


def compute_selection_status(
    selection_summary: dict,
    original_target_count: int = None,
) -> str:
    """Determine selection quality status from summary."""
    n_before = selection_summary.get("n_cells_before_selection", 0)
    n_after = selection_summary.get("n_cells_after_selection", 0)
    unmatched = selection_summary.get("unmatched_target_cell_types", [])
    approximate = selection_summary.get("approximate_mapped_targets", [])
    low_conf = selection_summary.get("low_confidence_targets", [])

    if n_after == 0:
        return "no_target_cells_matched"
    if unmatched:
        return "partial_selection"
    if approximate or low_conf:
        return "full_coverage_with_approximation"
    return "selected"
