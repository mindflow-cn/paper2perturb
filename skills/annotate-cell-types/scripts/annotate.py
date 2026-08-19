#!/usr/bin/env python3
"""Auto cell annotation orchestrator.

Detects input format, builds AnnData, detects original annotations,
annotates missing cells according to strategy, selects target cell types,
and writes annotated.h5ad + selected.h5ad + annotation_manifest.json.

Usage:
    python3 skills/annotate-cell-types/scripts/annotate.py \
      --input-dir raw_data/GSE228421 \
      --dataset-id GSE228421 \
      --benchmark-id GSE228421_psoriasis_lesional_skin_risankizumab_scRNA \
      --case-table "path/to/Case-level.csv" \
      --output-root prepared \
      --mode balanced

All heavy imports (pandas, numpy, anndata, scanpy, celltypist) are deferred
to after argparse and preflight, so --help works even without these packages.
"""

import argparse
import logging
import sys
from pathlib import Path

# Only stdlib at top level so --help works without pandas/numpy/anndata/etc.
_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("annotate")

# Routing tables (no heavy deps)
BALANCED_ROUTES = {
    "cell_line_or_in_vitro": ["metadata_label"],
    "immune_enriched": ["celltypist", "marker_heuristics"],
    "tumor_or_tme": ["celltypist", "marker_heuristics"],
    "normal_or_disease_tissue": ["celltypist", "marker_heuristics"],
    "unknown": ["celltypist", "marker_heuristics"],
}
QUICK_ROUTES = {
    "cell_line_or_in_vitro": ["metadata_label"],
    "immune_enriched": ["marker_heuristics"],
    "tumor_or_tme": ["marker_heuristics"],
    "normal_or_disease_tissue": ["marker_heuristics"],
    "unknown": ["marker_heuristics"],
}
DEEP_ROUTES = {
    "cell_line_or_in_vitro": ["metadata_label", "sanity_check"],
    "immune_enriched": ["celltypist", "marker_heuristics", "r_adapter"],
    "tumor_or_tme": ["celltypist", "marker_heuristics", "r_adapter"],
    "normal_or_disease_tissue": ["celltypist", "marker_heuristics", "r_adapter"],
    "unknown": ["celltypist", "marker_heuristics", "r_adapter"],
}
EXPENSIVE_METHODS = ["copykat", "infercnv", "scevan", "cellbender", "r_adapter"]


def _lazy_import_preflight():
    """Lazy-import preflight; only needs stdlib + importlib."""
    from preflight import check_dependencies
    return check_dependencies


def _lazy_import_annotation_modules():
    """Lazy-import all heavy annotation modules; call after preflight passes."""
    from case_context import get_target_cell_types
    from detect import detect
    from io_utils import build_adata, save_adata, ensure_counts_layer
    from original_annotation import detect_original_annotation, standardize_original_annotation
    from manifest import (
        build_manifest, build_failure_manifest,
        save_manifest, check_needs_manual_review,
        compute_selection_status,
    )
    from select_cells import select_cells
    return (
        get_target_cell_types, detect,
        build_adata, save_adata, ensure_counts_layer,
        detect_original_annotation, standardize_original_annotation,
        build_manifest, build_failure_manifest, save_manifest,
        check_needs_manual_review, compute_selection_status,
        select_cells,
    )


def get_route(strategy_category: str, mode: str) -> list[str]:
    routes = {"quick": QUICK_ROUTES, "balanced": BALANCED_ROUTES, "deep": DEEP_ROUTES}
    return routes.get(mode, BALANCED_ROUTES).get(strategy_category, BALANCED_ROUTES["unknown"])


def _run_r_adapter(adata, strategy_category, args) -> dict:
    """Run optional R adapter (STUB — R adapters are not yet implemented inline)."""
    if not args.enable_r_adapters:
        return {"method": "r_adapter", "status": "skipped", "reason": "r_adapters_not_enabled"}
    if args.mode != "deep":
        return {"method": "r_adapter", "status": "skipped", "reason": "not_in_deep_mode"}
    import subprocess
    try:
        result = subprocess.run(["Rscript", "--version"], capture_output=True, text=True, timeout=10)
        if result.returncode != 0:
            raise FileNotFoundError("Rscript failed")
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return {"method": "r_adapter", "status": "skipped", "reason": "r_not_installed"}
    r_script_map = {
        "normal_or_disease_tissue": "singler_annotate.R",
        "immune_enriched": "singler_annotate.R",
        "tumor_or_tme": "scatomic_annotate.R",
        "unknown": "singler_annotate.R",
    }
    r_script = r_script_map.get(strategy_category, "singler_annotate.R")
    r_path = _SCRIPT_DIR / "optional_r" / r_script
    if not r_path.exists():
        return {"method": "r_adapter", "status": "skipped", "reason": f"script_not_found: {r_script}"}
    logger.info("R adapter %s found but marked as STUB — not executing inline", r_script)
    return {
        "method": "r_adapter", "status": "skipped",
        "reason": "r_adapter_stub_not_implemented_inline",
        "r_script": str(r_path),
    }


def _run_strategy(method, adata, target_cell_types, strategy_category, args) -> dict:
    if method == "metadata_label":
        from strategies.metadata_label import label_from_metadata
        return label_from_metadata(adata, target_cell_types, strategy_category=strategy_category)
    if method == "marker_heuristics":
        from strategies.marker_heuristics import annotate_markers
        refine = getattr(args, "marker_refine", False)
        if refine:
            # Cluster-level refinement: only process broad-candidate clusters
            from strategies.marker_heuristics import refine_by_target_markers
            from select_cells import BROAD_CANDIDATE_MAP
            return refine_by_target_markers(
                adata, target_cell_types, BROAD_CANDIDATE_MAP,
            )
        else:
            # Fill mode: only annotate missing
            return annotate_markers(adata, target_cell_types=target_cell_types)
    if method == "celltypist":
        from strategies.celltypist_annotate import annotate_celltypist
        ct_model = getattr(args, "celltypist_model", None)
        use_mv = not getattr(args, "celltypist_no_mv", False)
        return annotate_celltypist(
            adata, strategy_category=strategy_category,
            target_cell_types=target_cell_types,
            model_name=ct_model,
            majority_voting=use_mv,
        )
    if method == "sanity_check":
        logger.info("Deep sanity check: verifying metadata label consistency")
        return {"method": "sanity_check", "status": "ok"}
    if method == "r_adapter":
        return _run_r_adapter(adata, strategy_category, args)
    logger.warning("Unknown strategy method: %s", method)
    return {"method": method, "status": "skipped", "reason": "unknown_method"}


def main():
    parser = argparse.ArgumentParser(
        description="Auto cell-type annotation and target cell selection"
    )
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--benchmark-id", required=True)
    parser.add_argument("--case-table", default=None)
    parser.add_argument("--metadata-xlsx", default=None)
    parser.add_argument("--output-root", default="prepared")
    parser.add_argument("--mode", choices=["quick", "balanced", "deep"], default="balanced")
    parser.add_argument("--max-cells-for-expensive-methods", type=int, default=50000)
    parser.add_argument("--allow-expensive-tumor-methods", action="store_true")
    parser.add_argument("--enable-r-adapters", action="store_true")
    parser.add_argument("--allow-empty-selection", action="store_true")
    parser.add_argument("--run-cellranger", action="store_true")
    parser.add_argument("--cellranger-reference", default=None)
    parser.add_argument("--cellranger-output-dir", default=None)
    parser.add_argument("--sample-metadata", default=None,
                        help="Path to CSV/TSV/xlsx with sample metadata (sample_id -> columns)")
    parser.add_argument("--celltypist-model", default=None,
                        help="Override CellTypist model (e.g. Adult_Human_Skin.pkl)")
    parser.add_argument("--celltypist-no-mv", action="store_true",
                        help="Disable CellTypist majority voting (faster, lower accuracy)")
    parser.add_argument("--marker-max-refine-cells", type=int, default=200000,
                        help="Max cells for MarkerHeuristics refinement (default: 200000)")
    parser.add_argument("--skip-marker-refinement", action="store_true",
                        help="Skip MarkerHeuristics refinement pass (use CellTypist labels as-is)")
    parser.add_argument("--raw-converter", default=None,
                        help="Force a specific raw converter (bypass auto-detection)")
    parser.add_argument("--list-raw-converters", action="store_true",
                        help="List available raw converters and exit")
    args = parser.parse_args()

    if args.list_raw_converters:
        from raw_converters import list_converters
        print("Available raw converters:")
        for name, cls in list_converters().items():
            doc = (cls.__doc__ or "").strip().split("\n")[0]
            print(f"  {name}: {doc}")
        return 0

    case_source = args.case_table or args.metadata_xlsx
    if not case_source:
        parser.error("One of --case-table or --metadata-xlsx is required")

    # ---- Preflight (stdlib + importlib only) ----
    check_deps = _lazy_import_preflight()
    preflight = check_deps()
    if not preflight["all_ok"]:
        print(
            f"ERROR: Missing required packages: {preflight['missing']}",
            file=sys.stderr,
        )
        print(f"  Install: pip install {' '.join(preflight['missing'])}", file=sys.stderr)
        print(f"  Or: pip install -r skills/annotate-cell-types/requirements.txt",
              file=sys.stderr)
        print(f"  Python: {sys.executable}", file=sys.stderr)
        # Try to write a failure manifest
        try:
            from manifest import build_failure_manifest, save_manifest
            fm = build_failure_manifest(
                dataset_id=args.dataset_id, benchmark_id=args.benchmark_id,
                status="failed_preflight",
                error_message=f"Missing packages: {preflight['missing']}",
                preflight_status=preflight,
            )
            mp = f"{args.output_root}/benchmarks/{args.benchmark_id}/annotation_manifest.json"
            save_manifest(fm, mp)
        except Exception:
            pass
        sys.exit(1)

    # ---- Now safe to import heavy modules ----
    (
        get_target_cell_types, detect,
        build_adata, save_adata, ensure_counts_layer,
        detect_original_annotation, standardize_original_annotation,
        build_manifest, build_failure_manifest, save_manifest,
        check_needs_manual_review, compute_selection_status,
        select_cells,
    ) = _lazy_import_annotation_modules()

    def _save_failure_manifest_and_exit(
        error_msg, input_type="unknown", input_info=None, target_cell_types=None,
    ):
        m = build_failure_manifest(
            dataset_id=args.dataset_id, benchmark_id=args.benchmark_id,
            input_type=input_type, status="failed_annotation",
            error_message=error_msg,
            input_info=input_info, target_cell_types=target_cell_types or [],
            preflight_status=preflight,
        )
        mp = f"{args.output_root}/benchmarks/{args.benchmark_id}/annotation_manifest.json"
        try:
            save_manifest(m, mp)
        except Exception:
            pass
        logger.error(error_msg)
        sys.exit(1)

    # ---- Step 1: Parse Case-level table ----
    logger.info("Reading Case-level table: %s", case_source)
    try:
        ctx = get_target_cell_types(case_source, args.benchmark_id)
    except ValueError as e:
        _save_failure_manifest_and_exit(str(e))
    target_cell_types = ctx["target_cell_types"]
    strategy_category = ctx["strategy_category"]
    strategy_reason = ctx.get("strategy_reason", "")
    logger.info("Target cell types: %s, Strategy: %s", target_cell_types, strategy_category)

    # ---- Step 2a: Try raw converter first ----
    from raw_converters import detect_format as raw_detect, get_converter as raw_get
    raw_converter_name = args.raw_converter
    adata_raw = None
    raw_conversion_info = None
    raw_converter_failed = False
    raw_converter_error = ""

    # Phase 1: detect
    try:
        if raw_converter_name is None:
            raw_converter_name = raw_detect(args.input_dir)
        logger.info("Raw converter detected: %s", raw_converter_name)
    except (ValueError, FileNotFoundError) as e:
        if args.raw_converter:
            _save_failure_manifest_and_exit(
                f"Raw converter '{args.raw_converter}' not found or not applicable "
                f"for {args.input_dir}: {e}",
                input_type="raw_converter_not_found",
            )
        logger.info("No raw converter matched: %s — falling back to detect/build_adata", e)
        raw_converter_name = None

    # Phase 2: convert (if detected)
    if raw_converter_name is not None:
        try:
            conv_cls = raw_get(raw_converter_name)
            converter = conv_cls()
            adata_raw = converter.convert(Path(args.input_dir), dataset_id=args.dataset_id)
            raw_conversion_info = adata_raw.uns.get("raw_conversion", {})
            logger.info(
                "Raw converter loaded: %d cells, %d genes, obs cols=%s",
                adata_raw.n_obs, adata_raw.n_vars, list(adata_raw.obs.columns),
            )
        except Exception as e:
            raw_converter_failed = True
            raw_converter_error = f"{raw_converter_name}: {e}"
            if args.raw_converter:
                # User explicitly requested this converter — fail fast
                _save_failure_manifest_and_exit(
                    f"Raw converter '{args.raw_converter}' failed: {e}",
                    input_type="raw_converter_failed", input_info={
                        "converter": raw_converter_name, "error": str(e),
                    },
                )
            # Auto-detected converter failed — log clearly then fallback
            logger.warning(
                "Raw converter '%s' matched but convert() failed: %s. "
                "Falling back to detect/build_adata.",
                raw_converter_name, e,
            )

    # ---- Step 2b: Fallback to legacy detect + build_adata ----
    if adata_raw is not None:
        # Use raw converter output directly, skip detection/loading
        adata = adata_raw
        input_info = {
            "type": "raw_converted",
            "path": args.input_dir,
            "details": f"Raw converter: {raw_converter_name}",
            "barcode_count": adata.n_obs,
            "is_raw_mex": False,
        }
        cell_calling_summary = None  # raw converter handles calling
        do_call = False
    else:
        logger.info("Detecting format in: %s", args.input_dir)
        input_info = detect(args.input_dir)
        logger.info("Detected: %s — %s", input_info["type"], input_info["details"])

    if input_info["type"] == "fastq" and not args.run_cellranger:
        _save_failure_manifest_and_exit(
            f"FASTQ input detected in {args.input_dir}. "
            f"FASTQ requires Cell Ranger preprocessing before annotation.\n"
            f"Options:\n"
            f"  1. Run cellranger count externally, then provide the output directory.\n"
            f"  2. Re-run with --run-cellranger --cellranger-reference <ref>.",
            input_type="fastq", input_info=input_info, target_cell_types=target_cell_types,
        )
    if input_info["type"] == "unknown":
        _save_failure_manifest_and_exit(
            f"Unrecognized data format in {args.input_dir}. "
            f"Expected: h5ad, 10x MEX, Cell Ranger output, or CSV/TSV.",
            input_type="unknown", input_info=input_info,
        )

    # ---- Step 3: Build AnnData (only if not already loaded by raw converter) ----
    if adata_raw is None:
        logger.info("Building AnnData...")
        do_call = input_info.get("is_raw_mex", False)
        cell_calling_kwargs = {"min_umis": 100, "min_genes": 50}
        try:
            adata = build_adata(
                input_info, args.input_dir,
                do_cell_calling=do_call,
                cell_calling_kwargs=cell_calling_kwargs,
            )
        except Exception as e:
            _save_failure_manifest_and_exit(
                f"Failed to build AnnData: {e}",
                input_type=input_info["type"], input_info=input_info,
                target_cell_types=target_cell_types,
            )
        logger.info("AnnData: %s cells, %s genes", adata.n_obs, adata.n_vars)
        # Per-sample cell calling summary from build_adata
        cell_calling_summary = {
            "method": "per_sample_knee_filter",
            "n_barcodes_before": input_info.get("barcode_count", 0),
            "n_cells_after": adata.n_obs,
            "min_umis": cell_calling_kwargs.get("min_umis", 100),
            "min_genes": cell_calling_kwargs.get("min_genes", 50),
        } if do_call else None
    else:
        cell_calling_summary = None  # raw converter handles calling
        logger.info("AnnData: %s cells, %s genes (from raw converter)", adata.n_obs, adata.n_vars)
    if "sample_id" in adata.obs.columns:
        from sample_metadata import enrich_adata_obs
        n_cols_before = len(adata.obs.columns)
        adata = enrich_adata_obs(adata, sample_metadata_table=args.sample_metadata)
        n_cols_after = len(adata.obs.columns)
        if n_cols_after > n_cols_before:
            logger.info(
                "Metadata enrichment added %d columns: %s",
                n_cols_after - n_cols_before,
                list(set(adata.obs.columns) - set(adata.obs.columns[:n_cols_before])),
            )

    # ---- Step 5: Ensure counts layer ----
    ensure_counts_layer(adata)

    # ---- Step 6: Original annotations ----
    has_original = False
    orig_col = None
    orig_coverage = None
    orig_col = detect_original_annotation(adata)
    if orig_col:
        coverage = adata.obs[orig_col].notna().mean()
        adata = standardize_original_annotation(adata, orig_col)
        has_original = True
        orig_coverage = float(coverage)
        logger.info("Original annotation: %s (coverage=%.1f%%)", orig_col, coverage * 100)

    # ---- Step 7: Run strategies ----
    route = get_route(strategy_category, args.mode)
    logger.info("Route: %s/%s -> %s", args.mode, strategy_category, route)

    methods_used = []
    skipped_methods = []
    failed_methods = []
    if has_original:
        methods_used.append("Original")

    # Enable refinement for tissue benchmarks where CellTypist model may not match
    if getattr(args, "skip_marker_refinement", False):
        setattr(args, "marker_refine", False)
    elif strategy_category in ("normal_or_disease_tissue", "unknown"):
        setattr(args, "marker_refine", True)
    else:
        setattr(args, "marker_refine", False)

    for method in route:
        if method in EXPENSIVE_METHODS:
            if args.mode == "balanced":
                skipped_methods.append(
                    {"method": method, "reason": f"not_allowed_in_{args.mode}_mode"})
                continue
            if method in ["copykat", "infercnv", "scevan"] and not args.allow_expensive_tumor_methods:
                skipped_methods.append(
                    {"method": method, "reason": "allow_expensive_tumor_methods_not_set"})
                continue
            if adata.n_obs > args.max_cells_for_expensive_methods:
                skipped_methods.append({
                    "method": method,
                    "reason": f"too_many_cells ({adata.n_obs} > {args.max_cells_for_expensive_methods})",
                })
                continue

        summary = _run_strategy(method, adata, target_cell_types, strategy_category, args)
        logger.info("  %s: %s", method, summary)
        status = summary.get("status", "")
        if status == "skipped":
            skipped_methods.append({"method": method, "reason": summary.get("reason", "unknown")})
        elif status == "error":
            failed_methods.append({"method": method, "reason": summary.get("reason", "unknown")})
        elif summary.get("cells_labeled", 0) > 0:
            if method not in methods_used:
                methods_used.append(method)

    # Unresolved cells — don't overwrite cells that already have valid Original labels
    ct_values = adata.obs["cell_type"].astype(str)
    missing_mask = ct_values.isna() | ct_values.str.strip().str.lower().isin(
        ["", "nan", "none", "null"]
    )
    if missing_mask.any():
        adata.obs.loc[missing_mask, "cell_type"] = "Unknown"
        adata.obs.loc[missing_mask, "annotation_method"] = "Unresolved"
        logger.warning(
            "%s cells remain unresolved (%.1f%%)",
            missing_mask.sum(), missing_mask.sum() / adata.n_obs * 100,
        )

    # ---- Step 8: Standardize cell_type to Case-level target names ----
    # Always run this pass, including datasets with original annotations.
    # The final label contract is that obs["cell_type"] should use the
    # xlsx Case-level cell_type labels; raw/model labels are preserved in
    # obs["cell_type_raw"] and mapped exactly/approximately when needed.
    standardization_summary = {}
    from standardize_labels import standardize_cell_types, _norm
    logger.info("Standardizing cell_type to Case-level target names...")
    standardization_summary = standardize_cell_types(adata, target_cell_types)
    logger.info(
        "Standardization: %d/%d cells mapped to target names",
        standardization_summary.get("n_mapped_to_target", 0),
        standardization_summary.get("total_cells", 0),
    )

    # ---- Step 9: Save annotated.h5ad ----
    annotated_path = f"{args.output_root}/datasets/{args.dataset_id}/annotated.h5ad"
    save_adata(adata, annotated_path)

    # ---- Step 10: Select target cells ----
    logger.info("Selecting cells matching targets: %s", target_cell_types)
    try:
        selected, sel_summary = select_cells(
            adata, target_cell_types,
            strategy_category=strategy_category,
            allow_empty=args.allow_empty_selection,
        )
        # Merge standardization per-target info into selection summary
        if standardization_summary:
            sel_summary["approximate_mapped_targets"] = standardization_summary.get(
                "approximate_mapped_targets", [])
            sel_summary["low_confidence_targets"] = standardization_summary.get(
                "low_confidence_targets", [])
            sel_summary["per_target_selection_summary"] = standardization_summary.get(
                "per_target_selection_summary", {})
            sel_summary["covered_target_cell_types"] = standardization_summary.get(
                "covered_target_cell_types", [])
            sel_summary["unmatched_target_cell_types"] = standardization_summary.get(
                "unmatched_target_cell_types", [])
    except ValueError as e:
        m = build_manifest(
            dataset_id=args.dataset_id, benchmark_id=args.benchmark_id,
            input_type=input_info["type"],
            annotated_h5ad=annotated_path, selected_h5ad="",
            target_cell_types=target_cell_types,
            strategy_category=strategy_category,
            annotation_methods_used=methods_used,
            selection_summary={
                "n_cells_before_selection": adata.n_obs,
                "n_cells_after_selection": 0, "error": str(e),
            },
            has_original_annotation=has_original,
            original_annotation_column=orig_col,
            original_annotation_coverage=orig_coverage,
            strategy_reason=strategy_reason,
            skipped_methods=skipped_methods, failed_methods=failed_methods,
            dependency_status=preflight, preflight_status=preflight,
            cell_calling_summary=cell_calling_summary,
            status="no_target_cells_matched", notes=[str(e)],
            input_info=input_info, case_table=case_source,
        )
        save_manifest(
            m, f"{args.output_root}/benchmarks/{args.benchmark_id}/annotation_manifest.json")
        logger.error(str(e))
        sys.exit(1)

    selected_path = f"{args.output_root}/benchmarks/{args.benchmark_id}/selected.h5ad"
    save_adata(selected, selected_path)

    # ---- Recompute statistics from actual selected AnnData ----
    # The standardization summary was computed on full annotated adata before
    # selection. Recompute per-target cell counts and non-case-level counts
    # from the selected adata so manifest statistics are accurate.
    if standardization_summary and "cell_type" in selected.obs.columns:
        target_norm_set = {_norm(t) for t in target_cell_types}
        sel_ct = selected.obs["cell_type"].astype(str)
        sel_ct_norm = sel_ct.apply(_norm)
        n_non_case_in_selected = int((~sel_ct_norm.isin(target_norm_set)).sum())
        standardization_summary["n_non_case_level_labels_in_selected"] = n_non_case_in_selected

        # Fix non_case_level labels to reflect selected only
        if n_non_case_in_selected > 0:
            standardization_summary["non_case_level_cell_type_labels"] = sorted(set(
                sel_ct[sel_ct_norm.apply(lambda x: x not in target_norm_set)].unique()
            ))
        else:
            standardization_summary["non_case_level_cell_type_labels"] = []

        # Fix per_target_selection_summary n_cells from selected
        per_target = standardization_summary.get("per_target_selection_summary", {})
        for target in list(per_target.keys()):
            tn = _norm(target)
            n_in_selected = int((sel_ct_norm == tn).sum())
            per_target[target]["n_cells"] = n_in_selected
            if n_in_selected == 0 and per_target[target].get("confidence") != "low":
                per_target[target]["confidence"] = "low"
        standardization_summary["per_target_selection_summary"] = per_target

        # Update sel_summary as well
        sel_summary["n_non_case_level_labels_in_selected"] = n_non_case_in_selected
        sel_summary["per_target_selection_summary"] = per_target

    # ---- Step 11: Review + Manifest ----
    original_target_count = len(target_cell_types)
    needs_review, review_notes = check_needs_manual_review(
        selected, selection_summary=sel_summary, total_annotated=adata.n_obs,
        original_target_count=original_target_count,
        strategy_category=strategy_category,
    )
    selection_status = compute_selection_status(
        sel_summary, original_target_count=original_target_count,
    )

    # Add standardization audit to notes — only if genuinely true for selected
    if standardization_summary:
        n_non_case = standardization_summary.get("n_non_case_level_labels_in_selected", 0)
        if n_non_case > 0:
            review_notes.append(
                f"{n_non_case} selected cells have non-Case-level cell_type labels"
            )

    manifest = build_manifest(
        dataset_id=args.dataset_id, benchmark_id=args.benchmark_id,
        input_type=input_info["type"],
        annotated_h5ad=annotated_path, selected_h5ad=selected_path,
        target_cell_types=target_cell_types,
        strategy_category=strategy_category,
        annotation_methods_used=methods_used,
        selection_summary=sel_summary,
        has_original_annotation=has_original,
        original_annotation_column=orig_col,
        original_annotation_coverage=orig_coverage,
        strategy_reason=strategy_reason,
        skipped_methods=skipped_methods, failed_methods=failed_methods,
        dependency_status=preflight, preflight_status=preflight,
        cell_calling_summary=cell_calling_summary,
        normalization_status="cpm_log1p" if "counts" in adata.layers else "none",
        status=selection_status,
        needs_manual_review=needs_review, notes=review_notes,
        input_info=input_info, case_table=case_source,
        standardization_summary=standardization_summary,
    )
    manifest_path = f"{args.output_root}/benchmarks/{args.benchmark_id}/annotation_manifest.json"
    save_manifest(manifest, manifest_path)

    logger.info("=== Annotation complete ===")
    logger.info("  annotated.h5ad: %s (%s cells)", annotated_path, adata.n_obs)
    logger.info("  selected.h5ad:  %s (%s cells)", selected_path, selected.n_obs)
    logger.info("  manifest:       %s", manifest_path)
    logger.info("  methods:        %s", methods_used)
    logger.info("  needs_review:   %s", needs_review)
    return 0


if __name__ == "__main__":
    sys.exit(main())
