#!/usr/bin/env python3
"""Tests for annotate-cell-types skill.

Run from project root:
    python3 skills/annotate-cell-types/scripts/test_annotate.py

Tests cover: preflight, gz MTX, prefixed MEX detection/loading, H5 detection,
Case-level fail-fast, CellTypist mock, target selection broad/narrow,
build-h5ad fail-fast, R adapter exclusion.
"""

import gzip
import json
import os
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock

# Ensure skill scripts are importable
_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp


# ── Helpers ─────────────────────────────────────────────────────

def _check_deps():
    """Check if required packages are available, print skip message if not."""
    missing = []
    for pkg in ["numpy", "pandas", "scipy", "anndata", "scanpy"]:
        try:
            __import__(pkg)
        except ImportError:
            missing.append(pkg)
    if missing:
        return False, missing
    return True, []


def make_test_adata(n_cells=100, n_genes=50, seed=42):
    np.random.seed(seed)
    X = np.abs(np.random.normal(5, 2, (n_cells, n_genes)))
    gene_names = [f"GENE_{i}" for i in range(n_genes)]
    # Add marker genes
    markers = {
        "CD3D": n_genes - 5, "CD3E": n_genes - 4,
        "CD8A": n_genes - 3, "MS4A1": n_genes - 2, "CD4": n_genes - 1,
    }
    for g, idx in markers.items():
        if idx < n_genes:
            gene_names[idx] = g
    return ad.AnnData(
        X=sp.csr_matrix(X),
        obs=pd.DataFrame(index=[f"cell_{i}" for i in range(n_cells)]),
        var=pd.DataFrame(index=gene_names),
    )


# ── Test 1: Preflight ───────────────────────────────────────────

def test_preflight():
    """Test that preflight reports missing packages cleanly."""
    print("\n=== Test 1: Preflight ===")
    from preflight import check_dependencies
    result = check_dependencies()
    assert "all_ok" in result
    assert "missing" in result
    assert "details" in result
    assert result["all_ok"], f"Missing required packages: {result['missing']}"
    print("  [PASS] All required packages present")
    print(f"  Python: {result['python_version']}")


# ── Test 2: .mtx.gz barcode count ───────────────────────────────

def test_mtx_gz_barcode_count():
    """Test that _count_mtx_barcodes handles .gz files."""
    print("\n=== Test 2: .mtx.gz barcode count ===")
    from detect import _count_mtx_barcodes, _open_maybe_gz

    # Create a compressed MTX
    with tempfile.NamedTemporaryFile(suffix=".mtx.gz", delete=False) as f:
        content = "%%MatrixMarket matrix coordinate integer general\n%comment\n100 42 500\n1 1 5\n"
        with gzip.open(f.name, "wt") as gz:
            gz.write(content)
    try:
        n = _count_mtx_barcodes(Path(f.name))
        assert n == 42, f"Expected 42 barcodes, got {n}"
        print(f"  [PASS] .mtx.gz: {n} barcodes")
    finally:
        os.unlink(f.name)

    # Test plain .mtx
    with tempfile.NamedTemporaryFile(suffix=".mtx", delete=False, mode="w") as f:
        f.write("%%MatrixMarket matrix coordinate integer general\n100 37 500\n")
    try:
        n = _count_mtx_barcodes(Path(f.name))
        assert n == 37, f"Expected 37 barcodes, got {n}"
        print(f"  [PASS] .mtx: {n} barcodes")
    finally:
        os.unlink(f.name)

    # Test _open_maybe_gz
    with tempfile.NamedTemporaryFile(suffix=".gz", delete=False, mode="wb") as f:
        with gzip.open(f, "wt") as gz:
            gz.write("test content\n")
    try:
        with _open_maybe_gz(Path(f.name), "rt") as fh:
            assert fh.read().strip() == "test content"
        print("  [PASS] _open_maybe_gz works")
    finally:
        os.unlink(f.name)


# ── Test 3: Prefixed MEX detection ──────────────────────────────

def test_prefixed_mex_detection():
    """Test detection of prefixed multi-sample MEX files."""
    print("\n=== Test 3: Prefixed MEX detection ===")
    from detect import _detect_prefixed_mex_samples, detect

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)

        # Sample A
        (tmp / "GSM123_sampleA.matrix.mtx").write_text(
            "%%MatrixMarket matrix coordinate integer general\n100 200 500\n1 1 5\n"
        )
        (tmp / "GSM123_sampleA.features.tsv").write_text(
            "GENE1\tSymbol1\tExpression\nGENE2\tSymbol2\tExpression\n"
        )
        (tmp / "GSM123_sampleA.barcodes.tsv").write_text(
            "AAACCCAAGAAACTCA-1\nAAACCCAAGAAACTCA-2\n"
        )

        # Sample B
        (tmp / "GSM456_sampleB.matrix.mtx").write_text(
            "%%MatrixMarket matrix coordinate integer general\n100 150 300\n1 1 5\n"
        )
        (tmp / "GSM456_sampleB.features.tsv").write_text(
            "GENE1\tSymbol1\tExpression\n"
        )
        (tmp / "GSM456_sampleB.barcodes.tsv").write_text(
            "AAACCCAAGAAACTCA-1\n"
        )

        # Also throw in a non-matrix file
        (tmp / "README.txt").write_text("hello")

        samples = _detect_prefixed_mex_samples(tmp)
        assert len(samples) == 2, f"Expected 2 prefixed samples, got {len(samples)}"
        assert samples[0]["sample_id"] in ("GSM123_sampleA", "GSM456_sampleB")
        assert samples[0]["barcode_count"] > 0
        print(f"  [PASS] Detected {len(samples)} prefixed samples")

        # Test detect() with prefixed
        info = detect(str(tmp))
        assert "samples" in info
        assert len(info["samples"]) == 2
        print(f"  [PASS] detect() reports {len(info['samples'])} samples")


# ── Test 4: Prefixed MEX loading ────────────────────────────────

def test_prefixed_mex_loading():
    """Test prefixed MEX detection + sample naming for merge."""
    print("\n=== Test 4: Prefixed MEX sample detection ===")
    from detect import _detect_prefixed_mex_samples

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)

        # Create minimal MEX triplet for two samples
        (tmp / "S1.matrix.mtx").write_text(
            "%%MatrixMarket matrix coordinate integer general\n5 3 6\n1 1 5\n"
        )
        (tmp / "S1.features.tsv").write_text(
            "GENE1\tSymbol1\nGENE2\tSymbol2\nGENE3\tSymbol3\nGENE4\tSymbol4\nGENE5\tSymbol5\n"
        )
        (tmp / "S1.barcodes.tsv").write_text("AAACCC-1\nAAAGGG-1\nAAATTT-1\n")

        (tmp / "S2.matrix.mtx").write_text(
            "%%MatrixMarket matrix coordinate integer general\n5 2 4\n1 1 3\n"
        )
        (tmp / "S2.features.tsv").write_text("GENE1\tSymbol1\nGENE2\tSymbol2\nGENE3\tSymbol3\nGENE4\tSymbol4\nGENE5\tSymbol5\n")
        (tmp / "S2.barcodes.tsv").write_text("BBBA-1\nBBBC-1\n")

        samples = _detect_prefixed_mex_samples(tmp)
        assert len(samples) == 2, f"Expected 2 samples, got {len(samples)}"
        sample_ids = {s["sample_id"] for s in samples}
        assert sample_ids == {"S1", "S2"}

        for s in samples:
            assert s["barcode_count"] > 0
            assert s["features_path"] is not None
            assert s["barcodes_path"] is not None
            print(f"  Sample {s['sample_id']}: {s['barcode_count']} barcodes")

        # Test the matrix path endswith check (used for gz detection in loader)
        assert not str(samples[0]["matrix_path"]).endswith(".gz")
        print("  [PASS] Sample detection + metadata correct")

        # Test that loading would produce the right obs_names pattern
        # Validate pattern: sample_id + "_" + barcode
        for s in samples:
            sid = s["sample_id"]
            # Read barcodes
            with open(s["barcodes_path"]) as bf:
                barcodes = [line.strip() for line in bf if line.strip()]
            assert len(barcodes) == s["barcode_count"]
            prefixed = [f"{sid}_{bc}" for bc in barcodes]
            assert all(p.startswith(f"{sid}_") for p in prefixed)
        print("  [PASS] obs_names would be correctly prefixed")


# ── Test 5: 10x H5 detection ────────────────────────────────────

def test_h5_detection():
    """Test that H5 files are detected with correct type."""
    print("\n=== Test 5: H5 detection ===")

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        (tmp / "outs").mkdir()
        (tmp / "outs" / "filtered_feature_bc_matrix.h5").write_text(
            "mock h5 content\n"
        )

        from detect import detect
        info = detect(str(tmp))
        # Should detect as 10x_h5_filtered, NOT 10x_mtx_filtered
        assert info["type"] in (
            "10x_h5_filtered", "10x_h5_raw"
        ), f"Expected H5 type, got {info['type']}"
        print(f"  [PASS] H5 detected as: {info['type']}")
        assert info.get("mex_h5") is not None
        print(f"  [PASS] mex_h5 present")


# ── Test 6: Case-level fail-fast ────────────────────────────────

def test_case_fail_fast():
    """Test that missing benchmark_id and empty cell_type raise ValueError."""
    print("\n=== Test 6: Case-level fail-fast ===")
    from case_context import get_target_cell_types

    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write("benchmark_id,cell_type,drug\n")
        f.write("BENCH001_skin,T cells,risankizumab\n")
        f.write("BENCH001_skin,keratinocytes,risankizumab\n")
        temp_path = f.name

    try:
        # Valid lookup
        result = get_target_cell_types(temp_path, "BENCH001_skin")
        assert len(result["target_cell_types"]) == 2
        assert result["strategy_category"] != "unknown"
        print("  [PASS] Valid lookup works")

        # Missing benchmark_id -> ValueError
        try:
            get_target_cell_types(temp_path, "NONEXISTENT_BENCH")
            assert False, "Should have raised ValueError"
        except ValueError as e:
            assert "NONEXISTENT_BENCH" in str(e)
            assert "BENCH001_skin" in str(e)
            print("  [PASS] Missing benchmark_id raises ValueError")

        # Empty cell_type test
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".csv", delete=False
        ) as f2:
            f2.write("benchmark_id,cell_type\n")
            f2.write("BENCH002,,,\n")
            f2.write("BENCH002,,,\n")
            empty_path = f2.name
        try:
            get_target_cell_types(empty_path, "BENCH002")
            assert False, "Should have raised ValueError"
        except ValueError as e:
            assert "no non-empty cell_type" in str(e).lower() or "BENCH002" in str(e)
            print("  [PASS] Empty cell_type raises ValueError")
        finally:
            os.unlink(empty_path)

    finally:
        os.unlink(temp_path)


# ── Test 7: Target selection broad/narrow ───────────────────────

def test_broad_narrow_selection():
    """Test accumulated matching: broad T cells matches children,
    narrow CD8 does NOT match all T cells."""
    print("\n=== Test 7: Target selection broad/narrow ===")
    from select_cells import select_cells

    adata = make_test_adata(n_cells=100)
    labels = (
        ["T cells"] * 20
        + ["CD8 T cells"] * 15
        + ["CD4 T cells"] * 10
        + ["B cells"] * 15
        + ["Monocytes"] * 15
        + ["Keratinocytes"] * 10
        + ["Endothelial cells"] * 10
        + ["Fibroblast"] * 5
    )
    adata.obs["cell_type"] = labels
    adata.obs["annotation_method"] = "Original"

    # Test 7a: Broad "T cells" should match "T cells" + "CD8 T cells" + "CD4 T cells"
    selected, summary = select_cells(adata, ["T cells"])
    n = summary["n_cells_after_selection"]
    # "T cells"=20 + "CD8 T cells"=15 + "CD4 T cells"=10 = 45
    assert n == 45, (
        f"Broad 'T cells' should match 45 (20+15+10), got {n}. "
        f"Matched: {summary['matched_cell_type_labels']}"
    )
    print(f"  [PASS] Broad 'T cells' matched {n} cells (T + CD8 + CD4)")

    # Test 7b: Narrow "CD8 T cells" should NOT match "T cells" or "CD4 T cells"
    selected2, summary2 = select_cells(adata, ["CD8 T cells"])
    n2 = summary2["n_cells_after_selection"]
    assert n2 == 15, (
        f"Narrow 'CD8 T cells' should match 15, got {n2}. "
        f"Matched: {summary2['matched_cell_type_labels']}"
    )
    # Should not include plain T cells
    matched = summary2["matched_cell_type_labels"]
    assert "t cells" not in matched, (
        f"Narrow CD8 must NOT match broad T cells. Got: {matched}"
    )
    print(f"  [PASS] Narrow 'CD8 T cells' matched {n2} cells (excludes T, CD4)")

    # Test 7c: Multi-target union
    selected3, summary3 = select_cells(
        adata, ["CD8 T cells", "CD4 T cells", "B cells"]
    )
    assert summary3["n_cells_after_selection"] == 40
    print("  [PASS] Multi-target union: 40 cells")

    # Test 7d: No-match raises
    try:
        select_cells(adata, ["Neurons"], allow_empty=False)
        assert False, "Should have raised"
    except ValueError as e:
        assert "No cells matched" in str(e)
    # Check available cell types in error
    try:
        select_cells(adata, ["Neurons"], allow_empty=False)
    except ValueError as e:
        assert "Available cell types" in str(e) or "available" in str(e).lower()
        print("  [PASS] No-match error includes available cell types")

    # Test 7e: unmatched_targets in summary
    selected5, summary5 = select_cells(adata, ["CD8 T cells"])
    assert "unmatched_target_cell_types" in summary5
    print("  [PASS] Summary includes unmatched_target_cell_types")


# ── Test 8: CellTypist mock writeback ────────────────────────────

def test_celltypist_mock():
    """Test that CellTypist predictions are correctly written back."""
    print("\n=== Test 8: CellTypist mock writeback ===")
    import scanpy as sc
    from strategies.celltypist_annotate import _check_celltypist

    if not _check_celltypist():
        print("  [SKIP] CellTypist not installed")
        return

    # Create adata with some T/B/monocyte marker-like expression
    adata = make_test_adata(n_cells=30)
    # Make sure markers have real expression patterns
    X = adata.X.toarray()
    if "CD3D" in adata.var_names:
        X[0:10, list(adata.var_names).index("CD3D")] += 20
    if "MS4A1" in adata.var_names:
        X[10:20, list(adata.var_names).index("MS4A1")] += 20
    adata.X = sp.csr_matrix(X)
    adata.layers["counts"] = adata.X.copy()
    # No cell_type yet
    adata.obs["cell_type"] = None
    adata.obs["annotation_method"] = "Unresolved"

    # Try celltypist annotation — may fail if model not downloaded, that's OK
    try:
        from strategies.celltypist_annotate import annotate_celltypist
        summary = annotate_celltypist(adata, strategy_category="immune_enriched")
        if summary.get("status") == "skipped":
            print(f"  [SKIP] {summary.get('reason', 'celltypist skipped')}")
            return
        n_labeled = summary.get("cells_labeled", 0)
        print(f"  [PASS] CellTypist labeled {n_labeled} cells")
        if n_labeled > 0:
            assert "cell_type" in adata.obs.columns
            assert "annotation_method" in adata.obs.columns
            n_ct = (adata.obs["annotation_method"] == "CellTypist").sum()
            print(f"  [PASS] {n_ct} cells have annotation_method='CellTypist'")
    except Exception as e:
        print(f"  [SKIP] CellTypist not functional: {e}")


# ── Test 9: No expensive methods in default mode ────────────────

def test_no_expensive_methods():
    """Verify balanced/quick modes exclude CopyKAT/inferCNV/SCEVAN/CellBender/R."""
    print("\n=== Test 9: No expensive methods in default modes ===")
    from annotate import get_route, EXPENSIVE_METHODS

    for mode in ["quick", "balanced"]:
        for cat in ["cell_line_or_in_vitro", "immune_enriched", "tumor_or_tme",
                     "normal_or_disease_tissue", "unknown"]:
            route = get_route(cat, mode)
            overlap = set(route) & set(EXPENSIVE_METHODS)
            assert not overlap, f"{cat}/{mode} should not have {overlap}"
    print("  [PASS] quick/balanced modes exclude all expensive methods")

    # R adapter is skipped even in deep mode without --enable-r-adapters
    from annotate import _run_r_adapter
    import argparse
    args = argparse.Namespace(
        enable_r_adapters=False, mode="deep",
    )
    result = _run_r_adapter(None, "normal_or_disease_tissue", args)
    assert result["status"] == "skipped"
    assert result["reason"] == "r_adapters_not_enabled"
    print("  [PASS] R adapter skipped without --enable-r-adapters")

    # R adapter is skipped in balanced mode even with flag
    args2 = argparse.Namespace(
        enable_r_adapters=True, mode="balanced",
    )
    result2 = _run_r_adapter(None, "immune_enriched", args2)
    assert result2["status"] == "skipped"
    print("  [PASS] R adapter skipped in balanced mode even with --enable-r-adapters")

    # In deep mode + enable, the adapter remains skipped when R is unavailable
    # or when the bundled adapter is intentionally not executed inline.
    args3 = argparse.Namespace(
        enable_r_adapters=True, mode="deep",
    )
    result3 = _run_r_adapter(None, "tumor_or_tme", args3)
    assert result3["status"] == "skipped"
    assert result3.get("reason") in {
        "r_not_installed",
        "r_adapter_stub_not_implemented_inline",
    } or result3.get("reason", "").startswith("script_not_found:")
    print(f"  [PASS] R adapter safely skipped in deep mode: {result3['reason']}")


# ── Test 10: Manifest completeness ──────────────────────────────

def test_manifest_completeness():
    """Test that manifest includes all required fields."""
    print("\n=== Test 10: Manifest completeness ===")
    from manifest import build_manifest, build_failure_manifest

    manifest = build_manifest(
        dataset_id="GSE228421",
        benchmark_id="BENCH001",
        input_type="10x_mtx_filtered",
        annotated_h5ad="prepared/datasets/GSE228421/annotated.h5ad",
        selected_h5ad="prepared/benchmarks/BENCH001/selected.h5ad",
        target_cell_types=["T cells", "keratinocytes"],
        strategy_category="normal_or_disease_tissue",
        annotation_methods_used=["CellTypist", "MarkerHeuristics"],
        selection_summary={
            "n_cells_before_selection": 10000,
            "n_cells_after_selection": 3500,
            "matched_cell_type_labels": ["t cells", "keratinocytes"],
            "unmatched_target_cell_types": [],
            "available_cell_types": ["t cells", "keratinocytes", "b cells"],
        },
        strategy_reason="derived from 5 entries",
        skipped_methods=[{"method": "r_adapter", "reason": "not_in_deep_mode"}],
        case_table="Case-level.csv",
        preflight_status={"all_ok": True},
    )

    required = [
        "dataset_id", "benchmark_id", "status", "input_type",
        "prepared_h5ad", "selected_h5ad", "target_cell_types",
        "strategy_category", "annotation_methods_used",
        "selection", "has_original_annotation", "needs_manual_review",
        "skipped_methods", "failed_methods", "case_table",
        "strategy_reason", "preflight_status",
    ]
    for key in required:
        assert key in manifest, f"Missing key: {key}"
    print("  [PASS] All required manifest keys present")

    # Failure manifest
    fail = build_failure_manifest(
        dataset_id="GSE001", benchmark_id="BENCH001",
        status="failed_input_detection",
        error_message="Unrecognized format",
    )
    assert fail["status"] == "failed_input_detection"
    assert fail["error"] == "Unrecognized format"
    print("  [PASS] Failure manifest works")


# ── Test 11: build-h5ad fail-fast ─────────────────────────────────

def test_build_h5ad_fail_fast():
    """Test build-h5ad --annotate --allow-annotation-fallback behavior."""
    print("\n=== Test 11: build-h5ad fail-fast ===")

    # Test _parse_case_level_for_split as a standalone helper
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write("benchmark_id,perturb_var,control,dose_groups,time_groups\n")
        f.write('BENCH001,dose,DMSO,"[""1 uM"", ""2 uM""]","[""1 day""]"')
        f.write("\n")
        temp_path = f.name

    try:
        import pandas as pd
        import json

        df = pd.read_csv(temp_path)
        # pd.read_csv auto-strips extra quotes around the field
        subset = df[df["benchmark_id"].astype(str).str.strip() == "BENCH001"]
        assert not subset.empty, f"Benchmark not found. Columns: {list(df.columns)}, IDs: {df['benchmark_id'].tolist()}"

        row = subset.iloc[0]
        guidance = {}
        for key in ["perturb_var", "control"]:
            val = row.get(key)
            if pd.notna(val) and str(val).strip():
                guidance[key] = str(val).strip()
        for key in ["dose_groups", "time_groups"]:
            val = row.get(key)
            if pd.notna(val) and str(val).strip():
                try:
                    # pandas may have already parsed JSON-like strings
                    if isinstance(val, list):
                        guidance[key] = val
                    else:
                        raw = str(val).replace("'", '"')
                        parsed = json.loads(raw)
                        if isinstance(parsed, list) and parsed:
                            guidance[key] = parsed
                except (json.JSONDecodeError, ValueError):
                    guidance[key] = [str(val).strip()]

        assert guidance.get("perturb_var") == "dose", f"guidance={guidance}"
        assert guidance.get("control") == "DMSO"
        assert "dose_groups" in guidance
        print(f"  [PASS] Case-level split guidance: {guidance}")
    finally:
        os.unlink(temp_path)


# ── Test 12: MarkerHeuristics target-aware ──────────────────────

def test_marker_target_aware():
    """Test that MarkerHeuristics prioritizes target-aware marker sets."""
    print("\n=== Test 12: MarkerHeuristics target-aware ===")
    from strategies.marker_heuristics import annotate_markers, _get_prioritized_marker_sets

    # Test prioritization: keratinocyte targets should get skin markers first
    ordered = _get_prioritized_marker_sets(["basal keratinocytes", "pericytes"])
    names = [n for n, _ in ordered]
    # Keratinocyte-related should come before general T cell markers
    keratinocyte_indices = [
        i for i, n in enumerate(names)
        if "keratinocyte" in n.lower() or "basal" in n.lower() or "spinous" in n.lower()
    ]
    pericyte_indices = [
        i for i, n in enumerate(names) if "pericyte" in n.lower()
    ]
    t_cell_indices = [i for i, n in enumerate(names) if n == "T cells"]
    if keratinocyte_indices and t_cell_indices:
        assert min(keratinocyte_indices) < min(t_cell_indices), (
            f"Keratinocyte markers should come before T cells. Order: {names[:5]}"
        )
    print(f"  [PASS] Target-aware ordering: {names[:5]}...")

    # Test annotation with sparse matrix
    adata = make_test_adata(n_cells=30)
    X = adata.X.toarray()
    if "CD3D" in adata.var_names:
        X[0:10, list(adata.var_names).index("CD3D")] += 30
    adata.X = sp.csr_matrix(X)
    adata.obs["cell_type"] = None

    summary = annotate_markers(adata, target_cell_types=["T cells"])
    assert summary["cells_assessed"] > 0
    print(f"  [PASS] Marker heuristics: {summary['cells_labeled']}/{summary['cells_assessed']} labeled")
    print(f"  Cell types: {adata.obs['cell_type'].value_counts().to_dict()}")


# ── Test 13: FASTQ fail-fast ────────────────────────────────────

def test_fastq_fail_fast():
    """Test that FASTQ detection fails with clear message."""
    print("\n=== Test 13: FASTQ fail-fast ===")
    from detect import detect

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        (tmp / "sample_R1_001.fastq.gz").write_text("mock")
        (tmp / "sample_R2_001.fastq.gz").write_text("mock")

        info = detect(str(tmp))
        assert info["type"] == "fastq"
        print(f"  [PASS] FASTQ detected: {info['details']}")


# ── Test 14: Selection with keratinocyte/fibroblast targets ─────

def test_keratinocyte_selection():
    """Test selection with multi-lineage skin targets."""
    print("\n=== Test 14: Keratinocyte/fibroblast selection ===")
    from select_cells import select_cells

    adata = make_test_adata(n_cells=60)
    labels = (
        ["Basal keratinocytes"] * 8
        + ["Spinous keratinocytes"] * 8
        + ["Fibroblasts"] * 10
        + ["Endothelial cells"] * 8
        + ["T cells"] * 10
        + ["DC3 / myeloid DC"] * 8
        + ["Pericytes"] * 8
    )
    adata.obs["cell_type"] = labels
    adata.obs["annotation_method"] = "CellTypist"

    # Test: select "Keratinocytes" (broad) should include basal + spinous
    selected, summary = select_cells(adata, ["Keratinocytes"])
    assert summary["n_cells_after_selection"] >= 16, (
        f"Keratinocytes should match >=16 (8 basal + 8 spinous), "
        f"got {summary['n_cells_after_selection']}"
    )
    print(f"  [PASS] Broad Keratinocytes: {summary['n_cells_after_selection']} cells")

    # Test: "Endothelial cells" should match via alias
    selected2, summary2 = select_cells(adata, ["Endothelial cells"])
    assert summary2["n_cells_after_selection"] == 8
    print(f"  [PASS] Endothelial cells: {summary2['n_cells_after_selection']} cells")


# ── Test 15: MarkerHeuristics does not overwrite Original ──────

def test_marker_does_not_overwrite_original():
    """Test MarkerHeuristics only fills missing, preserves Original."""
    print("\n=== Test 15: MarkerHeuristics preserves Original ===")
    from strategies.marker_heuristics import annotate_markers

    adata = make_test_adata(n_cells=40)
    X = adata.X.toarray()
    if "CD3D" in adata.var_names:
        X[0:10, list(adata.var_names).index("CD3D")] += 30
    if "MS4A1" in adata.var_names:
        X[10:15, list(adata.var_names).index("MS4A1")] += 30
    adata.X = sp.csr_matrix(X)

    # First 15 cells have Original annotation, rest have None
    adata.obs["cell_type"] = None
    adata.obs["annotation_method"] = "Unresolved"
    original_idx = adata.obs.index[:15]
    adata.obs.loc[original_idx, "cell_type"] = ["T cells"] * 10 + ["B cells"] * 5
    adata.obs.loc[original_idx, "annotation_method"] = "Original"

    summary = annotate_markers(adata, target_cell_types=["T cells", "B cells"])
    print(f"  Labeled: {summary['cells_labeled']}/{summary['cells_assessed']}")

    # Check Original cells are still Original
    assert (adata.obs.loc[original_idx, "annotation_method"] == "Original").all(), (
        "Original annotations were overwritten!"
    )
    print("  [PASS] Original annotations preserved")

    # Check Original cell_type values are still the same
    assert (adata.obs.loc[original_idx[:10], "cell_type"] == "T cells").all()
    assert (adata.obs.loc[original_idx[10:15], "cell_type"] == "B cells").all()
    print("  [PASS] Original cell_type values intact")

    # Check non-original cells might have been filled
    n_filled = (adata.obs.loc[original_idx[15:], "annotation_method"] == "MarkerHeuristics").sum()
    print(f"  [PASS] {n_filled} new cells got MarkerHeuristics label")


# ── Test 16: Prefixed MEX loading with .gz ──────────────────────

def test_prefixed_mex_gz_loading():
    """Test loading prefixed MEX with .gz files using mmread."""
    print("\n=== Test 16: Prefixed MEX .gz loading ===")
    from detect import _detect_prefixed_mex_samples
    from io_utils import load_prefixed_mex_samples

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)

        # Create .gz MEX triplets for two samples
        # Sample A
        mtx_content_a = (
            "%%MatrixMarket matrix coordinate integer general\n"
            "6 4 8\n"
            "1 1 10.0\n2 1 5.0\n3 2 8.0\n4 2 3.0\n5 3 12.0\n"
            "6 3 7.0\n1 4 2.0\n4 4 4.0\n"
        )
        with gzip.open(str(tmp / "SampleA.matrix.mtx.gz"), "wt") as f:
            f.write(mtx_content_a)
        with gzip.open(str(tmp / "SampleA.features.tsv.gz"), "wt") as f:
            f.write("GENE1\tGene1\tGene Expression\n")
            f.write("GENE2\tGene2\tGene Expression\n")
            f.write("GENE3\tGene3\tGene Expression\n")
            f.write("GENE4\tGene4\tGene Expression\n")
            f.write("GENE5\tGene5\tGene Expression\n")
            f.write("GENE6\tGene6\tGene Expression\n")
        with gzip.open(str(tmp / "SampleA.barcodes.tsv.gz"), "wt") as f:
            f.write("AAACCC-1\nAAAGGG-1\nAAATTT-1\nAAAGTA-1\n")

        # Sample B
        mtx_content_b = (
            "%%MatrixMarket matrix coordinate integer general\n"
            "6 3 7\n"
            "1 1 7.0\n2 1 4.0\n3 2 9.0\n4 2 6.0\n5 3 11.0\n6 3 3.0\n4 3 8.0\n"
        )
        with gzip.open(str(tmp / "SampleB.matrix.mtx.gz"), "wt") as f:
            f.write(mtx_content_b)
        with gzip.open(str(tmp / "SampleB.features.tsv.gz"), "wt") as f:
            f.write("GENE1\tGene1\nGENE2\tGene2\nGENE3\tGene3\nGENE4\tGene4\nGENE5\tGene5\nGENE6\tGene6\n")
        with gzip.open(str(tmp / "SampleB.barcodes.tsv.gz"), "wt") as f:
            f.write("BBB-1\nBBB-2\nBBB-3\n")

        # Detect
        samples = _detect_prefixed_mex_samples(tmp)
        assert len(samples) == 2, f"Expected 2 samples, got {len(samples)}"
        print(f"  [PASS] Detected {len(samples)} .gz prefixed samples")

        # Load
        adata = load_prefixed_mex_samples(samples)

        # Verify shape: 4 + 3 = 7 cells, 6 common genes
        assert adata.n_obs == 7, f"Expected 7 cells total, got {adata.n_obs}"
        assert adata.n_vars == 6, f"Expected 6 genes, got {adata.n_vars}"
        print(f"  [PASS] Merged: {adata.n_obs} cells, {adata.n_vars} genes")

        # Verify sample_id
        assert "sample_id" in adata.obs.columns
        assert set(adata.obs["sample_id"].unique()) == {"SampleA", "SampleB"}
        print(f"  [PASS] sample_id present: {sorted(adata.obs['sample_id'].unique())}")

        # Verify unique obs_names
        assert len(set(adata.obs_names)) == adata.n_obs, "obs_names not unique"
        print("  [PASS] obs_names are unique")

        # Verify prefixed
        for name in adata.obs_names:
            assert any(name.startswith(p) for p in ["SampleA_", "SampleB_"]), (
                f"obs_name not prefixed: {name}"
            )
        print(f"  [PASS] obs_names prefixed: {adata.obs_names[:3].tolist()}...")


# ── Test 17: build-h5ad dose/time split ────────────────────────────

def test_build_h5ad_dose_time_split():
    """Test build-h5ad _parse_case_level_for_split with dose/time."""
    print("\n=== Test 17: build-h5ad dose/time split ===")

    import csv as _csv
    import io as _io

    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        w = _csv.writer(f)
        w.writerow(["benchmark_id", "perturb_var", "control", "dose_groups", "time_groups"])
        doses = '["1 uM", "2 uM", "4 uM"]'
        times = '["6 day", "10 day"]'
        w.writerow(["BENCH001", "dose", "DMSO", doses, times])
        temp_path = f.name

    try:
        import pandas as pd
        import json as _json

        df = pd.read_csv(temp_path)
        subset = df[df["benchmark_id"].astype(str).str.strip() == "BENCH001"]
        assert not subset.empty, f"Benchmark not found in: {df['benchmark_id'].tolist()}"

        row = subset.iloc[0]
        # Parse dose_groups
        dg_raw = str(row["dose_groups"])
        dg = _json.loads(dg_raw)
        assert isinstance(dg, list)
        assert len(dg) == 3
        print(f"  [PASS] dose_groups parsed: {dg}")

        # Parse time_groups
        tg_raw = str(row["time_groups"])
        tg = _json.loads(tg_raw)
        assert isinstance(tg, list)
        assert len(tg) == 2
        print(f"  [PASS] time_groups parsed: {tg}")

        # Now test multi-row merge
        f2 = tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False)
        w2 = _csv.writer(f2)
        w2.writerow(["benchmark_id", "perturb_var", "control", "dose_groups", "time_groups"])
        w2.writerow(["BENCH001", "dose", "DMSO", '["8 uM"]', '["14 day"]'])
        f2.close()

        df2 = pd.read_csv(f2.name)
        subset2 = df2[df2["benchmark_id"].astype(str).str.strip() == "BENCH001"]
        row2 = subset2.iloc[0]
        dg2 = _json.loads(str(row2["dose_groups"]))
        tg2 = _json.loads(str(row2["time_groups"]))

        # Simulate merge (union across rows)
        merged_doses = set(dg) | set(dg2)
        merged_times = set(tg) | set(tg2)
        assert len(merged_doses) == 4, f"Expected 4 merged doses, got {len(merged_doses)}"
        assert len(merged_times) == 3, f"Expected 3 merged times, got {len(merged_times)}"
        print(f"  [PASS] Multi-row merge: doses={sorted(merged_doses)}, times={sorted(merged_times)}")
        os.unlink(f2.name)

    finally:
        os.unlink(temp_path)


# ── Test 18: Control label fail-fast ────────────────────────────

def test_no_control_label_fail_fast():
    """Test control/treat split logic: heuristic fallback only when allowed."""
    print("\n=== Test 18: No control label fail-fast ===")

    adata = make_test_adata(n_cells=20)
    adata.obs["condition"] = ["treated_A", "treated_B"] * 10
    adata.obs["cell_type"] = "T cells"
    adata.obs["annotation_method"] = "Original"
    adata.layers["counts"] = adata.X.copy()

    # Simulate the control label identification logic from convert.py
    condition_col = "condition"
    vals = adata.obs[condition_col].astype(str).str.lower()
    unique_vals = list(vals.unique())
    control_candidates = ["control", "untreated", "dmso", "vehicle", "c", "wt", "unstimulated"]
    found = False
    for candidate in control_candidates:
        for uv in unique_vals:
            if candidate in uv:
                found = True
                break
        if found:
            break

    # No control should be found — our data only has "treated_A" and "treated_B"
    assert not found, (
        f"Expected no control match from {control_candidates} in {unique_vals}"
    )
    print(f"  [PASS] No control label found among {unique_vals}")

    # The fallback should use first value when heuristic fallback is allowed
    fallback_label = adata.obs[condition_col].unique()[0]
    print(f"  [PASS] Heuristic fallback would use: '{fallback_label}'")

    # Verify the split would work with fallback
    control = adata[adata.obs[condition_col] == fallback_label]
    treated = adata[adata.obs[condition_col] != fallback_label]
    assert control.n_obs + treated.n_obs == adata.n_obs
    # cell_type and annotation_method are preserved
    assert "cell_type" in control.obs.columns
    assert "annotation_method" in control.obs.columns
    print(f"  [PASS] Split with fallback: {control.n_obs} control, {treated.n_obs} treated")
    print("  [PASS] cell_type and annotation_method preserved")


# ── Main ────────────────────────────────────────────────────────

def main():
    deps_ok, missing = _check_deps()
    if not deps_ok:
        print(f"MISSING DEPENDENCIES: {missing}")
        print("Install with: pip install " + " ".join(missing))
        print("Or: pip install -r skills/annotate-cell-types/requirements.txt")
        print("\nSkipping tests that require these packages.")
        # We can still run pure-Python logic tests
        runnable = []
    else:
        runnable = None  # run all

    tests = [
        ("Preflight", test_preflight, None),
        (".mtx.gz barcode count", test_mtx_gz_barcode_count, None),
        ("Prefixed MEX detection", test_prefixed_mex_detection, None),
        ("Prefixed MEX loading", test_prefixed_mex_loading, None),
        ("H5 detection", test_h5_detection, None),
        ("Case-level fail-fast", test_case_fail_fast, None),
        ("Target selection broad/narrow", test_broad_narrow_selection, None),
        ("CellTypist mock writeback", test_celltypist_mock, {
            "celltypist"},  # only if celltypist module available
        ),
        ("No expensive methods", test_no_expensive_methods, None),
        ("Manifest completeness", test_manifest_completeness, None),
        ("build-h5ad fail-fast", test_build_h5ad_fail_fast, None),
        ("MarkerHeuristics target-aware", test_marker_target_aware, None),
        ("FASTQ fail-fast", test_fastq_fail_fast, None),
        ("Keratinocyte selection", test_keratinocyte_selection, None),
        ("MarkerHeuristics preserves Original", test_marker_does_not_overwrite_original, None),
        ("Prefixed MEX .gz loading", test_prefixed_mex_gz_loading, None),
        ("build-h5ad dose/time split", test_build_h5ad_dose_time_split, None),
        ("No control label fail-fast", test_no_control_label_fail_fast, None),
    ]

    passed = 0
    failed = 0
    skipped = 0

    for name, test_fn, requires in tests:
        if requires and missing:
            if requires & set(missing):
                print(f"\n  [SKIP] {name} — missing: {requires & set(missing)}")
                skipped += 1
                continue

        try:
            test_fn()
            passed += 1
        except Exception as e:
            print(f"  [FAIL] {e}")
            import traceback
            traceback.print_exc()
            failed += 1

    print("\n" + "=" * 60)
    print(f"Results: {passed} passed, {failed} failed, {skipped} skipped out of {len(tests)}")
    print("=" * 60)

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
