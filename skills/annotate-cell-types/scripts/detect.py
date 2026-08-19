"""Detect input data format from a directory or file.

Detects:
  - FASTQ (.fastq.gz, R1/R2)
  - Cell Ranger output (outs/filtered_feature_bc_matrix or .h5)
  - 10x MEX filtered / raw (direct or prefixed multi-sample)
  - h5ad (.h5ad)
  - CSV/TSV/TXT (expression tables)
"""

import gzip
import logging
import re
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)

INPUT_TYPES = [
    "fastq",
    "cellranger_output",
    "10x_mtx_filtered",
    "10x_mtx_raw",
    "10x_h5_filtered",
    "10x_h5_raw",
    "h5ad",
    "csv",
    "unknown",
]

# Regex to detect prefixed MEX: catches <prefix>.matrix.mtx, <prefix>.matrix.mtx.gz,
# <prefix>_matrix.mtx, and prefixed barcodes/features.
_PREFIXED_MTX_RE = re.compile(
    r"^(.+?)[\._\-]matrix\.mtx(?:\.gz)?$", re.IGNORECASE
)
# Also catch *_counts.mtx, *_RNA_counts.mtx etc. with matching features/barcodes
_PREFIXED_COUNTS_RE = re.compile(
    r"^(.+?)[\._\-]?(?:_RNA)?_counts\.mtx(?:\.gz)?$", re.IGNORECASE
)


def _open_maybe_gz(path: Path, mode: str = "rt"):
    """Open a file, handling gzip transparently."""
    path_str = str(path)
    if path_str.endswith(".gz"):
        return gzip.open(path_str, mode)
    else:
        return open(path_str, mode)


def _count_mtx_barcodes(mtx_path: Path) -> int:
    """Read MTX header to get dimensions; n_barcodes = matrix columns.

    Handles both plain .mtx and .mtx.gz.
    """
    try:
        with _open_maybe_gz(mtx_path, "rt") as f:
            for line in f:
                line = line.strip()
                if line.startswith("%"):
                    continue
                parts = line.split()
                if len(parts) == 3:
                    return int(parts[1])
                break
    except Exception as e:
        logger.debug("Failed to read MTX header from %s: %s", mtx_path, e)
    return 0


def _mex_triplet_exists(directory: Path, prefix: str = "") -> tuple:
    """Check if a MEX triplet (matrix, features/genes, barcodes) exists.

    Args:
        directory: where to look
        prefix: optional filename prefix (e.g. 'GSM123_sampleA')

    Returns:
        (matrix_path, features_path, barcodes_path) or (None, None, None)
    """
    if prefix:
        # Prefixed: <prefix>.matrix.mtx(.gz), <prefix>.features.tsv(.gz), etc.
        suffix_variants = ["", ".gz"]
        separators = [".", "_", "-"]

        for sep in separators:
            for sv in suffix_variants:
                mtx = directory / f"{prefix}{sep}matrix.mtx{sv}"
                if mtx.exists():
                    for fsv in suffix_variants:
                        for fsep in separators:
                            feats = directory / f"{prefix}{fsep}features.tsv{fsv}"
                            if not feats.exists():
                                feats = directory / f"{prefix}{fsep}genes.tsv{fsv}"
                            barcs = directory / f"{prefix}{fsep}barcodes.tsv{fsv}"
                            if feats.exists() and barcs.exists():
                                return (mtx, feats, barcs)
                    # If mtx found but features/barcodes not with same sep, try all pairs
                    for fsv in suffix_variants:
                        for bsv in suffix_variants:
                            feats = None
                            barcs = None
                            for fsep in separators:
                                fc = directory / f"{prefix}{fsep}features.tsv{fsv}"
                                if fc.exists():
                                    feats = fc
                                    break
                                gc = directory / f"{prefix}{fsep}genes.tsv{fsv}"
                                if gc.exists():
                                    feats = gc
                                    break
                            for bsep in separators:
                                bc = directory / f"{prefix}{bsep}barcodes.tsv{bsv}"
                                if bc.exists():
                                    barcs = bc
                                    break
                            if feats and barcs:
                                return (mtx, feats, barcs)
    else:
        # No prefix: matrix.mtx(.gz), features.tsv(.gz), barcodes.tsv(.gz)
        for sv in ["", ".gz"]:
            m = directory / f"matrix.mtx{sv}"
            if m.exists():
                for fsv in ["", ".gz"]:
                    f = directory / f"features.tsv{fsv}"
                    if not f.exists():
                        f = directory / f"genes.tsv{fsv}"
                    b = directory / f"barcodes.tsv{fsv}"
                    if f.exists() and b.exists():
                        return (m, f, b)
    return (None, None, None)


def _detect_prefixed_mex_samples(directory: Path) -> list[dict]:
    """Detect all prefixed multi-sample MEX triplets in a directory.

    Scans for *.matrix.mtx / *.matrix.mtx.gz and matches corresponding
    features/genes.tsv and barcodes.tsv files by prefix.

    Returns list of dicts, each with sample_id, matrix_path, features_path,
    barcodes_path, barcode_count, is_raw.
    """
    samples = []

    # Find all matrix.mtx AND *_counts.mtx files
    mtx_patterns = (
        list(directory.glob("*matrix.mtx")) +
        list(directory.glob("*matrix.mtx.gz")) +
        list(directory.glob("*_counts.mtx")) +
        list(directory.glob("*_counts.mtx.gz"))
    )
    # Skip subdirectory files
    mtx_patterns = [p for p in mtx_patterns if p.parent == directory]

    seen_prefixes = set()

    for mtx_path in sorted(mtx_patterns):
        # Try both matrix.mtx and _counts.mtx patterns
        m = _PREFIXED_MTX_RE.match(mtx_path.name)
        if not m:
            m = _PREFIXED_COUNTS_RE.match(mtx_path.name)
        if m:
            prefix = m.group(1)
            if prefix in seen_prefixes:
                continue
            seen_prefixes.add(prefix)

            # Find matching features and barcodes
            feats_path = None
            barcs_path = None
            base = mtx_path.name
            # Try both "matrix.mtx" and "counts.mtx" replacement
            for mtx_suffix in ["matrix.mtx", "counts.mtx"]:
                for sep in [".", "_", "-"]:
                    search_pat = f"{sep}{mtx_suffix}"
                    if search_pat not in base:
                        continue
                    candidate_feats = directory / base.replace(
                        search_pat, f"{sep}features.tsv"
                    )
                    candidate_genes = directory / base.replace(
                        search_pat, f"{sep}genes.tsv"
                    )
                    candidate_barcs = directory / base.replace(
                        search_pat, f"{sep}barcodes.tsv"
                    )
                    for ext in ["", ".gz"]:
                        fp = Path(str(candidate_feats) + ext)
                        gp = Path(str(candidate_genes) + ext)
                        bp = Path(str(candidate_barcs) + ext)
                        if fp.exists():
                            feats_path = fp
                        elif gp.exists():
                            feats_path = gp
                        if bp.exists():
                            barcs_path = bp
                        if feats_path and barcs_path:
                            break
                    if feats_path and barcs_path:
                        break
                if feats_path and barcs_path:
                    break

            if feats_path and barcs_path:
                n_barcodes = _count_mtx_barcodes(mtx_path)
                is_raw = n_barcodes > 100_000
                samples.append({
                    "sample_id": prefix,
                    "matrix_path": str(mtx_path),
                    "features_path": str(feats_path),
                    "barcodes_path": str(barcs_path),
                    "barcode_count": n_barcodes,
                    "is_raw": is_raw,
                })

    return samples


def detect(input_path: str) -> dict:
    """Detect input data format.

    Returns dict with keys:
        type: one of INPUT_TYPES
        path: resolved input path
        details: human-readable description
        barcode_count: estimate (0 if unknown)
        is_raw_mex: bool
        samples: list of sample dicts (for prefixed MEX)
    """
    path = Path(input_path)

    # Single h5ad file
    if path.is_file() and path.suffix.lower() == ".h5ad":
        import anndata as ad
        try:
            adata = ad.read_h5ad(path, backed="r")
            return {
                "type": "h5ad",
                "path": str(path),
                "details": f"h5ad file: {path.name} ({adata.n_obs} cells, {adata.n_vars} genes)",
                "barcode_count": adata.n_obs,
                "is_raw_mex": False,
            }
        except Exception:
            return {
                "type": "h5ad",
                "path": str(path),
                "details": f"h5ad file: {path.name}",
                "barcode_count": 0,
                "is_raw_mex": False,
            }

    # Single CSV/TSV/TXT file (but not .mtx or .mtx.gz)
    if path.is_file():
        suffix = path.suffix.lower()
        base_lower = path.name.lower()
        if suffix in (".csv", ".tsv", ".txt") or (
            suffix == ".gz" and not base_lower.endswith(".mtx.gz")
        ):
            return {
                "type": "csv",
                "path": str(path),
                "details": f"Expression table: {path.name}",
                "barcode_count": 0,
                "is_raw_mex": False,
            }

    # Directory: scan for known patterns
    if path.is_dir():
        return _detect_directory(path)

    return {
        "type": "unknown",
        "path": str(path),
        "details": f"Cannot detect format for: {path}",
        "barcode_count": 0,
        "is_raw_mex": False,
    }


def _detect_directory(path: Path) -> dict:
    """Detect format inside a directory."""
    # ---- FASTQ ----
    fastq_files = list(path.rglob("*.fastq.gz")) + list(path.rglob("*.fq.gz"))
    r1_files = list(path.rglob("*R1*.fastq.gz")) + list(path.rglob("*_1.fastq.gz"))
    r2_files = list(path.rglob("*R2*.fastq.gz")) + list(path.rglob("*_2.fastq.gz"))
    if r1_files and r2_files:
        return {
            "type": "fastq",
            "path": str(path),
            "details": f"FASTQ: {len(r1_files)} R1 + {len(r2_files)} R2 file pairs",
            "barcode_count": 0,
            "is_raw_mex": False,
        }
    if fastq_files:
        return {
            "type": "fastq",
            "path": str(path),
            "details": f"FASTQ: {len(fastq_files)} files",
            "barcode_count": 0,
            "is_raw_mex": False,
        }

    # ---- Cell Ranger H5 ----
    for h5_name in ["filtered_feature_bc_matrix.h5", "raw_feature_bc_matrix.h5"]:
        h5_path = path / "outs" / h5_name
        if h5_path.exists():
            is_raw = "raw" in h5_name
            return {
                "type": "10x_h5_raw" if is_raw else "10x_h5_filtered",
                "path": str(h5_path),
                "details": (
                    f"Cell Ranger {'raw' if is_raw else 'filtered'} H5: {h5_name}"
                ),
                "barcode_count": 0,
                "is_raw_mex": is_raw,
                "mex_h5": str(h5_path),
            }
        # Also check directly in the root
        h5_root = path / h5_name
        if h5_root.exists():
            is_raw = "raw" in h5_name
            return {
                "type": "10x_h5_raw" if is_raw else "10x_h5_filtered",
                "path": str(h5_root),
                "details": f"10x H5: {h5_name}",
                "barcode_count": 0,
                "is_raw_mex": is_raw,
                "mex_h5": str(h5_root),
            }

    # ---- Cell Ranger MEX output directories ----
    for mex_sub, is_raw in [
        ("filtered_feature_bc_matrix", False),
        ("raw_feature_bc_matrix", True),
    ]:
        mex_dir = path / "outs" / mex_sub
        if mex_dir.exists():
            mtx, feats, barcs = _mex_triplet_exists(mex_dir)
            if mtx:
                n_barcodes = _count_mtx_barcodes(mtx)
                return {
                    "type": "10x_mtx_raw" if is_raw else "10x_mtx_filtered",
                    "path": str(mtx),
                    "details": (
                        f"Cell Ranger {'raw' if is_raw else 'filtered'} MEX: "
                        f"~{n_barcodes} barcodes"
                    ),
                    "barcode_count": n_barcodes,
                    "is_raw_mex": is_raw,
                    "mex_dir": str(mex_dir),
                }

    # ---- Prefixed multi-sample MEX ----
    prefixed_samples = _detect_prefixed_mex_samples(path)
    if prefixed_samples:
        total_barcodes = sum(s["barcode_count"] for s in prefixed_samples)
        has_raw = any(s["is_raw"] for s in prefixed_samples)
        return {
            "type": "10x_mtx_raw" if has_raw else "10x_mtx_filtered",
            "path": str(path),
            "details": (
                f"Prefixed multi-sample MEX: {len(prefixed_samples)} samples, "
                f"~{total_barcodes} total barcodes"
            ),
            "barcode_count": total_barcodes,
            "is_raw_mex": has_raw,
            "samples": prefixed_samples,
        }

    # ---- Direct MEX triplet (no prefix) ----
    mtx, feats, barcs = _mex_triplet_exists(path)
    if mtx:
        n_barcodes = _count_mtx_barcodes(mtx)
        is_raw = n_barcodes > 100_000
        return {
            "type": "10x_mtx_raw" if is_raw else "10x_mtx_filtered",
            "path": str(mtx),
            "details": f"10x MEX: ~{n_barcodes} barcodes",
            "barcode_count": n_barcodes,
            "is_raw_mex": is_raw,
            "mex_dir": str(path),
        }

    # ---- Any .h5ad files ----
    h5ad_files = list(path.rglob("*.h5ad"))
    if h5ad_files:
        return detect(str(h5ad_files[0]))

    # ---- Gzipped CSV/TSV/TXT (common in GEO) ----
    gz_tables = [f for f in path.iterdir()
                 if not f.name.startswith(".")
                 and f.name.lower().endswith((".csv.gz", ".tsv.gz", ".txt.gz", ".csv", ".tsv", ".txt"))
                 and "barcode" not in f.name.lower()
                 and "feature" not in f.name.lower()
                 and "gene" not in f.name.lower()]
    if gz_tables:
        return {
            "type": "csv",
            "path": str(gz_tables[0]),
            "details": f"Expression table(s): {len(gz_tables)} table(s), first={gz_tables[0].name}",
            "barcode_count": 0,
            "is_raw_mex": False,
        }

    return {
        "type": "unknown",
        "path": str(path),
        "details": "Directory with unrecognized contents",
        "barcode_count": 0,
        "is_raw_mex": False,
    }
