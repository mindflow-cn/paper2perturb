"""I/O utilities: build or load AnnData from detected input format.

Handles: h5ad, 10x MEX (direct & prefixed multi-sample), 10x H5, CSV/TSV.
For raw MEX: supports per-sample cell calling to avoid merging millions of
empty barcodes before filtering.
"""

import gzip
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread

logger = logging.getLogger(__name__)


def _parse_barcode_metadata(adata: "ad.AnnData") -> None:
    """Parse metadata from prefixed/annotated barcodes.

    Detects patterns like _d{digits} (timepoints) and _{alphanum}_
    (sample/patient IDs) embedded in barcode names after a format separator
    (e.g. GSE111014_AAACCTG-scRNA-seq_CLL6_d0).
    """
    barcodes = pd.Series(adata.obs_names.astype(str))
    sample_barcodes = barcodes[:min(100, len(barcodes))]

    # Check if barcodes have extra metadata after the core sequence.
    # Common pattern: {prefix}-{description} where description has `_` separated fields.
    has_suffix = sample_barcodes.str.contains(r'-.+_').mean()
    if has_suffix < 0.5:
        return

    # Parse timepoints: _d\d+ at end of suffix
    times = barcodes.str.extract(r'_(d\d+)$', expand=False)
    if times.notna().mean() > 0.3:
        def _norm_time(t):
            if pd.isna(t):
                return "unknown"
            t = str(t).lower()
            num = t.lstrip('d')
            if num == '0':
                return 'D0 untreated'
            return f'{num} day'
        time_col = times.apply(_norm_time)
        time_col.index = adata.obs.index
        adata.obs['time'] = time_col
        logger.info("Parsed time from barcodes: %d unique values",
                   adata.obs['time'].nunique())

    # Parse patient/sample ID: _{alphanum}_ pattern before timepoint
    patients = barcodes.str.extract(r'_([A-Za-z]+\d+)_d\d+', expand=False)
    if patients.notna().mean() > 0.3:
        patient_col = patients.fillna("unknown").astype(str)
        patient_col.index = adata.obs.index
        adata.obs['patient'] = patient_col
        logger.info("Parsed patient from barcodes: %d unique values",
                   adata.obs['patient'].nunique())


def _load_one_mex_sample(s: dict, do_cell_calling: bool = False, cell_calling_kwargs: dict = None) -> "ad.AnnData":
    """Load a single MEX sample. Optionally perform cell calling immediately.

    Args:
        s: sample dict from detect (matrix_path, features_path, barcodes_path, sample_id)
        do_cell_calling: if True, filter empty barcodes per-sample
        cell_calling_kwargs: kwargs for call_cells()

    Returns:
        AnnData with obs_names = sample_id + "_" + barcode, obs['sample_id'] = sample_id
    """
    import anndata as ad

    sample_id = s["sample_id"]

    # Read matrix
    mtx = mmread(s["matrix_path"])
    n_barcodes_from_file = s.get("barcode_count", 0)
    # Heuristic orientation: compare with barcode count, then dim ratio
    if n_barcodes_from_file > 0:
        if mtx.shape[0] == n_barcodes_from_file:
            pass
        elif mtx.shape[1] == n_barcodes_from_file:
            mtx = mtx.T.tocsr()
        elif mtx.shape[0] > mtx.shape[1] * 2:
            mtx = mtx.T.tocsr()
        else:
            mtx = mtx.tocsr()
    else:
        if mtx.shape[0] > mtx.shape[1] * 3:
            mtx = mtx.T.tocsr()
        else:
            mtx = mtx.tocsr()

    # Features
    feats_path = s["features_path"]
    features = pd.read_csv(
        feats_path, sep="\t", header=None,
    )
    if features.shape[1] == 1:
        gene_symbols = features.iloc[:, 0].astype(str)
    elif features.shape[1] == 2:
        gene_symbols = features.iloc[:, 1].fillna(features.iloc[:, 0]).astype(str)
    else:
        gene_symbols = features.iloc[:, 1].fillna(features.iloc[:, 0]).astype(str)

    # Barcodes
    barcs_path = s["barcodes_path"]
    barcodes = pd.read_csv(barcs_path, sep="\t", header=None)[0].astype(str)

    # Build
    prefixed_barcodes = [f"{sample_id}_{bc}" for bc in barcodes]
    adata_sample = ad.AnnData(
        X=mtx,
        obs=pd.DataFrame({"sample_id": sample_id}, index=prefixed_barcodes),
        var=pd.DataFrame(index=gene_symbols.values),
    )
    adata_sample.var_names_make_unique()

    # Parse metadata from barcodes (e.g. {barcode}-scRNA-seq_{patient}_{time})
    _parse_barcode_metadata(adata_sample)

    # Per-sample cell calling
    n_before = adata_sample.n_obs
    if do_cell_calling and n_before > 5000:
        from strategies.cell_calling import call_cells
        kwargs = cell_calling_kwargs or {}
        adata_sample = call_cells(adata_sample, **kwargs)
        logger.info(
            "  %s: %s -> %s cells (per-sample calling)",
            sample_id, n_before, adata_sample.n_obs,
        )

    return adata_sample


def load_prefixed_mex_samples(
    samples: list[dict],
    do_cell_calling: bool = False,
    cell_calling_kwargs: dict = None,
) -> "ad.AnnData":
    """Load multiple prefixed MEX samples and merge.

    Per-sample cell calling is supported to avoid merging millions of
    raw barcodes.

    Returns merged AnnData with obs['sample_id'] and unique obs_names.
    """
    import anndata as ad

    adatas = []
    for s in samples:
        adata_sample = _load_one_mex_sample(
            s, do_cell_calling=do_cell_calling,
            cell_calling_kwargs=cell_calling_kwargs,
        )
        adatas.append(adata_sample)

    if len(adatas) == 1:
        result = adatas[0]
    elif len(adatas) > 1:
        all_genes = set()
        for a in adatas:
            all_genes.update(a.var_names)
        if len(all_genes) > max(a.n_vars for a in adatas):
            logger.info(
                "Gene sets differ across samples — union of %d genes (outer join)",
                len(all_genes),
            )
        result = ad.concat(adatas, join="outer", fill_value=0)
        result.var_names_make_unique()
        logger.info(
            "Merged %d samples: %s cells, %s genes",
            len(adatas), result.n_obs, result.n_vars,
        )
    else:
        raise ValueError("No samples to load")

    return result


def load_10x_mex(mex_dir: str) -> "ad.AnnData":
    import scanpy as sc
    adata = sc.read_10x_mtx(mex_dir, var_names="gene_symbols", cache=False)
    adata.var_names_make_unique()
    return adata


def load_10x_h5(h5_path: str) -> "ad.AnnData":
    import scanpy as sc
    adata = sc.read_10x_h5(h5_path)
    adata.var_names_make_unique()
    return adata


def load_h5ad(h5ad_path: str) -> "ad.AnnData":
    import anndata as ad
    return ad.read_h5ad(h5ad_path)


def load_csv_table(table_path: str) -> "ad.AnnData":
    import anndata as ad
    path = Path(table_path)
    name_lower = path.name.lower()
    if name_lower.endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz")):
        sep = "\t"
    elif name_lower.endswith((".csv", ".csv.gz")):
        sep = ","
    else:
        sep = "\t" if path.suffix.lower() in (".tsv", ".txt") else ","
    try:
        df = pd.read_csv(table_path, sep=sep, index_col=0)
    except Exception:
        try:
            df = pd.read_csv(table_path, sep=sep)
        except Exception as e:
            raise ValueError(f"Failed to parse {table_path}: {e}")
    X = sp.csr_matrix(df.values)
    if df.shape[0] > df.shape[1]:
        # genes × cells → transpose to cells × genes
        adata = ad.AnnData(
            X=X.T, obs=pd.DataFrame(index=df.columns.astype(str)),
            var=pd.DataFrame(index=df.index.astype(str)),
        )
    else:
        # already cells × genes or square
        adata = ad.AnnData(
            X=X, obs=pd.DataFrame(index=df.index.astype(str)),
            var=pd.DataFrame(index=df.columns.astype(str)),
        )
    adata.var_names_make_unique()
    return adata


def build_adata(
    input_info: dict,
    input_dir: str,
    do_cell_calling: bool = False,
    cell_calling_kwargs: dict = None,
) -> "ad.AnnData":
    """Build or load AnnData from detected input.

    Args:
        input_info: dict from detect.detect()
        input_dir: root input directory
        do_cell_calling: if True and input is raw MEX, call cells per-sample
        cell_calling_kwargs: kwargs for cell calling
    """
    input_type = input_info["type"]

    if input_type == "h5ad":
        logger.info("Loading h5ad: %s", input_info["path"])
        return load_h5ad(input_info["path"])

    if input_type in ("10x_h5_filtered", "10x_h5_raw"):
        h5_path = input_info.get("mex_h5") or input_info.get("path")
        logger.info("Loading 10x H5: %s", h5_path)
        adata = load_10x_h5(h5_path)
        if input_type == "10x_h5_raw" and do_cell_calling:
            from strategies.cell_calling import call_cells
            n_before = adata.n_obs
            adata = call_cells(adata, **(cell_calling_kwargs or {}))
            logger.info("H5 cell calling: %s -> %s", n_before, adata.n_obs)
        return adata

    if "samples" in input_info and input_info["samples"]:
        logger.info(
            "Loading prefixed multi-sample MEX: %d samples", len(input_info["samples"]),
        )
        return load_prefixed_mex_samples(
            input_info["samples"],
            do_cell_calling=do_cell_calling,
            cell_calling_kwargs=cell_calling_kwargs,
        )

    if input_type in ("10x_mtx_filtered", "10x_mtx_raw", "cellranger_output"):
        mex_dir = input_info.get("mex_dir")
        if mex_dir:
            logger.info("Loading 10x MEX from: %s", mex_dir)
            adata = load_10x_mex(mex_dir)
            if input_type == "10x_mtx_raw" and do_cell_calling:
                from strategies.cell_calling import call_cells
                n_before = adata.n_obs
                adata = call_cells(adata, **(cell_calling_kwargs or {}))
                logger.info("MEX cell calling: %s -> %s", n_before, adata.n_obs)
            return adata
        h5_path = input_info.get("mex_h5")
        if h5_path:
            logger.info("Loading 10x H5 (fallback): %s", h5_path)
            return load_10x_h5(h5_path)
        raise ValueError(f"No MEX dir or H5 for type={input_type}: {input_info}")

    if input_type == "csv":
        logger.info("Loading CSV table: %s", input_info["path"])
        return load_csv_table(input_info["path"])

    if input_type == "fastq":
        raise ValueError(
            "FASTQ data requires Cell Ranger preprocessing before annotation.\n"
            "Options:\n"
            "  1. Run cellranger count externally, then provide the output directory.\n"
            "  2. Re-run with --run-cellranger --cellranger-reference <ref>.\n"
            "Input directory: " + input_dir
        )

    raise ValueError(
        f"Unknown input type: {input_type}, cannot build AnnData. Info: {input_info}"
    )


def save_adata(adata: "ad.AnnData", output_path: str):
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    adata.write(output_path, compression="gzip")
    logger.info("Saved: %s (%s cells, %s genes)", output_path, adata.n_obs, adata.n_vars)


def ensure_counts_layer(adata: "ad.AnnData"):
    if "counts" not in adata.layers:
        logger.info("Storing raw expression in layers['counts']")
        adata.layers["counts"] = adata.X.copy()
