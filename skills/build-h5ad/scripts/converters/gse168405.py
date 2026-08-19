"""Converter for GSE168405: Endometrial assembloids ± dasatinib (Drop-seq).

Experimental design:
  - Undifferentiated (Day0): ExM+E2 for 4 days — proliferation baseline
  - Decidualized (Day4): MDM for 4 days — differentiation control
  - Decidualized + dasatinib (Day4+dasatinib): MDM + 250 nM dasatinib — drug treated

Data format: CSV count matrix (gene x cell) + metadata CSV (Cell, SampleName,
Patient, Treatment, Cluster).

Cell types (from paper's scRNA-seq annotation, pre-labelled in Cluster column):
  - Epithelial: EpS1 (dividing), EpS2 (E2-responsive), EpS3 (ciliated),
                EpS4 (differentiated), EpS5 (senescent)
  - Stromal: SS1 (dividing), SS2 (E2-responsive), SS3 (pre-decidual),
             SS4 (decidual), SS5 (senescent decidual)
  - Transitional: TP (ambiguous epithelial+stromal markers)
"""

import logging
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from .base import BaseConverter

logger = logging.getLogger(__name__)

# Cluster to broad cell type mapping
_CLUSTER_TO_BROAD = {
    "SS1": "stromal", "SS2": "stromal", "SS3": "stromal",
    "SS4": "stromal", "SS5": "stromal",
    "EpS1": "epithelial", "EpS2": "epithelial", "EpS3": "epithelial",
    "EpS4": "epithelial", "EpS5": "epithelial",
    "TP": "transitional",
}

# Target clusters per cell_type_standard
_TARGET_CLUSTERS = {
    "stromal cell": {"SS1", "SS2", "SS3", "SS4", "SS5"},
    "transitional cell": {"TP"},
}


def _infer_target_clusters(metadata: dict) -> set[str]:
    """Determine which Clusters to keep based on benchmark cell_type_standard."""
    ct_std = str(metadata.get("cell_type_standard", "")).lower().strip()
    ct_orig = str(metadata.get("cell_type_original", "")).lower()

    # Direct lookup
    if ct_std in _TARGET_CLUSTERS:
        return _TARGET_CLUSTERS[ct_std]

    # Fallback: parse from cell_type_original "(SS1-SS5)" or "(TP)"
    if "stromal" in ct_std or "stromal" in ct_orig:
        return {"SS1", "SS2", "SS3", "SS4", "SS5"}
    if "transitional" in ct_std or "transitional" in ct_orig:
        return {"TP"}
    if "epithelial" in ct_std or "epithelial" in ct_orig:
        return {"EpS1", "EpS2", "EpS3", "EpS4", "EpS5"}

    # Try to extract cluster names from cell_type_original parens
    m = re.search(r'\(([^)]+)\)', str(metadata.get("cell_type_original", "")))
    if m:
        clusters = {c.strip() for c in m.group(1).replace("-", ",").split(",")}
        return {c for c in clusters if c in _CLUSTER_TO_BROAD}

    # Conservative: keep all stromal + transitional (most common for drug response)
    logger.warning(
        "Could not infer target clusters from cell_type_standard='%s' / "
        "cell_type_original='%s'. Keeping all stromal + transitional.",
        ct_std, ct_orig,
    )
    return {"SS1", "SS2", "SS3", "SS4", "SS5", "TP"}


class GSE168405Converter(BaseConverter):
    """GSE168405: dasatinib-treated endometrial assembloids (Drop-seq CSV)."""

    @staticmethod
    def detect(raw_dir: Path) -> bool:
        if not raw_dir.is_dir():
            return False
        has_counts = bool(list(raw_dir.glob("*assembloid_counts*")))
        has_metadata = bool(list(raw_dir.glob("*assembloid_metadata*")))
        return has_counts and has_metadata

    def convert(self, raw_dir: Path) -> tuple[ad.AnnData, ad.AnnData]:
        # Find input files
        counts_files = sorted(raw_dir.glob("*assembloid_counts*"))
        metadata_files = sorted(raw_dir.glob("*assembloid_metadata*"))

        if not counts_files or not metadata_files:
            raise ValueError(
                f"Missing expected files in {raw_dir}. "
                f"Counts: {counts_files}, Metadata: {metadata_files}"
            )

        counts_path = counts_files[0]
        meta_path = metadata_files[0]
        logger.info("Counts: %s", counts_path.name)
        logger.info("Metadata: %s", meta_path.name)

        # Read metadata
        meta = pd.read_csv(str(meta_path))
        logger.info("Metadata: %d cells, columns=%s", len(meta), list(meta.columns))

        # Read counts matrix (gene x cell, gene names in first column)
        # This matrix is dense CSV with genes as rows, cells as columns
        logger.info("Reading counts matrix...")
        counts_df = pd.read_csv(str(counts_path), index_col=0)
        logger.info("Counts: %d genes x %d cells", counts_df.shape[0], counts_df.shape[1])

        # Align cells between metadata and counts
        cell_ids_meta = set(meta.iloc[:, 0].astype(str).str.strip())
        cell_ids_counts = set(str(c) for c in counts_df.columns)

        common = cell_ids_meta & cell_ids_counts
        logger.info(
            "Cell alignment: %d in metadata, %d in counts, %d common",
            len(cell_ids_meta), len(cell_ids_counts), len(common),
        )

        # Filter to common cells
        meta = meta[meta.iloc[:, 0].astype(str).str.strip().isin(common)].copy()
        meta["_cell_id"] = meta.iloc[:, 0].astype(str).str.strip()
        meta = meta.set_index("_cell_id")
        meta = meta.loc[sorted(common)]

        # Build expression matrix (cells x genes)
        counts_df = counts_df[sorted(common)]
        expr = counts_df.values.T  # transpose: cells x genes
        if not isinstance(expr, np.ndarray):
            expr = np.array(expr)
        expr = sp.csr_matrix(expr.astype(np.float32))
        logger.info("Expression matrix: %d cells x %d genes", expr.shape[0], expr.shape[1])

        # Determine target clusters
        target_clusters = _infer_target_clusters(self.metadata)
        benchmark_id = self.metadata.get("benchmark_id", "unknown")
        logger.info(
            "Benchmark %s: target clusters=%s (cell_type_standard='%s')",
            benchmark_id, sorted(target_clusters),
            self.metadata.get("cell_type_standard", ""),
        )

        # Filter to target cell type
        cluster_col = None
        for col in meta.columns:
            if col.strip().lower() == "cluster":
                cluster_col = col
                break
        if cluster_col is None:
            raise ValueError(f"No 'Cluster' column found in metadata. Columns: {list(meta.columns)}")

        cluster_vals = meta[cluster_col].astype(str).str.strip()
        cell_type_mask = cluster_vals.isin(target_clusters)
        n_total = len(meta)
        n_target = cell_type_mask.sum()
        logger.info(
            "Cell type filter: %d/%d cells kept (%s). Cluster distribution: %s",
            n_target, n_total, sorted(target_clusters),
            cluster_vals[cell_type_mask].value_counts().to_dict(),
        )

        if n_target == 0:
            available = sorted(cluster_vals.unique())
            raise ValueError(
                f"No cells match target clusters {sorted(target_clusters)}. "
                f"Available clusters: {available}"
            )

        meta_filtered = meta[cell_type_mask]
        expr_filtered = expr[cell_type_mask.values]

        # Build obs metadata
        obs = pd.DataFrame(index=meta_filtered.index)
        obs["sample_id"] = "GSE168405"
        obs["source_file"] = counts_path.name
        obs["treatment"] = meta_filtered["Treatment"].astype(str).str.strip()
        obs["SampleName"] = meta_filtered["SampleName"].astype(str).str.strip()
        obs["Patient"] = meta_filtered["Patient"].astype(str).str.strip()
        obs["Cluster"] = cluster_vals[cell_type_mask].values
        obs["cell_type"] = obs["Cluster"]

        # Build AnnData
        gene_names = list(counts_df.index.astype(str))
        adata = ad.AnnData(
            X=expr_filtered,
            obs=obs,
            var=pd.DataFrame(index=gene_names),
        )
        adata.var_names = gene_names

        logger.info("Full AnnData: %d cells, %d genes", adata.n_obs, adata.n_vars)

        # Split by treatment condition
        # Control: Day4 (decidualized without dasatinib)
        # Treated: Day4+dasatinib (decidualized with dasatinib)
        treatment_vals = adata.obs["treatment"]
        control_mask = treatment_vals == "Day4"
        treated_mask = treatment_vals == "Day4+dasatinib"

        n_control = control_mask.sum()
        n_treated = treated_mask.sum()

        logger.info(
            "Split: control(Day4)=%d cells, treated(Day4+dasatinib)=%d cells. "
            "Treatment distribution: %s",
            n_control, n_treated,
            treatment_vals.value_counts().to_dict(),
        )

        if n_control == 0:
            available = sorted(treatment_vals.unique())
            raise ValueError(
                f"No Day4 (control) cells found. Available treatments: {available}"
            )
        if n_treated == 0:
            available = sorted(treatment_vals.unique())
            raise ValueError(
                f"No Day4+dasatinib (treated) cells found. Available treatments: {available}"
            )

        control = adata[control_mask].copy()
        treated = adata[treated_mask].copy()

        logger.info(
            "Control: %d cells (%s clusters), Treated: %d cells (%s clusters)",
            n_control, control.obs["Cluster"].value_counts().to_dict(),
            n_treated, treated.obs["Cluster"].value_counts().to_dict(),
        )

        # Normalize
        self.normalize(control)
        self.normalize(treated)

        return control, treated
