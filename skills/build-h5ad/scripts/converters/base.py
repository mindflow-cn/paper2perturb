"""Base converter class with shared normalization logic."""

import re
from pathlib import Path
from typing import Optional

import anndata as ad
import scanpy as sc
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.sparse import csr_matrix


def _normalize_label(text: str) -> str:
    """Normalize a label for case-insensitive comparison.

    Collapses whitespace, underscores, and hyphens so that
    "1.0_uM", "1.0-uM", and "1.0 uM" all match.
    """
    return re.sub(r"\s+", " ", str(text).lower().replace("_", " ").replace("-", " ")).strip()


class BaseConverter:
    """Abstract base class for dataset-specific converters.

    Subclasses must implement:
    - detect(raw_dir) -> bool (staticmethod)
    - convert(raw_dir) -> (control_adata, treated_adata)
    """

    def __init__(self, metadata: dict):
        self.metadata = metadata
        self.benchmark_id = metadata.get("benchmark_id", "unknown")
        self.control_type = metadata.get("control_type", "untreated")
        self.perturbation_name = metadata.get("perturbation_name", "")

    @staticmethod
    def detect(raw_dir: Path) -> bool:
        raise NotImplementedError

    def convert(self, raw_dir: Path) -> tuple[ad.AnnData, ad.AnnData]:
        raise NotImplementedError

    # ---- Shared normalization ----

    @staticmethod
    def normalize(adata: ad.AnnData) -> ad.AnnData:
        """CPM normalization (target_sum=1e4) + log1p.

        Stores raw counts in adata.layers['counts'] before normalizing.
        """
        adata.layers["counts"] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        # Use percent_top values that don't exceed n_vars (default [50,100,200,500]
        # fails for small gene panels like Tapestri with ~26 genes).
        max_top = min(adata.n_vars - 1, 50)
        pct_top = [50] if max_top >= 50 else ([max_top] if max_top > 0 else None)
        sc.pp.calculate_qc_metrics(adata, inplace=True, percent_top=pct_top)
        return adata

    @staticmethod
    def build_adata(
        expr_matrix: csr_matrix,
        var_names: list[str],
        obs_names: list[str],
        obs_metadata: dict[str, list] = None,
    ) -> ad.AnnData:
        """Build an AnnData object from a sparse expression matrix.

        expr_matrix: shape (n_cells, n_genes)
        """
        adata = ad.AnnData(X=expr_matrix)
        adata.obs_names = obs_names
        adata.var_names = var_names
        if obs_metadata:
            for key, values in obs_metadata.items():
                adata.obs[key] = values
        return adata

    @staticmethod
    def split_by_condition(
        adata: ad.AnnData, condition_col: str, control_value, treated_value
    ) -> tuple[ad.AnnData, ad.AnnData]:
        """Split AnnData into control and treated based on a condition column."""
        control = adata[adata.obs[condition_col] == control_value].copy()
        treated = adata[adata.obs[condition_col] == treated_value].copy()
        return control, treated

    @staticmethod
    def _to_seurat_format(
        adata: ad.AnnData,
        orig_ident_value: str,
        treated_values: set = None,
        perturb_var: str = None,
    ) -> ad.AnnData:
        """Format AnnData obs/var to Seurat-style columns matching data/example.

        Target obs columns: orig.ident, nCount_RNA, nFeature_RNA, percent.mt,
        cluster_id_res0.2
        Target var columns: adds gene_ids

        When treated_values is provided (from Case-level dose_groups / time_groups),
        orig.ident on ground_truth is STRICTLY validated: every unique value must
        be a member of the allowed set.  This enforces the contract that
        validate.py --sample-col orig.ident can map each group back to a
        Case-level-defined treatment group.
        """
        # Compute QC metrics from raw counts
        if "counts" in adata.layers:
            raw = adata.layers["counts"]
            raw_sum = np.array(raw.sum(axis=1)).ravel()
            adata.obs["nCount_RNA"] = raw_sum
            adata.obs["nFeature_RNA"] = np.array((raw > 0).sum(axis=1)).ravel()

            mt_mask = [str(g).startswith("MT-") for g in adata.var_names]
            if any(mt_mask):
                mt_counts = np.array(raw[:, mt_mask].sum(axis=1)).ravel()
                adata.obs["percent.mt"] = np.where(
                    raw_sum > 0, mt_counts / raw_sum * 100, 0.0
                )
            else:
                adata.obs["percent.mt"] = 0.0
        else:
            adata.obs["nCount_RNA"] = np.array(adata.X.sum(axis=1)).ravel()
            adata.obs["nFeature_RNA"] = np.array((adata.X > 0).sum(axis=1)).ravel()
            adata.obs["percent.mt"] = 0.0

        group_candidates = [
            "dose", "dose_label", "time", "timepoint",
            "Treatment_group", "treatment", "condition",
            "population", "group", "Group", "drug", "perturbation",
        ]

        # Set orig.ident
        if orig_ident_value == "C":
            # Control: always "C" (matching data/example)
            adata.obs["orig.ident"] = "C"
        elif treated_values and len(treated_values) > 0:
            # ── Guided mode: find column whose values match Case-level groups ──
            treated_norm = {_normalize_label(v) for v in treated_values}
            best_col = None
            # First pass: group_candidates in priority order
            for c in group_candidates:
                if c in adata.obs.columns:
                    col_vals = adata.obs[c].dropna().astype(str)
                    col_norm = set(col_vals.apply(_normalize_label))
                    if col_norm and col_norm.issubset(treated_norm):
                        best_col = c
                        break
            # Second pass: any obs column
            if best_col is None:
                for c in adata.obs.columns:
                    col_vals = adata.obs[c].dropna().astype(str)
                    col_norm = set(col_vals.apply(_normalize_label))
                    if (
                        col_norm and col_norm.issubset(treated_norm)
                        and col_norm != {"nan"}
                        and col_norm != {""}
                    ):
                        best_col = c
                        break
            # Third pass: fall back to existing heuristic
            if best_col is None:
                group_col = None
                for c in group_candidates:
                    if c in adata.obs.columns and adata.obs[c].nunique() > 1:
                        group_col = c
                        break
                if group_col is None:
                    for c in group_candidates:
                        if c in adata.obs.columns and adata.obs[c].nunique() >= 1:
                            group_col = c
                            break
                if group_col:
                    adata.obs["orig.ident"] = adata.obs[group_col].astype(str)
                else:
                    adata.obs["orig.ident"] = "T"
            else:
                adata.obs["orig.ident"] = adata.obs[best_col].astype(str)

            # ── Strict validation (only when best_col was found) ──
            if best_col is not None:
                actual_vals = adata.obs["orig.ident"].dropna().astype(str)
                actual_norm = set(actual_vals.apply(_normalize_label))
                if not actual_norm.issubset(treated_norm):
                    group_type = perturb_var if perturb_var else "dose/time"
                    raise ValueError(
                        f"orig.ident values in ground_truth must be members of "
                        f"Case-level {group_type}_groups.\n"
                        f"  Expected: {sorted(treated_values)}\n"
                        f"  Actual:   {sorted(actual_vals.unique())}\n"
                        f"  Column used: {best_col or group_col or 'N/A'}\n"
                        f"  Available obs columns: {sorted(adata.obs.columns)}"
                    )
        else:
            # Ground truth: try to preserve dose/time/treatment group labels
            group_col = None
            # First: multi-value columns (dose/time gradients)
            for c in group_candidates:
                if c in adata.obs.columns and adata.obs[c].nunique() > 1:
                    group_col = c
                    break
            # Second: single-value columns (single treatments)
            if group_col is None:
                for c in group_candidates:
                    if c in adata.obs.columns and adata.obs[c].nunique() >= 1:
                        group_col = c
                        break
            if group_col:
                adata.obs["orig.ident"] = adata.obs[group_col].astype(str)
            else:
                adata.obs["orig.ident"] = "T"

        # Clustering at resolution 0.2
        n_comps = min(30, adata.n_vars - 1, max(1, adata.n_obs - 1))
        LARGE_DATA = 100_000
        if n_comps >= 2 and adata.n_obs >= 3:
            sc.pp.pca(adata, n_comps=n_comps)
            if adata.n_obs > LARGE_DATA:
                # For large datasets, use mini-batch k-means on PCA space
                # (neighbors+leiden is O(n^2) and prohibitive)
                try:
                    from sklearn.cluster import MiniBatchKMeans
                    n_pcs = min(10, adata.obsm["X_pca"].shape[1])
                    n_clusters = min(20, max(2, adata.n_obs // 5000))
                    km = MiniBatchKMeans(
                        n_clusters=n_clusters, random_state=42,
                        n_init=3, batch_size=10000,
                    )
                    adata.obs["cluster_id_res0.2"] = (
                        km.fit_predict(adata.obsm["X_pca"][:, :n_pcs])
                    ).astype(str)
                except Exception:
                    adata.obs["cluster_id_res0.2"] = "0"
            else:
                try:
                    sc.pp.neighbors(adata)
                    sc.tl.leiden(adata, resolution=0.2, key_added="cluster_id_res0.2")
                except Exception:
                    n_pcs = min(10, adata.obsm["X_pca"].shape[1])
                    n_clusters = max(2, min(8, adata.n_obs // 10))
                    Z = linkage(adata.obsm["X_pca"][:, :n_pcs], method="ward")
                    adata.obs["cluster_id_res0.2"] = (
                        fcluster(Z, t=n_clusters, criterion="maxclust") - 1
                    )
        else:
            adata.obs["cluster_id_res0.2"] = "0"

        adata.var["gene_ids"] = adata.var_names

        # Ensure essential obs columns are present; preserve all existing obs columns.
        # This is especially important for the annotation path where cell_type,
        # annotation_method, cell_type_confidence, sample_id, and metadata columns
        # must survive through to control/ground_truth h5ad.
        always_keep = {
            "orig.ident", "nCount_RNA", "nFeature_RNA",
            "percent.mt", "cluster_id_res0.2",
        }
        annotation_keep = {
            "cell_type", "annotation_method", "cell_type_confidence",
            "sample_id", "sample",
        }
        split_keep = {
            "condition", "treatment", "dose", "time",
            "dose_label", "timepoint", "Time", "Dose",
            "Treatment_group", "drug", "perturbation",
            "population", "group",
        }
        keep_obs = set(adata.obs.columns)  # preserve all
        keep_obs |= always_keep
        keep_obs |= annotation_keep
        keep_obs |= split_keep

        for col in list(adata.obs.columns):
            if col not in keep_obs:
                del adata.obs[col]

        # Clear obsm and uns
        for key in list(adata.obsm.keys()):
            del adata.obsm[key]
        for key in list(adata.uns.keys()):
            del adata.uns[key]

        return adata

    @staticmethod
    def save_results(
        control: ad.AnnData,
        treated: ad.AnnData,
        output_dir: Path,
        treated_values: set = None,
        perturb_var: str = None,
    ) -> dict:
        """Save control.h5ad and ground_truth.h5ad to output_dir.

        Formats both to Seurat-style columns before writing.

        When treated_values is provided (from Case-level dose_groups /
        time_groups), orig.ident on ground_truth is strictly validated:
        every unique value must be in the allowed set.
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        control_path = output_dir / "control.h5ad"
        treated_path = output_dir / "ground_truth.h5ad"

        BaseConverter._to_seurat_format(control, "C")
        BaseConverter._to_seurat_format(
            treated, "T",
            treated_values=treated_values,
            perturb_var=perturb_var,
        )

        control.write(control_path)
        treated.write(treated_path)
        return {
            "control_file": str(control_path),
            "control_n_cells": control.n_obs,
            "treated_file": str(treated_path),
            "treated_n_cells": treated.n_obs,
            "n_genes": control.n_vars,
        }
