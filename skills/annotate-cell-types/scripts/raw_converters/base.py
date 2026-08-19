"""Base raw converter for turning raw data into annotation-ready h5ad.

Subclasses implement detect() and convert().
Shared utilities: AnnData building, orientation detection, gene dedup,
metadata joining, conversion manifest.
"""

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

logger = logging.getLogger(__name__)


class BaseRawConverter:
    """Convert raw expression data to annotation-ready AnnData.

    Subclass must implement:
      - detect(raw_dir) -> bool  (staticmethod)
      - convert(raw_dir, dataset_id=None) -> AnnData

    Output AnnData must have:
      - X: cells x genes (sparse CSR preferred)
      - layers["counts"]: raw counts when available
      - obs["sample_id"], obs["source_file"]
      - obs["condition"], treatment, drug, dose, time/time_label (when parsable)
      - obs["cell_type"] / obs["annotation_method"] (if original annotation exists)
      - var_names: gene symbols, unique
      - uns["raw_conversion"]: conversion manifest
    """

    def __init__(self):
        self.warnings = []

    @staticmethod
    def detect(raw_dir: Path) -> bool:
        raise NotImplementedError

    def convert(self, raw_dir: Path, dataset_id: str = None) -> "ad.AnnData":
        raise NotImplementedError

    # ── Shared utilities ───────────────────────────────────────

    @staticmethod
    def orient_cells_genes(mtx: sp.spmatrix, n_cells_expected: int = 0) -> sp.spmatrix:
        """Ensure matrix is cells × genes. Returns CSR."""
        if n_cells_expected > 0:
            if mtx.shape[1] == n_cells_expected:
                return mtx.T.tocsr()
            if mtx.shape[0] == n_cells_expected:
                return mtx.tocsr()
        if mtx.shape[0] > mtx.shape[1] * 2:
            return mtx.T.tocsr()
        return mtx.tocsr()

    @staticmethod
    def dedup_genes(gene_names: list) -> list:
        """Make gene names unique by appending suffix."""
        seen = {}
        result = []
        for g in gene_names:
            g = str(g).strip()
            if g in seen:
                seen[g] += 1
                result.append(f"{g}_{seen[g]}")
            else:
                seen[g] = 0
                result.append(g)
        return result

    @staticmethod
    def read_table(path: str, **kwargs) -> pd.DataFrame:
        """Read CSV/TSV with auto-detection."""
        p = Path(path)
        name_lower = p.name.lower()
        if name_lower.endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz")):
            sep = "\t"
        elif name_lower.endswith((".csv", ".csv.gz")):
            sep = ","
        else:
            sep = "\t" if any(p.name.endswith(x) for x in (".tsv", ".txt")) else ","
        try:
            return pd.read_csv(path, sep=sep, **kwargs)
        except Exception:
            return pd.read_csv(path, **kwargs)

    def validate_adata(self, adata: "ad.AnnData"):
        """Lightweight validation of converter output. Appends warnings."""
        if adata.n_obs == 0 or adata.n_vars == 0:
            self.warnings.append("Empty AnnData (0 cells or 0 genes)")
        if "counts" not in adata.layers:
            self.warnings.append("Missing layers['counts']")
        required_obs = ["sample_id", "source_file", "dataset_id"]
        for col in required_obs:
            if col not in adata.obs.columns:
                self.warnings.append(f"Missing required obs column: {col}")
        # Check metadata coverage
        for col in ["genotype", "treatment", "condition", "time_label", "drug"]:
            if col in adata.obs.columns:
                missing_rate = adata.obs[col].isna().mean()
                if missing_rate > 0.5:
                    self.warnings.append(
                        f"obs['{col}'] has {missing_rate:.1%} missing values"
                    )
        if hasattr(adata, "uns"):
            if "raw_conversion" in adata.uns:
                adata.uns["raw_conversion"]["warnings"] = self.warnings

    def build_manifest(self, converter_name: str, input_dir: str,
                       n_cells: int, n_genes: int, obs_cols: list,
                       files_loaded: list, metadata_parsed: dict = None) -> dict:
        return {
            "converter": converter_name,
            "input_dir": str(input_dir),
            "n_cells": int(n_cells),
            "n_genes": int(n_genes),
            "obs_columns": list(obs_cols),
            "files_loaded": files_loaded,
            "metadata_parsed": metadata_parsed or {},
            "warnings": self.warnings,
        }

    def try_join_metadata(self, adata: "ad.AnnData", meta_path: str,
                          join_col: str = "barcode") -> bool:
        """Try to join a metadata CSV to adata.obs by barcode.

        Handles prefixed barcodes (e.g. 'GSE296117_RA_1_...' matching 'RA_1_...').
        Returns True if successful.
        """
        if not Path(meta_path).exists():
            return False
        try:
            meta = pd.read_csv(meta_path, index_col=0)
        except Exception:
            try:
                meta = pd.read_csv(meta_path)
            except Exception as e:
                self.warnings.append(f"Failed to read {meta_path}: {e}")
                return False

        # Try direct join on obs_names
        common = set(adata.obs_names) & set(meta.index.astype(str))
        if len(common) > len(adata.obs_names) * 0.1:
            for col in meta.columns:
                if col not in adata.obs.columns:
                    adata.obs[col] = "unknown"
                mapped = meta[col].reindex(adata.obs_names)
                adata.obs[col] = mapped.fillna(adata.obs[col]).values
            logger.info("Joined metadata from %s (%d columns, %d matches)",
                       Path(meta_path).name, len(meta.columns), len(common))
            return True

        # Try stripping prefix from obs_names
        if adata.obs_names.str.contains("_").any():
            short_names = adata.obs_names.str.replace(r"^[^_]+_", "", regex=True)
            common2 = set(short_names) & set(meta.index.astype(str))
            if len(common2) > len(adata.obs_names) * 0.1:
                lookup = {s: o for s, o in zip(short_names, adata.obs_names)}
                for col in meta.columns:
                    if col not in adata.obs.columns:
                        adata.obs[col] = "unknown"
                    for short_n, orig_n in lookup.items():
                        if short_n in meta.index:
                            adata.obs.loc[orig_n, col] = str(meta.loc[short_n, col])
                logger.info("Joined metadata from %s via prefix-stripped barcodes (%d matches)",
                           Path(meta_path).name, len(common2))
                return True

        self.warnings.append(f"No barcode match for {meta_path}")
        return False

    def _parse_barcode_metadata(self, adata: "ad.AnnData") -> None:
        """Parse metadata from prefixed barcodes (e.g. AAA-scRNA-seq_CLL6_d0).

        Detects patterns like _d{digits} (timepoints) and _{alphanum}_
        (sample/patient IDs) embedded in barcode names after a format separator.
        """
        import re

        barcodes = pd.Series(adata.obs_names.astype(str))
        sample_barcodes = barcodes[:min(100, len(barcodes))]

        # Check if barcodes have extra metadata after the core sequence.
        # Common pattern: {seq}-{description} where description has `_` separated fields.
        has_suffix = sample_barcodes.str.contains(r'-.+_').mean()
        if has_suffix < 0.5:
            return

        suffixes = sample_barcodes.str.extract(r'-(.+)$', expand=False)
        if suffixes.isna().all():
            return

        # Parse timepoints: _d\d+ at end of suffix
        times = barcodes.str.extract(r'_(d\d+)$', expand=False)
        if times.notna().mean() > 0.3:
            # Normalise: d0 -> D0 untreated, d30 -> 30 day, etc.
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
