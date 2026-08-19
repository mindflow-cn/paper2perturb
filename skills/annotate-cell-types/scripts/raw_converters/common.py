"""Generic raw converters for common data formats.

- TenXMEXConverter: standard 10x MEX triplet
- TenXH5Converter: 10x HDF5 / CellRanger output
- SingleCSVConverter: single TSV/CSV expression matrix
- MultiCSVConverter: multiple TSV/CSV files (one per sample)
- H5adPassthroughConverter: existing h5ad files
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread

from .base import BaseRawConverter

logger = logging.getLogger(__name__)


class TenXMEXConverter(BaseRawConverter):
    """Standard 10x MEX (matrix.mtx + features.tsv + barcodes.tsv)."""

    @staticmethod
    def detect(raw_dir: Path) -> bool:
        if not raw_dir.is_dir():
            return False
        mtx = raw_dir / "matrix.mtx"
        mtx_gz = raw_dir / "matrix.mtx.gz"
        if not (mtx.exists() or mtx_gz.exists()):
            return False
        has_feat = any((raw_dir / f).exists() for f in
                       ["features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"])
        has_barc = any((raw_dir / f).exists() for f in
                       ["barcodes.tsv", "barcodes.tsv.gz"])
        return has_feat and has_barc

    def convert(self, raw_dir: Path, dataset_id: str = None) -> "ad.AnnData":
        import anndata as ad

        mtx_file = raw_dir / "matrix.mtx.gz" if (raw_dir / "matrix.mtx.gz").exists() else raw_dir / "matrix.mtx"
        mtx = mmread(str(mtx_file))
        # Features
        for fn in ["features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"]:
            fp = raw_dir / fn
            if fp.exists():
                feats = self.read_table(str(fp), header=None)
                break
        if feats.shape[1] >= 2:
            genes = feats.iloc[:, 1].fillna(feats.iloc[:, 0]).astype(str).rename(None)
        else:
            genes = feats.iloc[:, 0].astype(str).rename(None)
        # Barcodes
        for fn in ["barcodes.tsv.gz", "barcodes.tsv"]:
            fp = raw_dir / fn
            if fp.exists():
                barcs = pd.read_csv(str(fp), sep="\t", header=None)[0].astype(str).rename(None)
                break
        # Orient
        mtx = self.orient_cells_genes(mtx, len(barcs))
        genes = list(self.dedup_genes(genes))
        sid = dataset_id or raw_dir.parent.name
        adata = ad.AnnData(
            X=mtx, obs=pd.DataFrame({"sample_id": sid, "source_file": str(raw_dir.name)},
                                    index=barcs),
            var=pd.DataFrame(index=genes),
        )
        adata.layers["counts"] = adata.X.copy()
        adata.var_names = self.dedup_genes(adata.var_names)
        adata.uns["raw_conversion"] = self.build_manifest(
            "TenXMEXConverter", str(raw_dir), adata.n_obs, adata.n_vars,
            list(adata.obs.columns), [str(mtx_file)],
        )
        # Try metadata join
        for meta_name in ["*_meta.csv", "*_metadata.csv"]:
            for mp in raw_dir.glob(meta_name):
                self.try_join_metadata(adata, str(mp))
        # Parse metadata from prefixed barcodes (e.g. {barcode}-scRNA-seq_{patient}_{time})
        self._parse_barcode_metadata(adata)
        return adata


class TenXH5Converter(BaseRawConverter):
    """10x HDF5 or CellRanger outs/ directory."""

    @staticmethod
    def detect(raw_dir: Path) -> bool:
        if not raw_dir.is_dir():
            return False
        h5_direct = list(raw_dir.glob("*.h5"))
        h5_outs = raw_dir / "outs"
        if h5_outs.exists():
            h5_direct += list(h5_outs.glob("*.h5"))
        return len(h5_direct) > 0

    def convert(self, raw_dir: Path, dataset_id: str = None) -> "ad.AnnData":
        import anndata as ad
        import scanpy as sc

        h5_files = list(raw_dir.glob("*.h5"))
        outs = raw_dir / "outs"
        if outs.exists():
            h5_files += list(outs.glob("filtered_feature_bc_matrix.h5"))
        if not h5_files:
            raise FileNotFoundError(f"No .h5 file found in {raw_dir}")

        h5_path = h5_files[0]
        adata = sc.read_10x_h5(str(h5_path))
        adata.var_names = self.dedup_genes(adata.var_names)
        adata.obs["sample_id"] = dataset_id or raw_dir.parent.name
        adata.obs["source_file"] = h5_path.name
        if "counts" not in adata.layers:
            adata.layers["counts"] = adata.X.copy()
        adata.uns["raw_conversion"] = self.build_manifest(
            "TenXH5Converter", str(raw_dir), adata.n_obs, adata.n_vars,
            list(adata.obs.columns), [str(h5_path)],
        )
        return adata


class SingleCSVConverter(BaseRawConverter):
    """Single CSV/TSV expression matrix."""

    @staticmethod
    def detect(raw_dir: Path) -> bool:
        if raw_dir.is_file():
            name_lower = raw_dir.name.lower()
            if raw_dir.suffix.lower() in (".csv", ".tsv", ".txt"):
                return True
            if name_lower.endswith((".csv.gz", ".tsv.gz", ".txt.gz")):
                return True
            return False
        if raw_dir.is_dir():
            tables = [f for f in raw_dir.iterdir()
                      if not f.name.startswith(".")
                      and (f.suffix.lower() in (".csv", ".tsv", ".txt")
                           or f.name.lower().endswith((".csv.gz", ".tsv.gz", ".txt.gz")))]
            return len(tables) == 1
        return False

    def convert(self, raw_dir: Path, dataset_id: str = None) -> "ad.AnnData":
        import anndata as ad

        if raw_dir.is_file():
            table_path = raw_dir
            source = raw_dir.name
        else:
            tables = [f for f in raw_dir.iterdir()
                      if not f.name.startswith(".")
                      and (f.suffix.lower() in (".csv", ".tsv", ".txt")
                           or f.name.lower().endswith((".csv.gz", ".tsv.gz", ".txt.gz")))]
            table_path = tables[0]
            source = table_path.name

        df = self.read_table(str(table_path), index_col=0)
        X = sp.csr_matrix(df.values)
        if df.shape[0] > df.shape[1]:
            X = X.T
            obs_names = df.columns.astype(str)
            var_names = df.index.astype(str)
        else:
            obs_names = df.index.astype(str)
            var_names = df.columns.astype(str)

        sid = dataset_id or raw_dir.parent.name
        adata = ad.AnnData(
            X=X.tocsr(),
            obs=pd.DataFrame({"sample_id": sid, "source_file": source},
                             index=obs_names),
            var=pd.DataFrame(index=list(self.dedup_genes(var_names))),
        )
        adata.layers["counts"] = adata.X.copy()
        adata.uns["raw_conversion"] = self.build_manifest(
            "SingleCSVConverter", str(raw_dir), adata.n_obs, adata.n_vars,
            list(adata.obs.columns), [str(table_path)],
        )
        return adata


class MultiCSVConverter(BaseRawConverter):
    """Multiple CSV/TSV expression matrices — one per sample."""

    @staticmethod
    def detect(raw_dir: Path) -> bool:
        if not raw_dir.is_dir():
            return False
        # Exclude metadata, barcodes, features, and single-column files
        excluded = {".idf.txt", ".sdrf.txt", "metadata", "meta",
                     "_barcodes", "_features", "_genes", "barcodes.tsv",
                     "features.tsv", "genes.tsv"}
        tables = [f for f in raw_dir.iterdir()
                  if f.suffix.lower() in (".csv", ".tsv", ".txt")
                  and not f.name.startswith(".")
                  and not any(ex in f.name.lower() for ex in excluded)]
        return len(tables) >= 2

    @staticmethod
    def _is_expression_matrix(df: pd.DataFrame) -> bool:
        """Heuristic: does a DataFrame look like an expression matrix?"""
        if df.shape[0] < 5 or df.shape[1] < 5:
            return False
        # Check if most values are numeric
        numeric_cols = sum(1 for c in df.columns
                          if pd.api.types.is_numeric_dtype(df[c]))
        if numeric_cols < df.shape[1] * 0.5:
            return False
        # If single column, it's not an expression matrix (likely barcodes/metadata)
        if df.shape[1] == 1:
            return False
        return True

    def convert(self, raw_dir: Path, dataset_id: str = None) -> "ad.AnnData":
        import anndata as ad

        excluded = {".idf.txt", ".sdrf.txt", "metadata", "meta",
                     "_barcodes", "_features", "_genes", "barcodes.tsv",
                     "features.tsv", "genes.tsv"}
        tables = sorted([f for f in raw_dir.iterdir()
                         if f.suffix.lower() in (".csv", ".tsv", ".txt")
                         and not f.name.startswith(".")
                         and not any(ex in f.name.lower() for ex in excluded)])
        if not tables:
            raise FileNotFoundError(f"No expression tables in {raw_dir}")

        # Filter out non-expression tables
        valid_tables = []
        for tp in tables:
            try:
                df = pd.read_csv(str(tp), sep="\t" if tp.suffix == ".tsv" else ",", index_col=0, nrows=10)
                if self._is_expression_matrix(df):
                    valid_tables.append(tp)
                else:
                    logger.info("  Skipping %s: does not look like expression matrix", tp.name)
            except Exception:
                continue
        tables = valid_tables
        if not tables:
            raise FileNotFoundError(f"No valid expression tables in {raw_dir}")

        adatas = []
        for tp in tables:
            sample_name = tp.stem
            df = self.read_table(str(tp), index_col=0)
            X = sp.csr_matrix(df.values)
            if df.shape[0] > df.shape[1]:
                X = X.T
                obs_names = [f"{sample_name}_{c}" for c in df.columns.astype(str)]
                var_names = df.index.astype(str)
            else:
                obs_names = [f"{sample_name}_{c}" for c in df.index.astype(str)]
                var_names = df.columns.astype(str)

            a = ad.AnnData(
                X=X.tocsr(),
                obs=pd.DataFrame({"sample_id": sample_name, "source_file": tp.name},
                                 index=obs_names),
                var=pd.DataFrame(index=list(self.dedup_genes(var_names))),
            )
            adatas.append(a)
            logger.info("  Loaded %s: %d cells, %d genes", sample_name, a.n_obs, a.n_vars)

        merged = ad.concat(adatas, join="outer", fill_value=0)
        merged.var_names = self.dedup_genes(merged.var_names)
        merged.layers["counts"] = merged.X.copy()
        merged.uns["raw_conversion"] = self.build_manifest(
            "MultiCSVConverter", str(raw_dir), merged.n_obs, merged.n_vars,
            list(merged.obs.columns), [str(t) for t in tables],
        )
        # Try metadata join
        for pat in ["*_meta.csv", "*_metadata.csv", "*.sdrf.txt", "sample_sheet.tsv"]:
            for mp in raw_dir.glob(pat):
                self.try_join_metadata(merged, str(mp))
        self.validate_adata(merged)
        logger.info("Merged %d samples: %d cells, %d genes", len(tables), merged.n_obs, merged.n_vars)
        return merged


class H5adPassthroughConverter(BaseRawConverter):
    """Existing .h5ad file — pass through."""

    @staticmethod
    def detect(raw_dir: Path) -> bool:
        if raw_dir.is_file() and raw_dir.suffix.lower() == ".h5ad":
            return True
        if raw_dir.is_dir():
            return len(list(raw_dir.glob("*.h5ad"))) > 0
        return False

    def convert(self, raw_dir: Path, dataset_id: str = None) -> "ad.AnnData":
        import anndata as ad

        if raw_dir.is_file():
            h5ad_path = raw_dir
        else:
            h5ad_path = list(raw_dir.glob("*.h5ad"))[0]

        adata = ad.read_h5ad(str(h5ad_path))
        if dataset_id:
            adata.obs["dataset_id"] = dataset_id
        if "sample_id" not in adata.obs.columns:
            adata.obs["sample_id"] = dataset_id or h5ad_path.stem
        if "source_file" not in adata.obs.columns:
            adata.obs["source_file"] = h5ad_path.name
        adata.uns["raw_conversion"] = self.build_manifest(
            "H5adPassthroughConverter", str(raw_dir), adata.n_obs, adata.n_vars,
            list(adata.obs.columns), [str(h5ad_path)],
        )
        return adata
