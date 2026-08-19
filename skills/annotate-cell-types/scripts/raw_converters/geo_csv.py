"""GEO CSV raw converter.

Handles GEO supplementary file format: counts.csv(.gz) + metadata.csv(.gz).
Expression matrix: genes × cells (first column = gene names).
Metadata: cell barcodes with condition, cell_type, patient, etc.
"""

import gzip
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .base import BaseRawConverter

logger = logging.getLogger(__name__)

# Cell type column candidates to copy from metadata
_CELL_TYPE_COLS = {
    "cell_type", "celltype", "cell.types", "cell.type",
    "annotation", "cell_annotation", "cluster_name",
}


class GeoCSVConverter(BaseRawConverter):
    """GEO CSV format: counts.csv + metadata.csv (genes x cells matrix)."""

    @staticmethod
    def detect(raw_dir: Path) -> bool:
        if not raw_dir.is_dir():
            return False

        def _is_table(p: Path) -> bool:
            ext = p.suffix.lower()
            if ext in (".csv", ".tsv", ".txt"):
                return True
            if ext == ".gz":
                stem = Path(p.stem).suffix.lower()
                return stem in (".csv", ".tsv", ".txt")
            return False

        tables = [f for f in raw_dir.iterdir() if _is_table(f) and not f.name.startswith(".")]

        has_counts = any(
            any(kw in f.name.lower() for kw in ("count", "matrix", "expression", "expr", "tpm", "fpkm", "rpkm", "processed", "dge"))
            for f in tables
        )
        has_meta = any(
            any(kw in f.name.lower() for kw in ("metadata", "meta", "annotation"))
            for f in tables
        )

        if has_counts and has_meta:
            return True

        # Multiple expression matrices (no metadata) — e.g. one TPM file per sample
        if has_counts and len(tables) >= 2:
            return True

        # Fallback: exactly 2 tables, one much smaller than the other
        sizes = [(f, f.stat().st_size) for f in tables]
        if len(sizes) == 2 and (
            sizes[0][1] > 5 * sizes[1][1] or sizes[1][1] > 5 * sizes[0][1]
        ):
            return True

        return False

    @staticmethod
    def _detect_genes_as_rows(df: pd.DataFrame, n_sample: int = 100) -> bool:
        """Detect whether the DataFrame has genes as rows and cells as columns.

        Returns True if the index looks like gene names and columns look like
        cell barcodes, meaning the matrix should be transposed so cells are obs.
        """
        import re
        idx_sample = [str(v) for v in df.index[:n_sample]]
        col_sample = [str(v) for v in df.columns[:n_sample]]

        # Gene name: short (<12 chars), starts with uppercase letter or
        # common gene prefix, mostly alphanumeric, no barcode-like patterns.
        _re_gene = re.compile(r'^[A-Z][A-Za-z0-9.\-]{0,11}$')
        # Cell barcode: long (≥8 chars), ACGT-rich.
        _re_bc = re.compile(r'^[ACGT]{8,}')

        n_gene_idx = sum(1 for v in idx_sample if _re_gene.match(v))
        n_bc_idx = sum(1 for v in idx_sample if _re_bc.match(v))
        n_bc_col = sum(1 for v in col_sample if _re_bc.match(v))
        n_gene_col = sum(1 for v in col_sample if _re_gene.match(v))

        # Strong signal: index is gene-like AND columns are barcode-like
        if n_gene_idx > n_sample * 0.8 and n_bc_col > n_sample * 0.8:
            return True
        # Strong signal: index is barcode-like AND columns are gene-like
        # (already oriented as cells×genes)
        if n_bc_idx > n_sample * 0.8 and n_gene_col > n_sample * 0.8:
            return False
        # Ambiguous: fall back to shape heuristic (genes usually > cells)
        if df.shape[0] > df.shape[1]:
            return True
        return False

    def convert(self, raw_dir: Path, dataset_id: str = None) -> "ad.AnnData":
        import anndata as ad

        def _read_table(path: Path) -> pd.DataFrame:
            pstr = str(path)
            opener = gzip.open(pstr, "rt") if pstr.endswith(".gz") else open(pstr, "r")
            with opener as f:
                first = f.read(500)
            # Detect separator: tab > comma > whitespace
            if "\t" in first:
                sep = "\t"
            elif "," in first:
                sep = ","
            else:
                sep = r"\s+"
            return pd.read_csv(pstr, sep=sep, index_col=0, compression="gzip" if pstr.endswith(".gz") else None)

        # ── Identify expression and metadata files ──
        tables = [f for f in raw_dir.iterdir() if not f.name.startswith(".")]
        csv_files = [
            f for f in tables
            if f.suffix.lower() in (".csv", ".tsv", ".txt")
            or (f.suffix.lower() == ".gz" and Path(f.stem).suffix.lower() in (".csv", ".tsv", ".txt"))
        ]
        if not csv_files:
            raise FileNotFoundError(f"No table files found in {raw_dir}")

        expr_keywords = ("count", "matrix", "expression", "expr", "tpm", "fpkm", "rpkm", "processed", "dge")
        meta_keywords = ("metadata", "meta", "annotation")

        expr_files = []
        meta_file = None
        for f in csv_files:
            name = f.name.lower()
            if any(kw in name for kw in meta_keywords):
                meta_file = f
            elif any(kw in name for kw in expr_keywords):
                expr_files.append(f)

        # If no expression files matched by keyword, treat all as expression
        if not expr_files:
            expr_files = list(csv_files)
            # If there's a clear size disparity, the smallest might be metadata
            if len(expr_files) >= 2 and meta_file is None:
                expr_files.sort(key=lambda f: f.stat().st_size, reverse=True)
                sizes = [f.stat().st_size for f in expr_files]
                if len(sizes) == 2 and (sizes[0] > 5 * sizes[1] or sizes[1] > 5 * sizes[0]):
                    meta_file = expr_files.pop(-1)

        # Sort expression files for deterministic ordering
        expr_files.sort(key=lambda f: f.name)
        files_loaded = []

        # ── Load all expression matrices and merge ──
        if len(expr_files) == 1:
            # Single expression matrix
            counts_file = expr_files[0]
            logger.info("Loading expression: %s", counts_file.name)
            df = _read_table(counts_file)
            files_loaded.append(str(counts_file))

            X = sp.csr_matrix(df.values)
            # Detect orientation: if index looks like gene names (short,
            # alphanumeric, uppercase) and columns look like cell barcodes
            # (long ACGT strings), transpose even if genes < cells.
            needs_T = self._detect_genes_as_rows(df)
            if needs_T:
                X = X.T.tocsr()
                obs_names = df.columns.astype(str).tolist()
                var_names = df.index.astype(str).tolist()
            else:
                obs_names = df.index.astype(str).tolist()
                var_names = df.columns.astype(str).tolist()

            # Parse barcode metadata (e.g. "MGH170-P7-A01-CD45pos" → patient, plate, well, marker)
            obs_index = obs_names
        else:
            # Multiple expression matrices — merge by union of genes
            logger.info("Loading %d expression matrices for merge", len(expr_files))
            adatas = []
            for f in expr_files:
                logger.info("  Loading: %s", f.name)
                df = _read_table(f)
                files_loaded.append(str(f))

                X = sp.csr_matrix(df.values)
                if self._detect_genes_as_rows(df):
                    X = X.T.tocsr()
                    cell_names = [f"{f.stem}_{c}" for c in df.columns.astype(str)]
                    gene_names = df.index.astype(str).tolist()
                else:
                    cell_names = [f"{f.stem}_{c}" for c in df.index.astype(str)]
                    gene_names = df.columns.astype(str).tolist()

                a = ad.AnnData(
                    X=X,
                    obs=pd.DataFrame(
                        {"sample_id": f.stem, "source_file": f.name},
                        index=cell_names,
                    ),
                    var=pd.DataFrame(index=list(self.dedup_genes(gene_names))),
                )
                adatas.append(a)
                logger.info("    %d cells, %d genes", a.n_obs, a.n_vars)

            merged = ad.concat(adatas, join="outer", fill_value=0)
            merged.var_names = self.dedup_genes(merged.var_names)
            obs_index = list(merged.obs_names)
            logger.info("Merged: %d cells, %d genes", merged.n_obs, merged.n_vars)
            # Assign merged to adata for metadata join below
            adata = merged
            adata.layers["counts"] = adata.X.copy()
            # Skip the single-file AnnData creation below
            obs_names = None

        if obs_names is not None:
            sid = dataset_id or raw_dir.parent.name
            adata = ad.AnnData(
                X=X,
                obs=pd.DataFrame(
                    {"sample_id": sid, "source_file": expr_files[0].name if expr_files else ""},
                    index=obs_names,
                ),
                var=pd.DataFrame(index=list(self.dedup_genes(var_names))),
            )
            adata.layers["counts"] = adata.X.copy()
            logger.info("Expression: %d cells, %d genes", adata.n_obs, adata.n_vars)
            obs_index = obs_names

        # ── Parse barcode metadata from cell names ──
        self._parse_barcode_metadata(adata, obs_index)

        # ── Join external metadata if available ──
        if meta_file is not None:
            logger.info("Loading metadata: %s", meta_file.name)
            df_meta = _read_table(meta_file)
            files_loaded.append(str(meta_file))

            # Determine barcode column
            # If the current index looks like cell barcodes (ACGT-rich strings),
            # the metadata was already read with index_col=0 and no re-indexing
            # is needed. Otherwise, look for a barcode column in the data.
            _idx_sample = [str(v) for v in df_meta.index[:100]]
            _n_acgt = sum(
                1 for v in _idx_sample
                if __import__("re").match(r'^[ACGT]{8,}', v)
            )
            if _n_acgt > len(_idx_sample) * 0.8:
                # Index is already cell barcodes — keep it
                pass
            elif df_meta.index.name in ("barcode", "cell_barcode", "cell"):
                pass
            else:
                for bc in ["barcode", "cell_barcode", "cell"]:
                    if bc in df_meta.columns:
                        df_meta = df_meta.set_index(bc)
                        break
                else:
                    # Last resort: use first non-index column
                    pass

            meta_barcodes = set(df_meta.index.astype(str))
            obs_barcodes = list(adata.obs_names)

            # Normalize barcode format: expression matrices sometimes use
            # ".N" suffix (e.g. "AAACCTG.1") while metadata uses "-N"
            # (e.g. "AAACCTG-1"). Build a mapping from normalized form
            # to original obs barcode.
            import re as _re_bc
            _bc_norm_map = {}
            for bc in obs_barcodes:
                norm = _re_bc.sub(r"\.(\d+)$", r"-\1", bc)
                _bc_norm_map[norm] = bc
            _meta_norm_set = set(
                _re_bc.sub(r"\.(\d+)$", r"-\1", b) for b in meta_barcodes
            )

            matched = 0
            for col in df_meta.columns:
                if col == "source_file":
                    continue
                if col not in adata.obs.columns:
                    adata.obs[col] = "unknown"
                # Bin matching by normalized barcode
                meta_val_by_bc = {}
                for mb in meta_barcodes:
                    meta_val_by_bc[mb] = str(df_meta.loc[mb, col])
                for norm_bc, obs_bc in _bc_norm_map.items():
                    if norm_bc in meta_barcodes:
                        adata.obs.loc[obs_bc, col] = meta_val_by_bc[norm_bc]
                        matched += 1

            if matched == 0:
                short_map = {}
                for bc in obs_barcodes:
                    short = bc
                    if "_" in bc:
                        short = bc.rsplit("_", 1)[0]
                    short_map[short] = bc

                meta_short_barcodes = {
                    b.rsplit("_", 1)[0]: b
                    for b in meta_barcodes
                    if "_" in b
                }

                for col in df_meta.columns:
                    if col == "source_file":
                        continue
                    if col not in adata.obs.columns:
                        adata.obs[col] = "unknown"
                    for short_bc, full_bc in short_map.items():
                        if short_bc in meta_short_barcodes:
                            meta_idx = meta_short_barcodes[short_bc]
                            if meta_idx in df_meta.index:
                                adata.obs.loc[full_bc, col] = str(df_meta.loc[meta_idx, col])
                                matched += 1

            logger.info("Metadata joined: %d matches for %d columns", matched, len(df_meta.columns))

        # ── Build conversion manifest ──
        adata.uns["raw_conversion"] = self.build_manifest(
            "GeoCSVConverter", str(raw_dir),
            adata.n_obs, adata.n_vars,
            list(adata.obs.columns), files_loaded,
        )
        self.validate_adata(adata)
        logger.info(
            "GeoCSV: %d cells, %d genes (expr_files=%s, meta=%s)",
            adata.n_obs, adata.n_vars,
            [f.name for f in expr_files],
            meta_file.name if meta_file else "none",
        )
        return adata

    def _parse_barcode_metadata(self, adata: "ad.AnnData", obs_names: list):
        """Extract metadata from cell barcode patterns.

        Handles patterns like:
          - "MGH170-P7-A01-CD45pos" → patient=MGH170, plate=P7, well=A01, marker=CD45pos
          - "DMSO_P1_A03_S99" → condition=DMSO, plate=P1, well=A03, sample=S99
          - "AGI881_P1_A01_S1" → condition=AGI881, plate=P1, well=A01, sample=S1

        Barcodes may be embedded within longer obs_names (e.g. with file stem prefix).
        Uses re.search to find the pattern anywhere in the string.
        """
        import re

        # Pattern: PATIENT-PLATE-WELL-MARKER (e.g. MGH170-P7-A01-CD45pos)
        _re_patient = re.compile(
            r'(?P<patient>[A-Za-z]+\d+)-(?P<plate>P\d+)-(?P<well>[A-H]\d+)-(?P<marker>.+)$'
        )
        # Pattern: CONDITION_PLATE_WELL_SAMPLE (e.g. DMSO_P1_A03_S99)
        _re_cond = re.compile(
            r'(?P<condition>[A-Za-z0-9]+)_(?P<plate>P\d+)_(?P<well>[A-H](?:\d+)?)_(?P<sample>S\d+)$'
        )

        for obs_name in obs_names:
            m = _re_patient.search(obs_name)
            if m:
                for k, v in m.groupdict().items():
                    if k not in adata.obs.columns:
                        adata.obs[k] = "unknown"
                    adata.obs.loc[obs_name, k] = v
                continue

            m = _re_cond.search(obs_name)
            if m:
                for k, v in m.groupdict().items():
                    if k not in adata.obs.columns:
                        adata.obs[k] = "unknown"
                    adata.obs.loc[obs_name, k] = v
                continue

        n_patient = (adata.obs.get("patient", None) != "unknown").sum() if "patient" in adata.obs.columns else 0
        n_condition = (adata.obs.get("condition", None) != "unknown").sum() if "condition" in adata.obs.columns else 0
        if n_patient or n_condition:
            logger.info("Barcode metadata parsed: %d patient, %d condition labels", n_patient, n_condition)

        # ── Derive treatment column ──
        # Human samples (MGH170, MGH229, MGH202, BWH*): all treated with IDH inhibitor
        # Mouse samples: AGI881 = treated, DMSO = untreated
        if "treatment" not in adata.obs.columns:
            adata.obs["treatment"] = "unknown"
        if "condition" in adata.obs.columns:
            # Mouse-style barcodes (condition parsed from name)
            adata.obs.loc[adata.obs["condition"].str.upper() == "DMSO", "treatment"] = "untreated"
            adata.obs.loc[adata.obs["condition"].str.upper().str.startswith("AGI"), "treatment"] = "treated"
        if "patient" in adata.obs.columns:
            # Human-style barcodes — all are from treated patients
            human_mask = (
                (adata.obs["patient"] != "unknown") &
                ~adata.obs["patient"].str.lower().str.contains("mouse")
            )
            adata.obs.loc[human_mask, "treatment"] = "treated"
        if "marker" in adata.obs.columns:
            # Also set treatment from marker context (all human samples are treated)
            human_marker_mask = (
                (adata.obs["treatment"] == "unknown") &
                (adata.obs["marker"] != "unknown")
            )
            adata.obs.loc[human_marker_mask, "treatment"] = "treated"
