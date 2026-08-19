"""Raw converter for GSE172138: retinal organoid drug toxicity scRNA-seq.

PMID 35298655. 10x Chromium scRNA-seq of human retinal organoids (day 200)
treated with 6 drugs + PBS control. Data provided as a combined gene×cell
count matrix with cell barcodes suffixed by sample index.

GSM order from GEO series matrix (GSE172138_series_matrix.txt.gz):
  GSM5242603 = AD4_Control (PBS)        → _1
  GSM5242604 = AD4_Digoxin (40nM)       → _2
  GSM5242605 = AD4_Ethanol (500mM)      → _3
  GSM5242606 = AD4_Ketorolac (2.5mM)    → _4
  GSM5242607 = AD4_Methanol (0.4%)      → _5
  GSM5242608 = AD4_Sildenafil (225uM)   → _6
  GSM5242609 = AD4_Thioridazine (135uM)  → _7

Each cell barcode in the counts matrix ends with "_N" where N ∈ {1..7}
identifies the sample it came from.
"""

import gzip
import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .base import BaseRawConverter

logger = logging.getLogger(__name__)

# Suffix → (GSM ID, drug name, dose, condition, treatment)
_SAMPLE_MAP = {
    1: {"gsm_id": "GSM5242603", "drug": "Control", "dose": "PBS",
        "condition": "untreated (PBS)", "treatment": "untreated"},
    2: {"gsm_id": "GSM5242604", "drug": "Digoxin", "dose": "40 nM",
        "condition": "Digoxin", "treatment": "treated"},
    3: {"gsm_id": "GSM5242605", "drug": "Ethanol", "dose": "500 mM",
        "condition": "Ethanol", "treatment": "treated"},
    4: {"gsm_id": "GSM5242606", "drug": "Ketorolac", "dose": "2.5 mM",
        "condition": "Ketorolac", "treatment": "treated"},
    5: {"gsm_id": "GSM5242607", "drug": "Methanol", "dose": "32 mM",
        "condition": "Methanol", "treatment": "treated"},
    6: {"gsm_id": "GSM5242608", "drug": "Sildenafil", "dose": "225 uM",
        "condition": "Sildenafil", "treatment": "treated"},
    7: {"gsm_id": "GSM5242609", "drug": "Thioridazine", "dose": "135 uM",
        "condition": "Thioridazine", "treatment": "treated"},
}


def _parse_suffix(barcode: str) -> int | None:
    """Extract numeric suffix from barcode like 'AAACCCAAGCAAATCA_1'."""
    m = re.search(r"_(\d+)$", str(barcode))
    return int(m.group(1)) if m else None


def _detect_counts_file(raw_dir: Path) -> Path | None:
    """Find the counts matrix file."""
    for pattern in ("*counts*", "*matrix*", "*expression*"):
        for f in sorted(raw_dir.glob(pattern)):
            if f.suffix in (".gz", ".csv", ".tsv", ".txt"):
                return f
    return None


def _detect_umap_file(raw_dir: Path) -> Path | None:
    """Find UMAP embeddings file."""
    for pattern in ("*embedding*", "*umap*"):
        for f in sorted(raw_dir.glob(pattern)):
            if f.suffix in (".gz", ".csv", ".tsv", ".txt"):
                return f
    return None


class GSE172138Converter(BaseRawConverter):
    """GSE172138: retinal organoid drug toxicity — combined counts CSV."""

    @staticmethod
    def detect(raw_dir: Path) -> bool:
        if not raw_dir.is_dir():
            return False
        counts = _detect_counts_file(raw_dir)
        if counts is None:
            return False
        # Quick check: first line should have cell barcodes with _N suffixes
        pstr = str(counts)
        opener = gzip.open(pstr, "rt") if pstr.endswith(".gz") else open(pstr, "r")
        try:
            with opener as f:
                first_line = f.readline()
        except Exception:
            return False
        # Look for barcode pattern: ACGT characters with _N suffix
        barcode_pattern = re.compile(r"[ACGT]{10,}_\d")
        return bool(barcode_pattern.search(first_line))

    def convert(self, raw_dir: Path, dataset_id: str = None) -> "ad.AnnData":
        import anndata as ad

        # ── Load counts matrix ──
        counts_file = _detect_counts_file(raw_dir)
        if counts_file is None:
            raise FileNotFoundError(f"No counts file found in {raw_dir}")

        logger.info("Loading counts: %s (%.1f MB)", counts_file.name,
                     counts_file.stat().st_size / 1e6)

        df = pd.read_csv(
            str(counts_file),
            sep=",",
            index_col=0,
            compression="gzip" if str(counts_file).endswith(".gz") else None,
        )
        # Strip quotes from index and columns
        df.index = df.index.str.strip('"')
        df.columns = df.columns.str.strip('"')

        logger.info("Loaded matrix: %d genes × %d cells", df.shape[0], df.shape[1])

        # Matrix is genes × cells → transpose to cells × genes
        X = sp.csr_matrix(df.values.T)
        cell_barcodes = df.columns.astype(str).tolist()
        gene_names = df.index.astype(str).tolist()

        # ── Build sample metadata from barcode suffixes ──
        obs = pd.DataFrame(index=cell_barcodes)
        obs["cell_barcode"] = cell_barcodes

        for col in ["gsm_id", "drug", "dose", "condition", "treatment"]:
            obs[col] = "unknown"

        for bc in cell_barcodes:
            suffix = _parse_suffix(bc)
            if suffix is not None and suffix in _SAMPLE_MAP:
                meta = _SAMPLE_MAP[suffix]
                for k, v in meta.items():
                    obs.loc[bc, k] = v

        obs["sample_id"] = obs["gsm_id"] + "_" + obs["drug"]
        obs["source_file"] = counts_file.name

        # Set orig.ident for downstream split logic
        obs["orig.ident"] = obs["drug"]

        n_mapped = (obs["treatment"] != "unknown").sum()
        n_drugs = obs["drug"].nunique()
        logger.info("Barcode suffix mapping: %d/%d cells mapped, %d conditions",
                     n_mapped, len(cell_barcodes), n_drugs)
        for drug in sorted(obs["drug"].unique()):
            n = (obs["drug"] == drug).sum()
            logger.info("  %20s: %d cells", drug, n)

        # ── Build AnnData ──
        adata = ad.AnnData(
            X=X,
            obs=obs,
            var=pd.DataFrame(index=list(self.dedup_genes(gene_names))),
        )
        adata.layers["counts"] = adata.X.copy()

        # ── Load UMAP embeddings (optional) ──
        umap_file = _detect_umap_file(raw_dir)
        if umap_file is not None:
            logger.info("Loading UMAP: %s", umap_file.name)
            df_umap = pd.read_csv(
                str(umap_file),
                compression="gzip" if str(umap_file).endswith(".gz") else None,
            )
            df_umap.columns = df_umap.columns.str.strip('"')
            if "cell_barcode" in df_umap.columns:
                df_umap["cell_barcode"] = df_umap["cell_barcode"].str.strip('"')
                # Match by barcode
                umap_map = {}
                for _, row in df_umap.iterrows():
                    bc = str(row["cell_barcode"])
                    umap_map[bc] = row.to_dict()
                for col in ["UMAP_1", "UMAP_2"]:
                    adata.obs[col] = np.nan
                for bc in cell_barcodes:
                    if bc in umap_map:
                        for col in ["UMAP_1", "UMAP_2"]:
                            if col in umap_map[bc]:
                                adata.obs.loc[bc, col] = float(umap_map[bc][col])
                n_matched = adata.obs["UMAP_1"].notna().sum()
                logger.info("UMAP embeddings matched: %d/%d cells", n_matched, len(cell_barcodes))
            # Store UMAP in obsm
            if "UMAP_1" in adata.obs.columns and "UMAP_2" in adata.obs.columns:
                umap_vals = adata.obs[["UMAP_1", "UMAP_2"]].values
                adata.obsm["X_umap"] = umap_vals.astype(np.float32)

        # ── Build conversion manifest ──
        adata.uns["raw_conversion"] = self.build_manifest(
            "GSE172138Converter", str(raw_dir),
            adata.n_obs, adata.n_vars,
            list(adata.obs.columns),
            [str(counts_file)] + ([str(umap_file)] if umap_file else []),
        )
        self.validate_adata(adata)
        logger.info(
            "GSE172138: %d cells, %d genes, %d sample groups",
            adata.n_obs, adata.n_vars, adata.obs["drug"].nunique(),
        )
        return adata
