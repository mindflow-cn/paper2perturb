"""Raw converter for GSE157526: per-sample 10x MEX with GSM-prefixed filenames
+ marker-based cell type identification.

Cell types from Plasschaert et al 2018 airway epithelium markers:
  - Ciliated cells: FOXJ1+, SNTN+, TUBB4B+, SiR-Tubulin+
  - Secretory cells: SCGB1A1+, SCGB3A1+, CEACAM6+, CD66c+
  - Basal cells: KRT5+, TP63+, NGFR+, CD271+
  - Ionocyte cells: FOXI1+, CFTR+, ASCL3+

The paper (PMID 33507952) uses flow cytometry + scRNA-seq to resolve these
subtypes and reports cell-type-specific drug responses (e.g. IL-6 in secretory
cells, MT1F in ciliated cells).  Accurate cell_type labels are essential for
the downstream annotation pipeline to select the correct population.
"""

import gzip
import logging
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread

from .base import BaseRawConverter

logger = logging.getLogger(__name__)

# ── Sample metadata ──────────────────────────────────────────────
_SAMPLE_META = {
    "GSM4769386": {"sample": "S1", "treatment": "untreated", "infection": "uninfected",
                   "time": "N/A", "condition": "exclude"},
    "GSM4769387": {"sample": "S2", "treatment": "untreated", "infection": "SARS-CoV-2",
                   "time": "24h", "condition": "control"},
    "GSM4769388": {"sample": "S3", "treatment": "remdesivir", "infection": "SARS-CoV-2",
                   "time": "24h", "condition": "treated"},
    "GSM4769389": {"sample": "S5", "treatment": "untreated", "infection": "SARS-CoV-2",
                   "time": "48h", "condition": "control"},
    "GSM4769390": {"sample": "S6", "treatment": "remdesivir", "infection": "SARS-CoV-2",
                   "time": "48h", "condition": "treated"},
}

# ── Airway epithelial cell type marker genes ─────────────────────
# Based on Plasschaert et al (2018) Nature + paper-specific markers.
_CELL_TYPE_MARKERS = {
    "ciliated cells": [
        "FOXJ1", "SNTN", "TPPP3", "RSPH1", "RSPH4A", "RSPH9",
        "DNAH5", "DNAH9", "DNAH11", "CCDC153", "CFAP45", "CFAP54",
        "DYNLRB2", "TEKT1", "ROPN1L", "TUBB4B", "CAPS",
    ],
    "secretory cells": [
        "SCGB1A1", "SCGB3A1", "MUC5AC", "MUC5B", "CEACAM6",
        "SLPI", "MSMB", "PIGR", "BPIFA1", "SCGB1C1", "LCN2",
        "WFDC2", "CYP2F1",
    ],
    "basal cells": [
        "KRT5", "KRT14", "KRT15", "TP63", "NGFR",
        "PDPN", "DAPL1", "KRT17",
    ],
    "ionocyte cells": [
        "FOXI1", "CFTR", "ASCL3",
    ],
}


# ── Helpers ──────────────────────────────────────────────────────

def _parse_sample_name(filename: str) -> str | None:
    m = re.match(r"(GSM\d+)", filename)
    return m.group(1) if m else None


def _find_sample_files(raw_dir: Path) -> dict[str, dict[str, Path]]:
    samples = {}
    for f in raw_dir.iterdir():
        gsm = _parse_sample_name(f.name)
        if gsm is None:
            continue
        if gsm not in samples:
            samples[gsm] = {}
        fname = f.name.lower()
        if "matrix" in fname:
            samples[gsm]["matrix"] = f
        elif "barcode" in fname:
            samples[gsm]["barcodes"] = f
        elif "feature" in fname or "gene" in fname:
            samples[gsm]["features"] = f
    return {k: v for k, v in samples.items()
            if all(key in v for key in ("matrix", "barcodes", "features"))}


def _read_mtx_gz(path: Path):
    if path.suffix == ".gz":
        with gzip.open(path, "rb") as fh:
            return mmread(fh)
    return mmread(str(path))


def _read_tsv_gz_col0(path: Path) -> pd.Series:
    if path.suffix == ".gz":
        df = pd.read_csv(path, sep="\t", header=None, compression="gzip")
    else:
        df = pd.read_csv(path, sep="\t", header=None)
    return df.iloc[:, 0].astype(str)


def _read_features(path: Path) -> pd.Series:
    if path.suffix == ".gz":
        df = pd.read_csv(path, sep="\t", header=None, compression="gzip")
    else:
        df = pd.read_csv(path, sep="\t", header=None)
    if df.shape[1] >= 2:
        symbols = df.iloc[:, 1].astype(str)
        n_ensg_col0 = df.iloc[:, 0].astype(str).str.startswith("ENSG").sum()
        n_ensg_col1 = symbols.str.startswith("ENSG").sum()
        if n_ensg_col1 > n_ensg_col0 and n_ensg_col0 == 0:
            return df.iloc[:, 0].astype(str)
        return symbols
    return df.iloc[:, 0].astype(str)


def _dedup_genes(gene_names: list) -> list:
    seen = {}
    result = []
    for g in gene_names:
        g_str = str(g).strip()
        if g_str in seen:
            seen[g_str] += 1
            result.append(f"{g_str}_{seen[g_str]}")
        else:
            seen[g_str] = 0
            result.append(g_str)
    return result


def _assign_cell_types(adata: "ad.AnnData") -> pd.DataFrame:
    """Assign cell types via marker gene scoring.

    For each cell type, computes a mean-expression score across its marker
    genes (using log1p(CPM) normalised data), then assigns each cell to the
    highest-scoring type.  Cells where no score exceeds a minimum threshold
    are labelled "undefined".
    """
    import scanpy as sc

    # Work on a temporary normalised copy for scoring
    temp = adata.copy()
    temp.layers["counts"] = temp.X.copy()
    sc.pp.normalize_total(temp, target_sum=1e4)
    sc.pp.log1p(temp)

    # Score each cell type
    scores = {}
    present_markers = {}
    for ct_name, marker_list in _CELL_TYPE_MARKERS.items():
        found = [g for g in marker_list if g in temp.var_names]
        if not found:
            logger.warning("No markers found for %s — skipping", ct_name)
            continue
        present_markers[ct_name] = found
        marker_idx = [temp.var_names.get_loc(g) for g in found]
        if sp.issparse(temp.X):
            scores[ct_name] = np.array(temp.X[:, marker_idx].mean(axis=1)).ravel()
        else:
            scores[ct_name] = temp.X[:, marker_idx].mean(axis=1)

    if not scores:
        logger.warning("No cell type markers found in data — all cells labelled 'undefined'")
        result = pd.DataFrame({
            "cell_type": "undefined",
            "annotation_method": "unassigned",
        }, index=adata.obs.index)
        return result

    # Build score matrix: cells × cell_types
    ct_names = list(scores.keys())
    score_matrix = np.column_stack([scores[ct] for ct in ct_names])

    # Assign each cell to the highest-scoring type, with a minimum threshold
    best_idx = np.argmax(score_matrix, axis=1)
    best_score = score_matrix[np.arange(len(best_idx)), best_idx]

    MIN_SCORE = 0.15  # minimum mean log1p(CPM) to assign; below → "undefined"
    assigned = np.where(best_score >= MIN_SCORE, best_idx, -1)

    cell_types = []
    for i in range(len(assigned)):
        if assigned[i] >= 0:
            cell_types.append(ct_names[assigned[i]])
        else:
            cell_types.append("undefined")

    result = pd.DataFrame({
        "cell_type": cell_types,
        "annotation_method": "marker_scoring",
    }, index=adata.obs.index)

    # Log distribution
    counts = pd.Series(cell_types).value_counts()
    logger.info("Cell type assignment via marker scoring: %s", counts.to_dict())
    for ct_name in ct_names:
        markers_used = present_markers.get(ct_name, [])
        logger.info("  %s markers (%d): %s", ct_name, len(markers_used),
                     ", ".join(markers_used[:8]))

    return result


# ── Converter ────────────────────────────────────────────────────

class GSE157526RawConverter(BaseRawConverter):
    """Per-sample 10x MEX format with marker-based airway epithelial cell typing."""

    @staticmethod
    def detect(raw_dir: Path) -> bool:
        if not raw_dir.is_dir():
            return False
        samples = _find_sample_files(raw_dir)
        if not samples:
            return False
        gsm_ids = set(samples.keys())
        expected = {"GSM4769386", "GSM4769387", "GSM4769388", "GSM4769389", "GSM4769390"}
        if len(gsm_ids & expected) >= 3:
            return True
        return len(samples) >= 2

    def convert(self, raw_dir: Path, dataset_id: str = None) -> ad.AnnData:
        samples = _find_sample_files(raw_dir)
        logger.info("GSE157526 raw: found %d sample(s): %s",
                     len(samples), sorted(samples.keys()))

        adatas = []

        for gsm_id in sorted(samples.keys()):
            meta = _SAMPLE_META.get(gsm_id)
            if meta is None:
                meta = {"sample": gsm_id, "treatment": "unknown",
                        "infection": "unknown", "time": "unknown",
                        "condition": "unknown"}

            if meta["condition"] == "exclude":
                logger.info("  Skipping %s (%s): %s, %s",
                             meta["sample"], gsm_id, meta["treatment"], meta["infection"])
                continue

            files = samples[gsm_id]
            logger.info("  Loading %s (%s): %s, %s, %s",
                         meta["sample"], gsm_id, meta["treatment"],
                         meta["infection"], meta["time"])

            mtx = _read_mtx_gz(files["matrix"])
            barcodes = _read_tsv_gz_col0(files["barcodes"])
            features = _read_features(files["features"])
            deduped = _dedup_genes(features)

            # Orient to cells × genes
            mtx = sp.csc_matrix(mtx)
            if mtx.shape[0] == len(deduped) and mtx.shape[1] == len(barcodes):
                mtx = mtx.T.tocsr()
            else:
                mtx = mtx.tocsr()

            cell_ids = [f"{meta['sample']}_{b}" for b in barcodes]

            obs = pd.DataFrame({
                "sample_id": meta["sample"],
                "source_file": gsm_id,
                "treatment": meta["treatment"],
                "infection": meta["infection"],
                "time": meta["time"],
                "condition": meta["condition"],
            }, index=cell_ids)

            a = ad.AnnData(X=mtx, obs=obs,
                           var=pd.DataFrame(index=deduped))
            a.var_names = deduped
            a.layers["counts"] = a.X.copy()
            a.obs["dataset_id"] = dataset_id or "GSE157526"
            adatas.append(a)

        if not adatas:
            raise ValueError("No samples found after filtering")

        adata = ad.concat(adatas, join="outer", fill_value=0)
        logger.info("Combined: %d cells, %d genes", adata.n_obs, adata.n_vars)

        # ── Marker-based cell type assignment ──
        ct_labels = _assign_cell_types(adata)
        adata.obs["cell_type"] = ct_labels["cell_type"]
        adata.obs["annotation_method"] = ct_labels["annotation_method"]
        adata.obs["cell_type_raw"] = adata.obs["cell_type"]
        adata.obs["cell_type_mapping_method"] = "marker_scoring"

        # ── Build conversion manifest ──
        adata.uns["raw_conversion"] = {
            "converter": "GSE157526RawConverter",
            "input_dir": str(raw_dir),
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "obs_columns": list(adata.obs.columns),
            "files_loaded": sorted(
                str(f) for f in raw_dir.iterdir() if f.suffix in (".gz", ".tar")
            ),
            "cell_type_markers": {
                ct: genes for ct, genes in _CELL_TYPE_MARKERS.items()
            },
            "metadata_parsed": {},
            "warnings": [],
        }

        return adata
