"""Converter for GSE157526: SARS-CoV-2 infected primary airway epithelium ± remdesivir.

Experimental design:
  S1 (GSM4769386): untreated, uninfected → excluded
  S2 (GSM4769387): untreated, SARS-CoV-2, 24h → control
  S3 (GSM4769388): remdesivir, SARS-CoV-2, 24h → treated
  S5 (GSM4769389): untreated, SARS-CoV-2, 48h → control
  S6 (GSM4769390): remdesivir, SARS-CoV-2, 48h → treated

Data format: 10x MEX per sample (matrix.mtx.gz + barcodes.tsv.gz + features.tsv.gz).

Cell types are identified via airway epithelial marker scoring:
  - Ciliated: FOXJ1, SNTN, TUBB4B, etc.
  - Secretory: SCGB1A1, SCGB3A1, CEACAM6, etc.
  - Basal: KRT5, TP63, NGFR, etc.
  - Ionocyte: FOXI1, CFTR, ASCL3, etc.
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

from .base import BaseConverter

logger = logging.getLogger(__name__)

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

# Airway epithelial marker genes (Plasschaert et al 2018 + paper markers)
_CELL_TYPE_MARKERS = {
    "ciliated cells": [
        "FOXJ1", "SNTN", "TPPP3", "RSPH1", "DNAH5", "DNAH9",
        "CCDC153", "CFAP45", "CFAP54", "DYNLRB2", "TUBB4B", "CAPS",
    ],
    "secretory cells": [
        "SCGB1A1", "SCGB3A1", "MUC5AC", "MUC5B", "CEACAM6",
        "SLPI", "MSMB", "PIGR", "BPIFA1", "WFDC2",
    ],
    "basal cells": [
        "KRT5", "KRT14", "KRT15", "TP63", "NGFR", "PDPN", "KRT17",
    ],
    "ionocyte cells": [
        "FOXI1", "CFTR", "ASCL3",
    ],
}


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


def _read_tsv_gz(path: Path) -> pd.Series:
    if path.suffix == ".gz":
        df = pd.read_csv(path, sep="\t", header=None, compression="gzip")
    else:
        df = pd.read_csv(path, sep="\t", header=None)
    return df.iloc[:, 0]


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


def _assign_cell_types(adata: "ad.AnnData") -> pd.DataFrame:
    """Assign cell types via marker gene mean-expression scoring."""
    import scanpy as sc

    temp = adata.copy()
    temp.layers["counts"] = temp.X.copy()
    sc.pp.normalize_total(temp, target_sum=1e4)
    sc.pp.log1p(temp)

    scores = {}
    present_markers = {}
    for ct_name, markers in _CELL_TYPE_MARKERS.items():
        found = [g for g in markers if g in temp.var_names]
        if not found:
            continue
        present_markers[ct_name] = found
        idx = [temp.var_names.get_loc(g) for g in found]
        if sp.issparse(temp.X):
            scores[ct_name] = np.array(temp.X[:, idx].mean(axis=1)).ravel()
        else:
            scores[ct_name] = temp.X[:, idx].mean(axis=1)

    if not scores:
        return pd.DataFrame({
            "cell_type": "undefined",
            "annotation_method": "unassigned",
        }, index=adata.obs.index)

    ct_names = list(scores.keys())
    score_matrix = np.column_stack([scores[ct] for ct in ct_names])
    best_idx = np.argmax(score_matrix, axis=1)
    best_score = score_matrix[np.arange(len(best_idx)), best_idx]
    assigned = np.where(best_score >= 0.15, best_idx, -1)

    cell_types = [
        ct_names[assigned[i]] if assigned[i] >= 0 else "undefined"
        for i in range(len(assigned))
    ]

    counts = pd.Series(cell_types).value_counts()
    logger.info("Cell type assignment: %s", counts.to_dict())
    for ct_name in ct_names:
        logger.info("  %s markers (%d): %s",
                     ct_name, len(present_markers.get(ct_name, [])),
                     ", ".join(present_markers.get(ct_name, [])[:6]))

    return pd.DataFrame({
        "cell_type": cell_types,
        "annotation_method": "marker_scoring",
    }, index=adata.obs.index)


class GSE157526Converter(BaseConverter):
    """GSE157526: remdesivir-treated SARS-CoV-2 infected airway epithelium."""

    @staticmethod
    def detect(raw_dir: Path) -> bool:
        if not raw_dir.is_dir():
            return False
        samples = _find_sample_files(raw_dir)
        if not samples:
            return False
        gsm_ids = set(samples.keys())
        expected = {"GSM4769386", "GSM4769387", "GSM4769388", "GSM4769389", "GSM4769390"}
        return len(gsm_ids & expected) >= 3

    def convert(self, raw_dir: Path) -> tuple[ad.AnnData, ad.AnnData]:
        samples = _find_sample_files(raw_dir)
        logger.info("Found %d sample(s) with MEX triplets: %s",
                     len(samples), sorted(samples.keys()))

        adatas = []
        gene_set = None

        for gsm_id in sorted(samples.keys()):
            meta = _SAMPLE_META.get(gsm_id)
            if meta is None:
                logger.warning("  Skipping unknown sample %s", gsm_id)
                continue
            if meta["condition"] == "exclude":
                logger.info("  Skipping %s (%s): %s, %s",
                            meta["sample"], gsm_id, meta["treatment"], meta["infection"])
                continue

            files = samples[gsm_id]
            logger.info("  Loading %s (%s): %s, %s, %s",
                         meta["sample"], gsm_id, meta["treatment"],
                         meta["infection"], meta["time"])

            mtx = _read_mtx_gz(files["matrix"])
            barcodes = _read_tsv_gz(files["barcodes"]).astype(str)
            features = _read_features(files["features"]).astype(str)

            seen = {}
            deduped = []
            for g in features:
                g_str = str(g)
                if g_str in seen:
                    seen[g_str] += 1
                    deduped.append(f"{g_str}_{seen[g_str]}")
                else:
                    seen[g_str] = 1
                    deduped.append(g_str)

            mtx = sp.csc_matrix(mtx)
            if mtx.shape[0] == len(features) and mtx.shape[1] == len(barcodes):
                mtx = mtx.T.tocsr()

            cell_ids = [f"{meta['sample']}_{b}" for b in barcodes]

            obs = pd.DataFrame({
                "sample_id": meta["sample"],
                "gsm_id": gsm_id,
                "treatment": meta["treatment"],
                "infection": meta["infection"],
                "time": meta["time"],
                "condition": meta["condition"],
            }, index=cell_ids)

            a = ad.AnnData(X=mtx, obs=obs, var=pd.DataFrame(index=deduped))
            a.var_names = deduped

            if gene_set is None:
                gene_set = set(deduped)
            elif set(deduped) != gene_set:
                logger.warning(
                    "Gene set mismatch for %s: %d vs %d genes",
                    gsm_id, len(deduped), len(gene_set),
                )

            adatas.append(a)

        if not adatas:
            raise ValueError("No samples to process after filtering")

        adata = ad.concat(adatas, join="outer", fill_value=0)
        logger.info("Combined: %d cells, %d genes", adata.n_obs, adata.n_vars)

        # Marker-based cell type assignment
        ct_labels = _assign_cell_types(adata)
        adata.obs["cell_type"] = ct_labels["cell_type"]
        adata.obs["annotation_method"] = ct_labels["annotation_method"]
        adata.obs["cell_type_raw"] = adata.obs["cell_type"]
        adata.obs["cell_type_mapping_method"] = "marker_scoring"

        # Split by condition
        control_mask = adata.obs["condition"] == "control"
        treated_mask = adata.obs["condition"] == "treated"

        logger.info("Control: %d cells, Treated: %d cells, cell types in control: %s, treated: %s",
                     control_mask.sum(), treated_mask.sum(),
                     adata[control_mask].obs["cell_type"].value_counts().to_dict(),
                     adata[treated_mask].obs["cell_type"].value_counts().to_dict())

        control = adata[control_mask].copy()
        treated = adata[treated_mask].copy()

        self.normalize(control)
        self.normalize(treated)

        return control, treated
