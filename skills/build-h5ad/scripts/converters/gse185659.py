"""Converter for GSE185659: methylprednisolone-treated lung transplant BAL T cells.

Experimental design (PMID 35285873):
  - scRNA-seq + TCR-seq of recipient-derived CD3+ T cells from BAL
  - 3 participants (P1, P3, P8), each with ACR (rejection) and Treated (post-glucocorticoid)
  - P8 also has an Early (pre-ACR) sample
  - 10x Chromium Single Cell 5' v1, sequenced on Illumina HiSeq 2500
  - Alignment with CellRanger v5

Data format: 10x HDF5 per sample (*_raw_feature_bc_matrix.h5).

Benchmark split:
  - control: ACR samples (P1_ACR, P3_ACR, P8_ACR) — rejection, pre-treatment
  - ground_truth: Treated samples (P1_Treated, P3_Treated, P8_Treated) — post-treatment
  - P8_Early excluded (pre-ACR, not part of ACR vs post-treatment comparison)

perturb_var = time, time_groups = ["60 day"] (~2 months post-ACR primary follow-up).
"""

import logging
import re
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

from .base import BaseConverter

logger = logging.getLogger(__name__)

# Per-sample metadata derived from GEO records and paper
_SAMPLE_META = {
    "GSM5621344": {"participant": "P1", "condition": "ACR",     "label": "P1_ACR"},
    "GSM5621345": {"participant": "P1", "condition": "Treated", "label": "P1_Treated"},
    "GSM5621346": {"participant": "P3", "condition": "ACR",     "label": "P3_ACR"},
    "GSM5621347": {"participant": "P3", "condition": "Treated", "label": "P3_Treated"},
    "GSM5621348": {"participant": "P8", "condition": "Early",   "label": "P8_Early"},
    "GSM5621349": {"participant": "P8", "condition": "ACR",     "label": "P8_ACR"},
    "GSM5621350": {"participant": "P8", "condition": "Treated", "label": "P8_Treated"},
}

# Follow-up times per participant (paper line 35): P1=2mo, P3=6mo, P8=2mo post-ACR
_PARTICIPANT_FOLLOWUP = {"P1": "60 day", "P3": "180 day", "P8": "60 day"}


def _parse_gsm_id(filename: str) -> str | None:
    m = re.match(r"(GSM\d+)", filename)
    return m.group(1) if m else None


def _read_10x_h5(path: Path) -> tuple[sp.csr_matrix, list[str], list[str]]:
    """Read a 10x Genomics HDF5 file and return (matrix, barcodes, gene_symbols)."""
    with h5py.File(path, "r") as f:
        mg = f["matrix"]
        shape = tuple(mg["shape"][:])
        data = np.array(mg["data"], dtype=np.float32)
        indices = np.array(mg["indices"], dtype=np.int32)
        indptr = np.array(mg["indptr"], dtype=np.int32)

        # 10x HDF5 stores matrix in CSC orientation: (genes, barcodes)
        mtx = sp.csc_matrix((data, indices, indptr), shape=shape)

        barcodes = [b.decode("utf-8") for b in mg["barcodes"][:]]
        gene_symbols = [g.decode("utf-8") for g in f["matrix/features/name"][:]]

    return mtx, barcodes, gene_symbols


class GSE185659Converter(BaseConverter):
    """GSE185659: methylprednisolone-treated lung transplant BAL CD8+ TRM cells."""

    @staticmethod
    def detect(raw_dir: Path) -> bool:
        if not raw_dir.is_dir():
            return False
        h5_files = list(raw_dir.glob("*_filtered_feature_bc_matrix.h5"))
        if len(h5_files) < 3:
            return False
        # Check for characteristic GSM IDs
        gsm_ids = set()
        for f in h5_files:
            gsm = _parse_gsm_id(f.name)
            if gsm:
                gsm_ids.add(gsm)
        expected = {"GSM5621344", "GSM5621345", "GSM5621346", "GSM5621347",
                     "GSM5621348", "GSM5621349", "GSM5621350"}
        return len(gsm_ids & expected) >= 3

    def convert(self, raw_dir: Path) -> tuple[ad.AnnData, ad.AnnData]:
        h5_files = sorted(raw_dir.glob("*_filtered_feature_bc_matrix.h5"))
        logger.info("Found %d raw h5 file(s)", len(h5_files))

        adatas = []
        gene_set = None

        for h5_path in h5_files:
            gsm_id = _parse_gsm_id(h5_path.name)
            if gsm_id is None:
                logger.warning("  Skipping %s: cannot parse GSM ID", h5_path.name)
                continue

            meta = _SAMPLE_META.get(gsm_id)
            if meta is None:
                logger.warning("  Skipping unknown sample %s", gsm_id)
                continue

            logger.info("  Loading %s: %s (%s)", gsm_id, meta["label"], meta["condition"])

            mtx, barcodes, gene_symbols = _read_10x_h5(h5_path)

            # Deduplicate gene symbols
            seen = {}
            deduped = []
            for g in gene_symbols:
                g_str = str(g)
                if g_str in seen:
                    seen[g_str] += 1
                    deduped.append(f"{g_str}_{seen[g_str]}")
                else:
                    seen[g_str] = 1
                    deduped.append(g_str)

            # 10x HDF5 is (genes, barcodes) — transpose to (cells, genes)
            if mtx.shape[0] == len(gene_symbols) and mtx.shape[1] == len(barcodes):
                mtx = mtx.T.tocsr()

            cell_ids = [f"{meta['label']}_{b.rstrip('-1')}" for b in barcodes]

            # Build obs with metadata
            followup = _PARTICIPANT_FOLLOWUP.get(meta["participant"], "")
            obs = pd.DataFrame({
                "sample_id": meta["label"],
                "gsm_id": gsm_id,
                "participant": meta["participant"],
                "condition": meta["condition"],
                "time": followup if meta["condition"] == "Treated" else "",
            }, index=cell_ids)

            a = ad.AnnData(X=mtx, obs=obs, var=pd.DataFrame(index=deduped))
            a.var_names = deduped

            if gene_set is None:
                gene_set = set(deduped)
            elif set(deduped) != gene_set:
                logger.warning(
                    "Gene set mismatch for %s: %d vs expected %d",
                    gsm_id, len(deduped), len(gene_set),
                )

            adatas.append(a)

        if not adatas:
            raise ValueError("No samples to process")

        adata = ad.concat(adatas, join="outer", fill_value=0)
        logger.info("Combined: %d cells, %d genes", adata.n_obs, adata.n_vars)

        # Report cell counts per condition
        for cond in ["ACR", "Treated", "Early"]:
            n = (adata.obs["condition"] == cond).sum()
            logger.info("  %s: %d cells", cond, n)

        # Split: ACR → control, Treated → ground_truth, Early excluded
        control_mask = adata.obs["condition"] == "ACR"
        treated_mask = adata.obs["condition"] == "Treated"

        logger.info("Control (ACR): %d cells, Treated: %d cells",
                     control_mask.sum(), treated_mask.sum())

        control = adata[control_mask].copy()
        treated = adata[treated_mask].copy()

        self.normalize(control)
        self.normalize(treated)

        return control, treated
