"""GSE186814 raw converter.

Handles per-sample 10x MEX triplet files with GEO GSM naming:
  GSM5661545_u1_NC_matrix.mtx.gz
  GSM5661545_u1_NC_barcodes.tsv.gz
  GSM5661545_u1_NC_genes.tsv.gz

Each sample is a 10x MEX format with genes.tsv (2-column: gene_id, gene_name).
Metadata is embedded in the sample name: {u1|u2}_{NC|VPA}
"""

import gzip
import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread

from .base import BaseRawConverter

logger = logging.getLogger(__name__)


class GSE186814Converter(BaseRawConverter):
    """Per-sample GEO 10x MEX triplets with GSM prefix."""

    @staticmethod
    def detect(raw_dir: Path) -> bool:
        if not raw_dir.is_dir():
            return False
        # Look for GSM-prefixed matrix.mtx.gz files
        mtx_files = list(raw_dir.glob("GSM*_matrix.mtx.gz"))
        if not mtx_files:
            mtx_files = list(raw_dir.glob("GSM*_matrix.mtx"))
        if len(mtx_files) >= 2:
            return True
        return False

    def convert(self, raw_dir: Path, dataset_id: str = None) -> "ad.AnnData":
        import anndata as ad

        # Find all GSM-prefixed MEX triplets
        mtx_files = sorted(raw_dir.glob("GSM*_matrix.mtx.gz"))
        if not mtx_files:
            mtx_files = sorted(raw_dir.glob("GSM*_matrix.mtx"))

        if not mtx_files:
            raise FileNotFoundError(f"No GSM*_matrix.mtx(.gz) files in {raw_dir}")

        adatas = []
        for mtx_path in mtx_files:
            # Parse sample name from filename: GSM5661545_u1_NC_matrix.mtx.gz → u1_NC
            stem = mtx_path.name
            match = re.match(r'GSM\d+_(.+)_matrix\.mtx(?:\.gz)?$', stem)
            if not match:
                logger.warning("  Could not parse sample name from: %s", stem)
                continue
            sample_name = match.group(1)

            # Find corresponding barcodes and genes files
            barc_path = None
            genes_path = None
            for suffix in ["_barcodes.tsv.gz", "_barcodes.tsv"]:
                candidate = mtx_path.with_name(
                    mtx_path.name.replace("_matrix.mtx.gz", suffix)
                    .replace("_matrix.mtx", suffix)
                )
                if not candidate.exists():
                    # Also try with different base name
                    base = stem.replace("_matrix.mtx.gz", "").replace("_matrix.mtx", "")
                    candidate = mtx_path.parent / f"{base}{suffix}"
                if candidate.exists():
                    barc_path = candidate
                    break

            for suffix in ["_genes.tsv.gz", "_genes.tsv", "_features.tsv.gz", "_features.tsv"]:
                candidate = mtx_path.with_name(
                    mtx_path.name.replace("_matrix.mtx.gz", suffix)
                    .replace("_matrix.mtx", suffix)
                )
                if not candidate.exists():
                    base = stem.replace("_matrix.mtx.gz", "").replace("_matrix.mtx", "")
                    candidate = mtx_path.parent / f"{base}{suffix}"
                if candidate.exists():
                    genes_path = candidate
                    break

            if barc_path is None:
                raise FileNotFoundError(f"No barcodes file for sample {sample_name}")
            if genes_path is None:
                raise FileNotFoundError(f"No genes/features file for sample {sample_name}")

            logger.info("  Loading sample %s: %s", sample_name, mtx_path.name)

            # Read MEX matrix
            if str(mtx_path).endswith(".gz"):
                import io
                with gzip.open(str(mtx_path), "rb") as fh:
                    mtx = mmread(io.BytesIO(fh.read()))
            else:
                mtx = mmread(str(mtx_path))

            # Read barcodes
            barcs = pd.read_csv(str(barc_path), sep="\t", header=None, compression="gzip"
                               if str(barc_path).endswith(".gz") else None)[0].astype(str)
            barcs = [f"{sample_name}_{b}" for b in barcs]

            # Read genes
            genes_df = pd.read_csv(str(genes_path), sep="\t", header=None, compression="gzip"
                                   if str(genes_path).endswith(".gz") else None)
            if genes_df.shape[1] >= 2:
                gene_names = genes_df.iloc[:, 1].fillna(genes_df.iloc[:, 0]).astype(str).tolist()
            else:
                gene_names = genes_df.iloc[:, 0].astype(str).tolist()

            # Ensure genes are unique
            seen = {}
            deduped = []
            for g in gene_names:
                if g not in seen:
                    seen[g] = 0
                    deduped.append(g)
                else:
                    seen[g] += 1
                    deduped.append(f"{g}_{seen[g]}")

            # Transpose MEX: scipy mmread treats it as genes×cells but we want cells×genes
            if mtx.shape[0] == len(gene_names) and mtx.shape[1] == len(barcs):
                X = mtx.T.tocsr()
            elif mtx.shape[1] == len(gene_names) and mtx.shape[0] == len(barcs):
                X = mtx.tocsr()
            else:
                # Guess: more rows → genes
                if mtx.shape[0] > mtx.shape[1]:
                    X = mtx.T.tocsr()
                else:
                    X = mtx.tocsr()

            # Parse metadata from sample name: u1_NC, u1_VPA, u2_NC, u2_VPA
            line = sample_name.split("_")[0]  # u1 or u2
            cond_raw = "_".join(sample_name.split("_")[1:])  # NC or VPA

            is_vpa = "vpa" in cond_raw.lower()
            treatment = "treated" if is_vpa else "untreated"
            drug = "valproic acid" if is_vpa else "none"
            # Map to values matching Case-level control/dose_groups for split logic
            condition = "1.0 mM" if is_vpa else "untreated"

            obs = pd.DataFrame({
                "sample_id": sample_name,
                "source_file": mtx_path.name,
                "cell_line": f"U{line[-1]}M" if line.lower().startswith("u1") else f"U{line[-1]}F",
                "condition": condition,
                "treatment": treatment,
                "drug": drug,
                "orig.ident": condition,
            }, index=barcs)

            a = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=deduped))
            a.layers["counts"] = a.X.copy()
            adatas.append(a)

            logger.info("    %s: %d cells, %d genes (line=%s, condition=%s)",
                        sample_name, a.n_obs, a.n_vars,
                        obs["cell_line"].iloc[0], condition)

        if not adatas:
            raise RuntimeError(f"No valid samples loaded from {raw_dir}")

        # Merge all samples by union of genes
        merged = ad.concat(adatas, join="outer", fill_value=0)
        merged.var_names = self.dedup_genes(merged.var_names)
        merged.layers["counts"] = merged.X.copy()
        merged.uns["raw_conversion"] = self.build_manifest(
            "GSE186814Converter", str(raw_dir), merged.n_obs, merged.n_vars,
            list(merged.obs.columns), [str(f) for f in mtx_files],
        )
        self.validate_adata(merged)
        logger.info("Merged %d samples: %d cells, %d genes",
                    len(adatas), merged.n_obs, merged.n_vars)
        return merged
