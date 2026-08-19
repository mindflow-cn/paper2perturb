"""E-MTAB-9154 raw converter.

KOLF-2 iPSC-derived dopamine neurons treated with rotenone / tunicamycin.
6 TSV files: WT_UNT, WT_ROT, WT_TUN, HE_UNT, HE_ROT, HE_TUN.
Parses genotype (WT/HE) and treatment (UNT/ROT/TUN) from filenames.
Attempts to parse SDRF metadata.
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .base import BaseRawConverter

logger = logging.getLogger(__name__)


class EMTAB9154RawConverter(BaseRawConverter):
    """E-MTAB-9154: 6 TSV expression matrices + SDRF metadata."""

    @staticmethod
    def detect(raw_dir: Path) -> bool:
        if not raw_dir.is_dir():
            return False
        # Check for characteristic TSV files
        required = {"WT_UNT.tsv", "WT_ROT.tsv", "WT_TUN.tsv",
                     "HE_UNT.tsv", "HE_ROT.tsv", "HE_TUN.tsv"}
        existing = {f.name for f in raw_dir.iterdir() if f.suffix == ".tsv"}
        return required.issubset(existing)

    def convert(self, raw_dir: Path, dataset_id: str = None) -> "ad.AnnData":
        import anndata as ad

        # Parse treatment/genotype mapping
        tsv_files = sorted([
            f for f in raw_dir.iterdir()
            if f.suffix == ".tsv" and not f.name.startswith(".")
        ])

        if not tsv_files:
            raise FileNotFoundError(f"No TSV files found in {raw_dir}")

        adatas = []
        files_loaded = []

        for tp in tsv_files:
            name = tp.stem  # e.g. "WT_UNT"
            parts = name.split("_")
            genotype = parts[0] if len(parts) >= 1 else "unknown"
            treatment = parts[1] if len(parts) >= 2 else "unknown"

            # Read expression matrix
            df = pd.read_csv(str(tp), sep="\t", index_col=0)
            X = sp.csr_matrix(df.values)
            if df.shape[0] > df.shape[1]:
                X = X.T
                obs_names = [f"{name}_{c}" for c in df.columns.astype(str)]
                var_names = df.index.astype(str)
            else:
                obs_names = [f"{name}_{c}" for c in df.index.astype(str)]
                var_names = df.columns.astype(str)

            a = ad.AnnData(
                X=X.tocsr(),
                obs=pd.DataFrame({
                    "sample_id": name,
                    "source_file": tp.name,
                    "genotype": genotype,
                    "treatment": treatment,
                    "condition": treatment,
                    "dataset_id": dataset_id or "E-MTAB-9154",
                }, index=obs_names),
                var=pd.DataFrame(index=list(self.dedup_genes(var_names))),
            )
            adatas.append(a)
            files_loaded.append(str(tp))
            logger.info("  Loaded %s: %d cells, %d genes (genotype=%s, treatment=%s)",
                       name, a.n_obs, a.n_vars, genotype, treatment)

        # Merge all
        merged = ad.concat(adatas, join="outer", fill_value=0)
        merged.var_names = self.dedup_genes(merged.var_names)
        merged.layers["counts"] = merged.X.copy()

        # Try SDRF metadata
        sdrf = raw_dir / "E-MTAB-9154.sdrf.txt"
        metadata_parsed = {}
        if sdrf.exists():
            try:
                sdrf_df = pd.read_csv(str(sdrf), sep="\t")
                metadata_parsed = {
                    "sdrf_columns": list(sdrf_df.columns),
                    "sdrf_rows": len(sdrf_df),
                }
                # Map SDRF source names to our sample_ids
                # A53T = HE (heterozygous A53T mutation = TRIO+/-)
                sdrf_to_sample = {
                    "a53t_rot": "HE_ROT", "a53t_tun": "HE_TUN", "a53t_unt": "HE_UNT",
                    "wt_rot": "WT_ROT", "wt_tun": "WT_TUN", "wt_unt": "WT_UNT",
                }
                if "Source Name" in sdrf_df.columns:
                    # First, create all SDRF columns in obs (initialized to empty string)
                    sdrf_columns_to_add = [
                        c for c in sdrf_df.columns
                        if c != "Source Name" and c not in merged.obs.columns
                    ]
                    for col in sdrf_columns_to_add:
                        merged.obs[col] = ""

                    # Then fill per-sample
                    for _, row in sdrf_df.iterrows():
                        src = str(row["Source Name"]).strip().lower()
                        # Find matching sample_id
                        matched_sid = None
                        for src_key, sid in sdrf_to_sample.items():
                            if src_key in src or src in src_key:
                                matched_sid = sid
                                break
                        if matched_sid is None:
                            for sid in merged.obs["sample_id"].unique():
                                if src in sid.lower() or sid.lower() in src:
                                    matched_sid = sid
                                    break
                        if matched_sid is not None:
                            mask = merged.obs["sample_id"] == matched_sid
                            for col in sdrf_columns_to_add:
                                merged.obs.loc[mask, col] = str(row[col])
                    logger.info("Parsed SDRF: %d rows, %d columns added",
                               len(sdrf_df), len(sdrf_columns_to_add))
                else:
                    self.warnings.append("SDRF has no 'Source Name' column")
            except Exception as e:
                self.warnings.append(f"SDRF parse failed: {e}")

        merged.uns["raw_conversion"] = self.build_manifest(
            "EMTAB9154RawConverter", str(raw_dir),
            merged.n_obs, merged.n_vars,
            list(merged.obs.columns), files_loaded, metadata_parsed,
        )
        self.validate_adata(merged)
        logger.info("E-MTAB-9154: %d cells, %d genes from %d files",
                   merged.n_obs, merged.n_vars, len(tsv_files))
        return merged
