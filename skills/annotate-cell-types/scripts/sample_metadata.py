"""Sample metadata enrichment.

Parses metadata from GSM/sample filenames and optionally user-provided
tables (GEO metadata / sample sheet). Injects columns into adata.obs
that downstream split logic needs: time, condition, dose, patient_id, etc.

Key principles:
- L/NL maps to lesion_status (is_lesional), NOT control/treated.
- is_control is only set when Case-level or user metadata explicitly defines it.
- Visit-to-real-time mapping is configurable via JSON file or CLI flag.
"""

import json
import logging
import re
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

# Default visit-time mappings. Override via --visit-time-mapping JSON file.
DEFAULT_VISIT_MAPPINGS = [
    (re.compile(r"V(?P<v>\d+)"), {
        1: "Day0", 2: "3 day", 3: "14 day",
        4: "28 day", 5: "56 day",
    }),
    (re.compile(r"D(?P<v>\d+)"), {
        0: "Day0", 1: "Day1", 3: "Day3",
        6: "Day6", 7: "Day7", 10: "Day10",
        14: "Day14", 21: "Day21", 28: "Day28",
    }),
    (re.compile(r"(?P<v>\d+)\s*(h|hr|hour)", re.IGNORECASE), None),
    (re.compile(r"(?P<v>\d+)\s*(d|day)", re.IGNORECASE), None),
    (re.compile(r"(?P<v>\d+)\s*(w|wk|week)", re.IGNORECASE), None),
]

_active_visit_mappings = list(DEFAULT_VISIT_MAPPINGS)
_mapping_source = "default"


def load_visit_mapping(config_path: str):
    """Load visit-time mapping from a JSON config file.

    Format: [{"pattern": "V(\\d+)", "map": {"1": "Day0", "2": "3 day"}}, ...]
    """
    global _active_visit_mappings, _mapping_source
    with open(config_path) as f:
        configs = json.load(f)
    new_mappings = []
    for cfg in configs:
        pat = re.compile(cfg["pattern"])
        mp = {int(k): v for k, v in cfg.get("map", {}).items()}
        new_mappings.append((pat, mp if mp else None))
    _active_visit_mappings = new_mappings
    _mapping_source = config_path
    logger.info("Loaded visit-time mapping from %s (%d patterns)", config_path, len(new_mappings))


def get_mapping_source() -> str:
    return _mapping_source


def _map_visit_to_time(visit_label: str) -> str:
    for pattern, mapping in _active_visit_mappings:
        m = pattern.match(str(visit_label))
        if m:
            try:
                v = int(m.group("v"))
            except (ValueError, IndexError):
                return str(visit_label)
            if mapping and v in mapping:
                return mapping[v]
            return str(visit_label)
    return str(visit_label)


LESION_KEYWORDS = {
    "l": True, "lesional": True, "lesion": True,
    "nl": False, "non-lesional": False, "nonlesional": False,
    "n": False,
}

FILENAME_PATTERNS = [
    re.compile(
        r"^(GSM\d+)_"
        r"P(?P<patient_id>\d+)"
        r"_V(?P<visit_num>\d+)"
        r"_(?P<lesion>[NL])L?"
        r"(?:_(?P<extra>.+))?$"
    ),
    re.compile(
        r"^(GSM\d+)"
        r"(?:_(?P<field_1>[^_]+))?"
        r"(?:_(?P<field_2>[^_]+))?"
        r"(?:_(?P<field_3>[^_]+))?"
        r"(?:_(?P<field_4>[^_]+))?"
        r"(?:_(?P<field_5>[^_]+))?$"
    ),
]


def _interpret_fields(match_dict: dict) -> dict:
    result = {}
    raw = {}
    for k, v in match_dict.items():
        if v is not None:
            raw[k] = v

    if "patient_id" in raw:
        result["patient_id"] = f"P{raw['patient_id']}"

    if "visit_num" in raw:
        try:
            vnum = int(raw["visit_num"])
        except ValueError:
            vnum = raw["visit_num"]
        visit_label = f"V{vnum}"
        result["visit_label"] = visit_label
        result["visit_num"] = vnum
        real_time = _map_visit_to_time(visit_label)
        result["time_label"] = real_time

    if "lesion" in raw:
        l = raw["lesion"].upper()
        is_lesional = LESION_KEYWORDS.get(l.lower(), None)
        if is_lesional is not None:
            result["lesion_status"] = "Lesional" if is_lesional else "Non-lesional"
            result["is_lesional"] = is_lesional
        if "time_label" in result:
            if is_lesional:
                result["condition_raw"] = f"{result['time_label']} lesional skin"
            else:
                result["condition_raw"] = f"{result['time_label']} non-lesional skin"

    if "gsm" in raw:
        result["gsm"] = raw["gsm"]

    for i in range(1, 6):
        fk = f"field_{i}"
        if fk in raw and raw[fk]:
            result[f"meta_{i}"] = raw[fk]

    return result


def parse_filename_metadata(sample_name: str) -> dict:
    for pattern in FILENAME_PATTERNS:
        m = pattern.match(sample_name)
        if m:
            return _interpret_fields(m.groupdict())
    return {}


def enrich_from_filenames(sample_ids: list[str]) -> pd.DataFrame:
    rows = []
    for sid in sample_ids:
        meta = parse_filename_metadata(sid)
        meta["sample_id"] = sid
        rows.append(meta)
    df = pd.DataFrame(rows).set_index("sample_id")

    # Only add standard metadata columns that have at least one non-None value.
    # This prevents downstream code from seeing all-empty columns (e.g. when
    # sample IDs don't match GSM patterns).
    for col in [
        "patient_id", "time_label", "visit_label", "visit_num",
        "lesion_status", "is_lesional", "condition_raw", "gsm",
    ]:
        if col not in df.columns:
            continue
        non_null = df[col].notna() & (df[col].astype(str).str.strip() != "")
        if not non_null.any():
            df.drop(columns=[col], inplace=True)
    return df


def enrich_from_table(table_path: str, sample_ids: list[str]) -> pd.DataFrame:
    path = Path(table_path)
    try:
        if path.suffix in (".csv", ".tsv", ".txt"):
            df = pd.read_csv(table_path)
        elif path.suffix in (".xlsx", ".xls"):
            df = pd.read_excel(table_path)
        else:
            logger.warning("Unrecognized metadata table format: %s", table_path)
            return pd.DataFrame()
    except Exception as e:
        logger.warning("Failed to read metadata table %s: %s", table_path, e)
        return pd.DataFrame()

    id_col = None
    for col in ["sample_id", "sample", "GSM", "geo_accession", "Sample", "sample_name"]:
        if col in df.columns:
            id_col = col
            break
    if id_col is None:
        logger.warning("Metadata table has no sample_id column. Columns: %s", list(df.columns)[:10])
        return pd.DataFrame()

    df = df.rename(columns={id_col: "sample_id"})
    df["sample_id"] = df["sample_id"].astype(str)
    df = df.set_index("sample_id")
    matched = df[df.index.isin(sample_ids)]
    return matched


def enrich_adata_obs(
    adata: "ad.AnnData",
    sample_metadata_table: str = None,
) -> "ad.AnnData":
    if "sample_id" not in adata.obs.columns:
        logger.warning("No sample_id column in obs")
        return adata

    sample_ids = sorted(adata.obs["sample_id"].unique())
    meta_df = None

    if sample_metadata_table:
        meta_df = enrich_from_table(sample_metadata_table, sample_ids)

    if meta_df is None or meta_df.empty:
        meta_df = enrich_from_filenames(sample_ids)

    if meta_df.empty:
        logger.info("No metadata could be inferred for %d samples", len(sample_ids))
        return adata

    obs = adata.obs.copy()
    obs = obs.join(meta_df, on="sample_id", how="left")

    cols_written = []
    for col in meta_df.columns:
        vals = obs[col].values
        # Skip columns that are entirely empty after join
        str_vals = pd.Series(vals).astype(str)
        has_content = (str_vals.notna() & (str_vals.str.strip() != "") &
                       (str_vals.str.strip() != "nan") & (str_vals.str.strip() != "none"))
        if not has_content.any():
            logger.debug("Skipping all-empty metadata column: %s", col)
            continue
        cols_written.append(col)
        if col in ("gsm", "patient_id", "time_label", "visit_label", "condition_raw", "lesion_status"):
            adata.obs[col] = pd.array(vals, dtype="string").fillna("")
        elif col in ("is_lesional",):
            adata.obs[col] = vals.astype(float)
        elif col == "visit_num":
            adata.obs[col] = pd.to_numeric(vals, errors="coerce")
        else:
            adata.obs[col] = pd.Categorical(vals.astype(str))

    # Record metadata inference method
    adata.uns["metadata_inference_method"] = _mapping_source
    adata.uns["metadata_visit_mapping"] = str([
        (p.pattern, mp) for p, mp in _active_visit_mappings if mp
    ])
    logger.info(
        "Enriched obs with metadata: %d columns (%s) [source: %s]",
        len(cols_written), ", ".join(cols_written), _mapping_source,
    )
    return adata
