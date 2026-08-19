"""Parse Case-level table to extract target cell types and determine strategy category.

Reads a Case-level CSV (or a metadata xlsx with a Case-level sheet).
For the requested benchmark_id, collects unique non-empty cell_type values
as target_cell_types, then chooses a strategy_category.

FAIL-FAST: missing benchmark_id, missing cell_type column, or empty cell_type
values all raise ValueError immediately.
"""

import ast
import logging
import re
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

# Cell-type keywords for strategy categorization
CELL_LINE_KEYWORDS = [
    "cell line", "organoid", "hesc", "ipsc", "induced pluripotent",
    "embryonic stem", "progenitor", "cultured", "crispr",
    "primary culture", "primary_culture",
    "hek293", "hela", "k562", "jurkat", "lncap", "mcf7", "a549",
    "ht29", "hct116", "hacat", "sh-sy5y", "u2os", "cho", "nih3t3",
    "raw264", "thp-1", "hl-60", "raji",
]

TUMOR_KEYWORDS = [
    "tumor", "cancer", "malignant", "carcinoma", "melanoma",
    "leukemia", "lymphoma", "metastasis", "tme", "caf",
    "cancer-associated fibroblast", "glioblastoma", "sarcoma",
    "adenocarcinoma", "neoplastic",
]

IMMUNE_KEYWORDS = [
    "cd4", "cd8", "t cell", "b cell", "nk cell", "natural killer",
    "plasma cell", "monocyte", "macrophage", "dendritic cell",
    "neutrophil", "pbmc", "lymphocyte", "myeloid",
    "treg", "th17", "th1", "th2", "ctl", "tcm", "tem", "temra",
    "mast cell", "basophil", "eosinophil", "ilc",
]

TISSUE_KEYWORDS = [
    "keratinocyte", "fibroblast", "endothelial", "epithelial",
    "neuronal", "glial", "hepatocyte", "pericyte",
    "smooth muscle", "mesenchymal", "chondrocyte", "osteoblast",
    "melanocyte", "adipocyte", "myocyte", "alveolar",
    "pneumocyte", "podocyte", "acinar", "ductal",
    "goblet cell", "ciliated", "enterocyte", "paneth",
]


def _normalize(s: str) -> str:
    """Lowercase, replace underscores/hyphens with spaces, collapse whitespace."""
    t = str(s).strip().lower()
    t = t.replace("_", " ").replace("-", " ")
    t = re.sub(r"\s+", " ", t)
    return t


def read_case_table(case_table_path: str) -> pd.DataFrame:
    """Read the Case-level table from CSV or xlsx."""
    path = Path(case_table_path)
    if path.suffix.lower() in (".xlsx", ".xls"):
        df = pd.read_excel(path, sheet_name="Case-level")
    else:
        try:
            df = pd.read_csv(path)
        except (UnicodeDecodeError, pd.errors.ParserError):
            try:
                df = pd.read_csv(path, encoding="utf-8-sig")
            except (UnicodeDecodeError, pd.errors.ParserError):
                df = pd.read_csv(path, encoding="latin-1")
    return df


def _find_cell_type_column(df: pd.DataFrame) -> str:
    """Find the cell_type column by name. Raises if not found."""
    candidates = [
        "cell_type", "celltype", "cell type",
        "cell_type_or_sample_system", "sample_system",
    ]
    for col in candidates:
        if col in df.columns:
            return col
    for col in df.columns:
        if "cell_type" in col.lower() or "celltype" in col.lower():
            return col
    raise ValueError(
        f"Cannot find cell_type column in Case-level table. "
        f"Available columns: {list(df.columns)}"
    )


def _find_benchmark_id_column(df: pd.DataFrame) -> str:
    """Find the benchmark_id column by name. Raises if not found."""
    candidates = ["benchmark_id", "benchmark id", "benchmark_id_sanitized"]
    for col in candidates:
        if col in df.columns:
            return col
    raise ValueError(
        f"Cannot find benchmark_id column in Case-level table. "
        f"Available columns: {list(df.columns)}"
    )


def classify_strategy(target_cell_types: list[str]) -> str:
    """Choose strategy_category from target_cell_types.

    Precedence: cell_line_or_in_vitro > tumor_or_tme >
                normal_or_disease_tissue > immune_enriched > unknown
    """
    all_text = " ".join(_normalize(t) for t in target_cell_types)

    # Check for specific immune cell type markers first — these indicate
    # a heterogeneous mixture, not a homogeneous cell line / primary culture.
    _specific_immune = {"cd4", "cd8", "b cell", "nk cell", "t cell", "plasma cell",
                        "monocyte", "macrophage", "dendritic cell", "neutrophil"}

    for kw in CELL_LINE_KEYWORDS:
        if kw in all_text:
            # "primary culture" / "primary cells" alone shouldn't override
            # specific immune cell type markers — the sample is heterogeneous.
            if kw in ("primary culture", "primary_culture"):
                if any(imm in all_text for imm in _specific_immune):
                    break  # fall through to immune/tissue check
            return "cell_line_or_in_vitro"

    for kw in TUMOR_KEYWORDS:
        if kw in all_text:
            return "tumor_or_tme"

    has_immune = any(kw in all_text for kw in IMMUNE_KEYWORDS)
    has_tissue = any(kw in all_text for kw in TISSUE_KEYWORDS)

    if has_tissue:
        return "normal_or_disease_tissue"

    if has_immune:
        all_immune = True
        for t in target_cell_types:
            t_norm = _normalize(t)
            if not any(kw in t_norm for kw in IMMUNE_KEYWORDS):
                all_immune = False
                break
        if all_immune:
            return "immune_enriched"
        else:
            return "normal_or_disease_tissue"

    return "unknown"


def get_target_cell_types(
    case_table_path: str,
    benchmark_id: str,
) -> dict:
    """Extract target cell types and strategy for a benchmark.

    FAILS FAST on:
      - benchmark_id not found in the table
      - cell_type column not found
      - no non-empty cell_type values for the benchmark

    Returns:
        dict with keys:
            target_cell_types: list[str]
            strategy_category: str
            target_cell_type_source: str
            strategy_reason: str
            benchmark_id: str
    """
    df = read_case_table(case_table_path)

    # Must find columns first
    bid_col = _find_benchmark_id_column(df)
    ct_col = _find_cell_type_column(df)

    available_bids = sorted(df[bid_col].dropna().astype(str).str.strip().unique())

    # Filter to this benchmark
    mask = df[bid_col].astype(str).str.strip() == str(benchmark_id).strip()
    subset = df[mask]

    if subset.empty:
        raise ValueError(
            f"benchmark_id '{benchmark_id}' not found in Case-level table "
            f"(column: '{bid_col}').\n"
            f"Available benchmark_ids (showing up to 20): {available_bids[:20]}"
        )

    # Collect unique non-empty cell_type values
    raw_types = subset[ct_col].dropna().astype(str).str.strip()

    # Handle JSON-list format
    parsed = []
    for val in raw_types:
        if val.startswith("[") and val.endswith("]"):
            try:
                items = ast.literal_eval(val)
                if isinstance(items, list):
                    parsed.extend(items)
                    continue
            except (ValueError, SyntaxError):
                pass
        parsed.append(val)

    # Filter empty/sentinel
    clean = [t for t in parsed if t and t.lower() not in ("nan", "none", "")]
    if not clean:
        raise ValueError(
            f"benchmark_id '{benchmark_id}' has {len(subset)} rows but "
            f"no non-empty cell_type values in column '{ct_col}'. "
            f"Raw values: {raw_types.tolist()[:10]}"
        )

    # Unique, preserving original casing
    unique_norm = sorted(set(_normalize(t) for t in clean))
    orig_map = {}
    for t in clean:
        n = _normalize(t)
        if n not in orig_map and n in unique_norm:
            orig_map[n] = str(t).strip()

    target_cell_types = [orig_map.get(t, t) for t in unique_norm]
    strategy = classify_strategy(target_cell_types)

    logger.info(
        "benchmark_id=%s: target_cell_types=%s, strategy_category=%s",
        benchmark_id, target_cell_types, strategy,
    )

    return {
        "target_cell_types": target_cell_types,
        "strategy_category": strategy,
        "target_cell_type_source": "case_table.cell_type",
        "strategy_reason": f"derived from {len(clean)} non-empty cell_type entries",
        "benchmark_id": benchmark_id,
    }
