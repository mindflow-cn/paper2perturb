#!/usr/bin/env python3
"""Rule-based QC for Paper2Perturb h5ad/test_case/xlsx consistency."""

from __future__ import annotations

import argparse
import ast
import json
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import pandas as pd


DEFAULT_TABLE = Path("result.xlsx")
PERTURB_COLS = [
    "orig.ident",
    "dose",
    "dose_label",
    "condition",
    "time",
    "treatment_time",
    "sample",
    "sample_id",
    "Sample",
    "batch",
]
REQUIRED_OBS_COLS = ["cell_type"]
TOP_LEVEL_FIELDS = [
    "benchmark_id",
    "tissue",
    "source_type",
    "perturbation_type",
    "perturbation_name",
    "smiles",
    "description",
]
CASE_FIELDS = [
    "test_id",
    "benchmark_id",
    "perturb_var",
    "control",
    "dose_groups",
    "time_groups",
    "relation",
]
ALLOWED_DOSE_UNITS_TEXT = "mM, uM, µM, nM, mg, mg/day, mg/kg"
DOSE_RE = re.compile(
    r"^\s*~?\s*[-+]?\d+(?:\.\d+)?\s*(?:mM|uM|µM|nM|mg|mg\s*/\s*(?:day|kg))(?:\s+[A-Za-z][A-Za-z0-9_-]*)*\s*$",
    re.I,
)
TIME_RE = re.compile(
    r"^\s*~?\s*[-+]?\d+(?:\.\d+)?\s*(?:min|mins|minute|minutes|h|hr|hrs|hour|hours|day|days|week|weeks)\s*$",
    re.I,
)
NA_STRINGS = {"", "nan", "none", "na", "n/a", "null"}


@dataclass
class Check:
    status: str
    item: str
    detail: str
    suggestion: str = ""


def clean_scalar(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    return str(value).strip()


def parse_array(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, float) and math.isnan(value):
        return []
    if isinstance(value, list):
        return [clean_scalar(v) for v in value if clean_scalar(v).lower() not in NA_STRINGS]
    if isinstance(value, tuple):
        return [clean_scalar(v) for v in value if clean_scalar(v).lower() not in NA_STRINGS]
    text = clean_scalar(value)
    if text.lower() in NA_STRINGS:
        return []
    try:
        parsed = json.loads(text)
    except Exception:
        try:
            parsed = ast.literal_eval(text)
        except Exception:
            parsed = None
    if isinstance(parsed, list):
        return [clean_scalar(v) for v in parsed if clean_scalar(v).lower() not in NA_STRINGS]
    if isinstance(parsed, tuple):
        return [clean_scalar(v) for v in parsed if clean_scalar(v).lower() not in NA_STRINGS]
    return [text]


def unique_nonempty(values: Iterable[Any], limit: int | None = None) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        text = clean_scalar(value)
        if text.lower() in NA_STRINGS or text in seen:
            continue
        seen.add(text)
        out.append(text)
        if limit and len(out) >= limit:
            break
    return out


def normalize_text(value: Any) -> str:
    return re.sub(r"\s+", " ", clean_scalar(value))


def normalize_label(value: str) -> str:
    text = normalize_text(value).replace("µ", "u").lower()
    text = re.sub(r"\bhrs?\b|\bhours?\b", "h", text)
    text = re.sub(r"\bminutes?\b|\bmins?\b", "min", text)
    text = re.sub(r"\bdays\b", "day", text)
    text = re.sub(r"\bweeks\b", "week", text)
    return text


def comparable_array(value: Any) -> list[str]:
    return [normalize_text(v) for v in parse_array(value)]


def is_dose(value: str) -> bool:
    return bool(DOSE_RE.match(value))


def is_time(value: str) -> bool:
    return bool(TIME_RE.match(value))


def invalid_dose_labels(values: Iterable[str]) -> list[str]:
    return [value for value in values if not is_dose(value)]


def json_equal(a: Any, b: Any) -> bool:
    aa = comparable_array(a)
    bb = comparable_array(b)
    if len(aa) != 1 or len(bb) != 1:
        return aa == bb
    return aa[0] == bb[0]


def add(checks: list[Check], status: str, item: str, detail: str, suggestion: str = "") -> None:
    checks.append(Check(status, item, detail, suggestion))


def read_h5ad_obs(path: Path) -> tuple[list[str], dict[str, list[str]], str | None]:
    if not path.exists():
        return [], {}, f"missing file: {path}"
    try:
        import anndata as ad
    except Exception as exc:
        return [], {}, f"cannot import anndata: {exc}"
    try:
        adata = ad.read_h5ad(path, backed="r")
        cols = list(adata.obs.columns)
        value_cols = list(dict.fromkeys(PERTURB_COLS + REQUIRED_OBS_COLS))
        values = {
            col: unique_nonempty(adata.obs[col].astype(str).values)
            for col in value_cols
            if col in adata.obs.columns
        }
        adata.file.close()
        return cols, values, None
    except Exception as exc:
        return [], {}, f"cannot read {path}: {exc}"


def infer_benchmark_id(test_id: str) -> str:
    parts = test_id.split("_")
    if len(parts) < 2:
        raise ValueError(f"cannot infer benchmark_id from test_id: {test_id}")
    return "_".join(parts[:2])


def pick_sheet(path: Path, preferred: list[str]) -> str:
    xl = pd.ExcelFile(path)
    for name in preferred:
        if name in xl.sheet_names:
            return name
    raise ValueError(f"{path} does not contain any of these sheets: {preferred}; found={xl.sheet_names}")


def load_table(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    benchmark_sheet = pick_sheet(path, ["Benchmark-level", "Benchmark"])
    case_sheet = pick_sheet(path, ["Case-level", "Case"])
    benchmark = pd.read_excel(path, sheet_name=benchmark_sheet, dtype=object)
    case = pd.read_excel(path, sheet_name=case_sheet, dtype=object)
    return benchmark, case


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def find_json_case(data: dict[str, Any], test_id: str) -> dict[str, Any] | None:
    for case in data.get("test_cases", []):
        if case.get("test_id") == test_id:
            return case
    return None


def resolve_benchmark_dir(package_dir: Path, benchmark_id: str) -> Path:
    candidates = [
        package_dir / benchmark_id,
        package_dir / "data" / benchmark_id,
        package_dir / "h5ad" / benchmark_id,
        Path("data") / benchmark_id,
        Path("h5ad") / benchmark_id,
    ]
    for path in candidates:
        if (path / "test_case.json").exists() or (path / "ground_truth.h5ad").exists() or (path / "control.h5ad").exists():
            return path
    return candidates[0]


def compare_json_table(
    checks: list[Check],
    data: dict[str, Any],
    json_case: dict[str, Any] | None,
    bench_row: pd.Series | None,
    case_row: pd.Series | None,
    benchmark_id: str,
    test_id: str,
) -> None:
    if bench_row is None:
        add(checks, "FAIL", "Benchmark row", f"{benchmark_id} not found in Benchmark sheet", "Add or correct the Benchmark row in the metadata table.")
    else:
        add(checks, "PASS", "Benchmark row", f"{benchmark_id} found in Benchmark sheet")
        for field in TOP_LEVEL_FIELDS:
            expected = bench_row.get(field, "")
            actual = data.get(field, "")
            if json_equal(actual, expected):
                add(checks, "PASS", f"JSON top-level {field}", f"JSON={clean_scalar(actual)!r}, table={clean_scalar(expected)!r}")
            else:
                add(
                    checks,
                    "FAIL",
                    f"JSON top-level {field}",
                    f"JSON={clean_scalar(actual)!r}, table={clean_scalar(expected)!r}",
                    f"Update test_case.json top-level `{field}` to match Benchmark.{field}.",
                )

    if json_case is None:
        add(checks, "FAIL", "JSON case", f"{test_id} not found in test_case.json", "Add this case from the Case sheet or correct the test_id.")
        return
    add(checks, "PASS", "JSON case", f"{test_id} found in test_case.json")

    if case_row is None:
        add(checks, "FAIL", "Case row", f"{test_id} not found in Case sheet", "Add or correct the Case row in the metadata table.")
        return
    add(checks, "PASS", "Case row", f"{test_id} found in Case sheet")

    for field in CASE_FIELDS:
        expected = case_row.get(field, "")
        actual = data.get("benchmark_id", "") if field == "benchmark_id" and field not in json_case else json_case.get(field, "")
        if json_equal(actual, expected):
            add(checks, "PASS", f"JSON case {field}", f"JSON={parse_array(actual) if field.endswith('_groups') else clean_scalar(actual)!r}, table={parse_array(expected) if field.endswith('_groups') else clean_scalar(expected)!r}")
        else:
            add(
                checks,
                "FAIL",
                f"JSON case {field}",
                f"JSON={parse_array(actual) if field.endswith('_groups') else clean_scalar(actual)!r}, table={parse_array(expected) if field.endswith('_groups') else clean_scalar(expected)!r}",
                f"Update test_case.json case `{test_id}` field `{field}` from the Case sheet.",
            )

    table_drug = clean_scalar(case_row.get("drug", ""))
    json_drug = clean_scalar(json_case.get("drug", data.get("perturbation_name", "")))
    if table_drug and normalize_text(table_drug) != normalize_text(json_drug):
        add(checks, "FAIL", "drug / perturbation_name", f"JSON/case drug={json_drug!r}, Case.drug={table_drug!r}", "Align case drug and top-level perturbation_name with the Case sheet.")
    elif table_drug:
        add(checks, "PASS", "drug / perturbation_name", f"{json_drug!r}")


def check_perturb_var(checks: list[Check], source_case: dict[str, Any], test_id: str) -> tuple[str, list[str]]:
    perturb_var = clean_scalar(source_case.get("perturb_var", ""))
    expected = parse_array(source_case.get("dose_groups" if perturb_var == "dose" else "time_groups", []))
    if perturb_var not in {"dose", "time"}:
        add(
            checks,
            "FAIL",
            "perturb_var legality",
            f"perturb_var={perturb_var!r}",
            "Change perturb_var to `dose` or `time`, or explicitly extend the schema/checker for this variable type.",
        )
        return perturb_var, expected
    add(checks, "PASS", "perturb_var legality", f"perturb_var={perturb_var!r}")
    if not expected:
        group_field = "dose_groups" if perturb_var == "dose" else "time_groups"
        add(
            checks,
            "FAIL",
            f"{group_field} completeness",
            f"{group_field} is empty for {test_id}",
            f"Fill `{group_field}` in the metadata table and test_case.json before validating h5ad.",
        )
    else:
        add(checks, "PASS", "expected perturbation groups", f"{expected}")
        if perturb_var == "dose":
            invalid = invalid_dose_labels(expected)
            if invalid:
                add(
                    checks,
                    "FAIL",
                    "dose_groups unit",
                    f"invalid dose_groups labels: {invalid}",
                    f"Use dose labels with units from this allowed set: {ALLOWED_DOSE_UNITS_TEXT}.",
                )
            else:
                add(checks, "PASS", "dose_groups unit", f"all dose_groups use allowed units: {ALLOWED_DOSE_UNITS_TEXT}")
    return perturb_var, expected


def classify_sets(observed: list[str], expected: list[str]) -> tuple[str, str]:
    o = set(observed)
    e = set(expected)
    if not o:
        return "FAIL", "orig.ident has no usable values"
    if o == e:
        return "PASS", "O == E"
    if e and e < o:
        return "FAIL", f"E < O; extra={sorted(o - e)}"
    if o and o < e:
        return "FAIL", f"O < E; missing={sorted(e - o)}"
    if o & e:
        return "FAIL", f"partial overlap; extra={sorted(o - e)}, missing={sorted(e - o)}"
    return "FAIL", "no overlap"


def check_required_cell_type(checks: list[Check], label: str, cols: list[str], obs_values: dict[str, list[str]]) -> None:
    if "cell_type" not in cols:
        add(
            checks,
            "FAIL",
            f"{label} cell_type column",
            "`cell_type` missing",
            "Add `obs['cell_type']` and populate it with the standardized cell type for every cell.",
        )
        return
    values = obs_values.get("cell_type", [])
    if not values:
        add(
            checks,
            "FAIL",
            f"{label} cell_type values",
            "`cell_type` is empty",
            "Populate `obs['cell_type']` with non-empty standardized cell type values.",
        )
        return
    add(checks, "PASS", f"{label} cell_type values", f"{values}")


def check_control_orig_ident(checks: list[Check], cols: list[str], obs_values: dict[str, list[str]]) -> None:
    """Validate control orig.ident against the canonical control h5ad label.

    The Case/control xlsx value can carry source-specific control metadata.
    It is intentionally not used for control.h5ad.obs['orig.ident'], which is
    normalized to the fixed label "C" across converted benchmark files.
    """
    expected = ["C"]
    if "orig.ident" not in cols:
        add(
            checks,
            "FAIL",
            "control.h5ad orig.ident column",
            "`orig.ident` missing",
            "Add `control.h5ad.obs['orig.ident']` and set every value to 'C'.",
        )
        return

    values = obs_values.get("orig.ident", [])
    if not values:
        add(
            checks,
            "FAIL",
            "control.h5ad orig.ident values",
            "`orig.ident` is empty",
            "Populate `control.h5ad.obs['orig.ident']` with 'C'.",
        )
        return

    status, detail = classify_sets(values, expected)
    suggestion = ""
    if status != "PASS":
        suggestion = "Set or remap `control.h5ad.obs['orig.ident']` to the fixed control label 'C'."
    add(checks, status, "control.h5ad orig.ident fixed C label", f"O={values}; E={expected}; {detail}", suggestion)


def check_h5ad(
    checks: list[Check],
    bench_dir: Path,
    benchmark_id: str,
    perturb_var: str,
    expected: list[str],
) -> None:
    for filename in ["test_case.json", "ground_truth.h5ad", "control.h5ad"]:
        path = bench_dir / filename
        if path.exists():
            add(checks, "PASS", f"required file {filename}", str(path))
        else:
            add(checks, "FAIL", f"required file {filename}", str(path), f"Create or restore `{filename}` in `{benchmark_id}`.")

    gt_path = bench_dir / "ground_truth.h5ad"
    cols, obs_values, error = read_h5ad_obs(gt_path)
    if error:
        add(checks, "FAIL", "read ground_truth.h5ad", error, "Make sure anndata can read the h5ad and the file is not corrupted.")
        return
    add(checks, "PASS", "read ground_truth.h5ad", f"{gt_path}; obs columns={len(cols)}")
    check_required_cell_type(checks, "ground_truth.h5ad", cols, obs_values)

    control_path = bench_dir / "control.h5ad"
    control_cols, control_obs_values, control_error = read_h5ad_obs(control_path)
    if control_error:
        add(checks, "FAIL", "read control.h5ad", control_error, "Make sure anndata can read the control h5ad and the file is not corrupted.")
    else:
        add(checks, "PASS", "read control.h5ad", f"{control_path}; obs columns={len(control_cols)}")
        check_required_cell_type(checks, "control.h5ad", control_cols, control_obs_values)
        check_control_orig_ident(checks, control_cols, control_obs_values)

    if "orig.ident" not in cols:
        add(
            checks,
            "FAIL",
            "orig.ident column",
            f"`orig.ident` missing; available perturbation-like columns={sorted(obs_values)}",
            "Add `obs['orig.ident']` and populate it with the case's expected dose_groups/time_groups values.",
        )
        return

    orig = obs_values.get("orig.ident", [])
    if not orig:
        add(checks, "FAIL", "orig.ident values", "`orig.ident` is empty", "Populate `orig.ident` with expected perturbation group labels.")
    else:
        add(checks, "PASS", "orig.ident values", f"{orig}")

    if expected:
        status, detail = classify_sets(orig, expected)
        suggestion = ""
        if status != "PASS":
            group_field = "dose_groups" if perturb_var == "dose" else "time_groups"
            suggestion = f"Set or remap `ground_truth.h5ad.obs['orig.ident']` to the current case `{group_field}` values: {expected}."
        add(checks, status, "orig.ident vs expected groups", f"O={orig}; E={expected}; {detail}", suggestion)

    dose_like = {col: vals for col, vals in obs_values.items() if col in {"orig.ident", "dose", "dose_label", "condition"}}
    time_like = {col: vals for col, vals in obs_values.items() if col in {"orig.ident", "time", "treatment_time", "condition"}}
    active = dose_like if perturb_var == "dose" else time_like
    expected_set = set(expected)
    for col, vals in active.items():
        val_set = set(vals)
        if not vals:
            continue
        if expected_set and val_set <= expected_set:
            add(checks, "PASS", f"h5ad obs.{col} maps to {perturb_var}", f"{vals}")
        elif col == "orig.ident":
            continue
        else:
            add(
                checks,
                "WARN",
                f"h5ad obs.{col} maps to {perturb_var}",
                f"values={vals}; expected={expected}",
                f"If `{col}` is intended to encode {perturb_var}, replace/remap it to {expected}; otherwise document it as non-perturbation metadata.",
            )

    check_type_units(checks, perturb_var, orig, expected)


def check_type_units(checks: list[Check], perturb_var: str, orig: list[str], expected: list[str]) -> None:
    values = orig or []
    if perturb_var == "dose":
        bad_type = [v for v in values if is_time(v)]
        missing_unit = [v for v in values if re.fullmatch(r"\s*[-+]?\d+(?:\.\d+)?\s*", v)]
        non_dose = [v for v in values if not is_dose(v)]
        if bad_type:
            add(checks, "FAIL", "dose/time type match", f"orig.ident contains time values: {bad_type}", "Replace time labels with dose labels from dose_groups.")
        if missing_unit:
            add(checks, "FAIL", "dose unit", f"naked numeric dose labels: {missing_unit}", f"Add explicit dose units from the allowed set: {ALLOWED_DOSE_UNITS_TEXT}.")
        if non_dose and not bad_type and not missing_unit:
            add(checks, "FAIL", "dose label format", f"invalid dose labels: {non_dose}", f"Use dose labels with units from this allowed set: {ALLOWED_DOSE_UNITS_TEXT}.")
        if not bad_type and not missing_unit and not non_dose:
            add(checks, "PASS", "dose label format", f"orig.ident dose labels use allowed units: {ALLOWED_DOSE_UNITS_TEXT}")
    elif perturb_var == "time":
        bad_type = [v for v in values if is_dose(v)]
        missing_unit = [v for v in values if re.fullmatch(r"\s*[-+]?\d+(?:\.\d+)?\s*", v)]
        non_time = [v for v in values if not is_time(v)]
        if bad_type:
            add(checks, "FAIL", "dose/time type match", f"orig.ident contains dose values: {bad_type}", "Replace dose labels with time labels from time_groups.")
        if missing_unit:
            add(checks, "FAIL", "time unit", f"naked numeric time labels: {missing_unit}", "Add explicit time units such as hour/day from the source table.")
        if non_time and not bad_type and not missing_unit:
            add(checks, "WARN", "time label format", f"non-standard time labels: {non_time}", f"Use explicit time labels matching time_groups: {expected}.")
        if not bad_type and not missing_unit and not non_time:
            add(checks, "PASS", "time label format", "orig.ident values look like time labels")


def choose_source_case(json_case: dict[str, Any] | None, case_row: pd.Series | None) -> dict[str, Any]:
    if case_row is None:
        return json_case or {}
    return {
        "perturb_var": clean_scalar(case_row.get("perturb_var", "")),
        "dose_groups": parse_array(case_row.get("dose_groups", [])),
        "time_groups": parse_array(case_row.get("time_groups", [])),
    }


def run(args: argparse.Namespace) -> tuple[list[Check], dict[str, str]]:
    table = Path(args.table)
    package_dir = Path(args.package_dir) if args.package_dir else table.parent
    test_id = args.test_id
    benchmark_id = args.benchmark_id or infer_benchmark_id(test_id)
    checks: list[Check] = []
    bench_dir = resolve_benchmark_dir(package_dir, benchmark_id)
    meta = {
        "test_id": test_id,
        "benchmark_id": benchmark_id,
        "table": str(table),
        "package_dir": str(package_dir),
        "benchmark_dir": str(bench_dir),
    }

    if not table.exists():
        add(checks, "FAIL", "table path", f"{table} does not exist", "Pass --table with the correct result.xlsx/cell_line.xlsx path.")
        return checks, meta

    try:
        benchmark_df, case_df = load_table(table)
    except Exception as exc:
        add(checks, "FAIL", "read metadata table", str(exc), "Use result.xlsx with Benchmark-level/Case-level sheets or cell_line.xlsx with Benchmark/Case sheets.")
        return checks, meta
    bench_matches = benchmark_df[benchmark_df["benchmark_id"].astype(str) == benchmark_id]
    case_matches = case_df[case_df["test_id"].astype(str) == test_id]
    bench_row = None if bench_matches.empty else bench_matches.iloc[0]
    case_row = None if case_matches.empty else case_matches.iloc[0]

    json_path = bench_dir / "test_case.json"
    if json_path.exists():
        data = load_json(json_path)
        json_case = find_json_case(data, test_id)
        add(checks, "PASS", "read test_case.json", str(json_path))
    else:
        data = {}
        json_case = None
        add(checks, "FAIL", "read test_case.json", f"{json_path} missing", "Restore test_case.json from the table metadata.")

    compare_json_table(checks, data, json_case, bench_row, case_row, benchmark_id, test_id)
    source_case = choose_source_case(json_case, case_row)
    perturb_var, expected = check_perturb_var(checks, source_case, test_id)
    check_h5ad(checks, bench_dir, benchmark_id, perturb_var, expected)
    return checks, meta


def render_markdown(checks: list[Check], meta: dict[str, str]) -> str:
    fail_count = sum(c.status == "FAIL" for c in checks)
    warn_count = sum(c.status == "WARN" for c in checks)
    pass_count = sum(c.status == "PASS" for c in checks)
    lines = [
        f"# h5ad orig.ident QC: `{meta['test_id']}`",
        "",
        f"- benchmark_id: `{meta['benchmark_id']}`",
        f"- table: `{meta['table']}`",
        f"- package_dir: `{meta['package_dir']}`",
        f"- benchmark_dir: `{meta['benchmark_dir']}`",
        f"- summary: PASS={pass_count}, WARN={warn_count}, FAIL={fail_count}",
        "",
        "| status | item | detail | suggestion |",
        "|--------|------|--------|------------|",
    ]
    for c in checks:
        lines.append(f"| {c.status} | {escape(c.item)} | {escape(c.detail)} | {escape(c.suggestion)} |")
    return "\n".join(lines)


def render_json(checks: list[Check], meta: dict[str, str]) -> str:
    return json.dumps(
        {
            "meta": meta,
            "summary": {
                "PASS": sum(c.status == "PASS" for c in checks),
                "WARN": sum(c.status == "WARN" for c in checks),
                "FAIL": sum(c.status == "FAIL" for c in checks),
            },
            "checks": [c.__dict__ for c in checks],
        },
        ensure_ascii=False,
        indent=2,
    )


def escape(value: str) -> str:
    return clean_scalar(value).replace("|", "\\|").replace("\n", "<br>")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--test-id", required=True, help="Case test_id, e.g. 32094658_01_01")
    parser.add_argument("--table", default=str(DEFAULT_TABLE), help="Path to result.xlsx or cell_line.xlsx")
    parser.add_argument("--package-dir", default="", help="Root directory for benchmark data; defaults to table parent and searches data/<benchmark_id>, h5ad/<benchmark_id>, and <benchmark_id>")
    parser.add_argument("--benchmark-id", default="", help="Override inferred benchmark_id")
    parser.add_argument("--format", choices=["markdown", "json"], default="markdown")
    args = parser.parse_args(argv)

    checks, meta = run(args)
    print(render_json(checks, meta) if args.format == "json" else render_markdown(checks, meta))
    return 1 if any(c.status == "FAIL" for c in checks) else 0


if __name__ == "__main__":
    sys.exit(main())
