#!/usr/bin/env python3
"""Rule-based CSV/XLSX field validators invoked by validate-metadata skill."""

import re
import json
import ast
import argparse
from pathlib import Path

import pandas as pd

VALID_SOURCE_TYPES = {"cell_line", "primary_culture", "patient_sample", "organoid", "PDX"}
VALID_PERTURB_VARS = {"dose", "time"}
VALID_RELATIONS = {"UP", "DOWN"}


def parse_list(val):
    """Parse a value that may be a JSON/Python list string or an actual list."""
    if isinstance(val, list):
        return val
    if pd.isna(val) or str(val).strip() == "":
        return []
    s = str(val).strip()
    for parser in (json.loads, ast.literal_eval):
        try:
            parsed = parser(s)
            if isinstance(parsed, list):
                return parsed
        except (json.JSONDecodeError, ValueError, SyntaxError):
            continue
    # Not parseable as list, return as-is
    return s


def check_benchmark_id(value, pmid, row_idx):
    """Rule 1: benchmark_id = <pmid>_<2-digit zero-padded index starting from 01>"""
    if pd.isna(value) or not isinstance(value, str):
        return f"Row {row_idx}: benchmark_id is missing or not a string (value={value})"

    m = re.match(r"^(\d+)_(\d{2,})$", str(value))
    if not m:
        return f"Row {row_idx}: benchmark_id '{value}' does not match pattern <pmid>_<NN>"

    actual_pmid, index_part = m.group(1), m.group(2)
    if actual_pmid != str(pmid):
        return f"Row {row_idx}: benchmark_id '{value}' pmid part '{actual_pmid}' != row pmid '{pmid}'"

    if len(index_part) < 2:
        return f"Row {row_idx}: benchmark_id '{value}' index '{index_part}' is not zero-padded to at least 2 digits"

    return None


def check_source_type(value, row_idx):
    """Rule 2: source_type must be a valid enum value."""
    if pd.isna(value) or str(value).strip() == "":
        return f"Row {row_idx}: source_type is empty"
    if str(value).strip() not in VALID_SOURCE_TYPES:
        return f"Row {row_idx}: source_type '{value}' not in {sorted(VALID_SOURCE_TYPES)}"
    return None


def check_smiles(value, row_idx, extra_columns):
    """Rule 3: smiles cannot be empty; if empty, must have a reason."""
    if pd.isna(value) or str(value).strip() == "":
        return f"Row {row_idx}: smiles is empty — must provide a reason or fill the field"
    return None


def check_test_id(value, benchmark_id, row_idx):
    """Rule 4: test_id = <benchmark_id>_<index>"""
    if pd.isna(value) or not isinstance(value, str):
        return f"Row {row_idx}: test_id is missing or not a string (value={value})"

    bid = str(benchmark_id) if not pd.isna(benchmark_id) else ""
    if not str(value).startswith(bid + "_"):
        return f"Row {row_idx}: test_id '{value}' does not start with benchmark_id '{bid}_'"
    return None


def check_dose_groups(value, row_idx):
    """Rule 5: dose_groups must be non-empty list; each entry must have mM, uM, or nM unit."""
    parsed = parse_list(value)
    if not isinstance(parsed, list) or len(parsed) == 0:
        return f"Row {row_idx}: dose_groups is empty or not a list (value={value})"

    for i, entry in enumerate(parsed):
        entry_str = str(entry).strip()
        if not re.search(r"\d+\.?\d*\s*(mM|uM|nM|mg|mg/kg|mg/day)", entry_str):
            return f"Row {row_idx}: dose_groups[{i}] '{entry_str}' missing mM/uM/nM/mg/mg/kg/mg/day unit"
    return None


def check_time_groups(value, perturb_var, row_idx):
    """Rule 6: if perturb_var == 'time', time_groups must be non-empty; each entry must have 'day' unit."""
    pv = str(perturb_var).strip() if not pd.isna(perturb_var) else ""
    if pv != "time":
        return None

    parsed = parse_list(value)
    if not isinstance(parsed, list) or len(parsed) == 0:
        return f"Row {row_idx}: perturb_var='time' but time_groups is empty or not a list (value={value})"

    for i, entry in enumerate(parsed):
        entry_str = str(entry).strip()
        if "day" not in entry_str.lower():
            return f"Row {row_idx}: time_groups[{i}] '{entry_str}' missing 'day' unit"
    return None


def check_perturb_var(value, row_idx):
    """Rule 7: perturb_var must be 'dose' or 'time'."""
    if pd.isna(value) or str(value).strip() == "":
        return f"Row {row_idx}: perturb_var is empty"
    if str(value).strip() not in VALID_PERTURB_VARS:
        return f"Row {row_idx}: perturb_var '{value}' not in {sorted(VALID_PERTURB_VARS)}"
    return None


def check_relation(value, row_idx):
    """Rule 8: relation must be 'UP' or 'DOWN'."""
    if pd.isna(value) or str(value).strip() == "":
        return f"Row {row_idx}: relation is empty"
    if str(value).strip() not in VALID_RELATIONS:
        return f"Row {row_idx}: relation '{value}' not in {sorted(VALID_RELATIONS)}"
    return None


def validate_benchmark_level(df):
    """Run benchmark-level checks (rules 1-3)."""
    errors = []
    for idx, row in df.iterrows():
        e = check_benchmark_id(row.get("benchmark_id"), row.get("pmid"), idx)
        if e:
            errors.append(e)
        e = check_source_type(row.get("source_type"), idx)
        if e:
            errors.append(e)
        e = check_smiles(row.get("smiles"), idx, None)
        if e:
            errors.append(e)
    return errors


def validate_case_level(df):
    """Run case-level checks (rules 4-8)."""
    errors = []
    for idx, row in df.iterrows():
        e = check_test_id(row.get("test_id"), row.get("benchmark_id"), idx)
        if e:
            errors.append(e)
        e = check_dose_groups(row.get("dose_groups"), idx)
        if e:
            errors.append(e)
        e = check_time_groups(row.get("time_groups"), row.get("perturb_var"), idx)
        if e:
            errors.append(e)
        e = check_perturb_var(row.get("perturb_var"), idx)
        if e:
            errors.append(e)
        e = check_relation(row.get("relation"), idx)
        if e:
            errors.append(e)
    return errors


def main():
    parser = argparse.ArgumentParser(description="Rule-based CSV/XLSX checks for validate-metadata skill")
    parser.add_argument("path", nargs="?", default="result.xlsx",
                        help="Path to the Excel file (default: result.xlsx)")
    parser.add_argument("--sheet-benchmark", default="Benchmark-level",
                        help="Benchmark-level sheet name")
    parser.add_argument("--sheet-case", default="Case-level",
                        help="Case-level sheet name")
    args = parser.parse_args()

    path = Path(args.path)
    if not path.exists():
        print(f"ERROR: {path} not found")
        return 1

    all_errors = []

    xl = pd.ExcelFile(path)

    if args.sheet_benchmark in xl.sheet_names:
        df_b = pd.read_excel(path, sheet_name=args.sheet_benchmark)
        all_errors += validate_benchmark_level(df_b)
        print(f"[Benchmark-level] {len(df_b)} rows checked, {len(all_errors)} errors found "
              f"(after checking this sheet)")
    else:
        print(f"[Benchmark-level] sheet '{args.sheet_benchmark}' not found, skipping")

    case_err_count_before = len(all_errors)

    if args.sheet_case in xl.sheet_names:
        df_c = pd.read_excel(path, sheet_name=args.sheet_case)
        all_errors += validate_case_level(df_c)
        case_errs = len(all_errors) - case_err_count_before
        print(f"[Case-level] {len(df_c)} rows checked, {case_errs} errors found")
    else:
        print(f"[Case-level] sheet '{args.sheet_case}' not found, skipping")

    if all_errors:
        print(f"\n=== {len(all_errors)} total violation(s) ===\n")
        for err in all_errors:
            print(f"  - {err}")
        return 1
    else:
        print("\nAll rule-based checks passed.")
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
