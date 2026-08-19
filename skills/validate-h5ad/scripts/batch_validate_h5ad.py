#!/usr/bin/env python3
"""Run validate_h5ad.py over every Case-level test_id row."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import pandas as pd


DEFAULT_TABLE = Path("result.xlsx")
DEFAULT_OUT_DIR = Path("tmp/validate_h5ad_batch")
CHECKER_PATH = Path(__file__).with_name("validate_h5ad.py")


def load_checker() -> Any:
    spec = importlib.util.spec_from_file_location("validate_h5ad", CHECKER_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import checker from {CHECKER_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def read_test_ids(table: Path) -> list[str]:
    xl = pd.ExcelFile(table)
    if "Case-level" in xl.sheet_names:
        sheet_name = "Case-level"
    elif "Case" in xl.sheet_names:
        sheet_name = "Case"
    else:
        raise ValueError(f"{table} does not contain Case-level or Case sheet; found={xl.sheet_names}")
    case_df = pd.read_excel(table, sheet_name=sheet_name, dtype=object)
    if "test_id" not in case_df.columns:
        raise ValueError(f"{table} {sheet_name} sheet does not contain a test_id column")
    return [
        str(value).strip()
        for value in case_df["test_id"].dropna().tolist()
        if str(value).strip()
    ]


def summarize_case(checker: Any, table: Path, package_dir: Path | None, test_id: str) -> dict[str, Any]:
    args = argparse.Namespace(
        test_id=test_id,
        table=str(table),
        package_dir=str(package_dir) if package_dir else "",
        benchmark_id="",
        format="json",
    )
    checks, meta = checker.run(args)
    status_counts = Counter(check.status for check in checks)
    failures = [check for check in checks if check.status == "FAIL"]
    warnings = [check for check in checks if check.status == "WARN"]
    overall = "FAIL" if failures else "WARN" if warnings else "PASS"
    return {
        "test_id": test_id,
        "benchmark_id": meta["benchmark_id"],
        "overall": overall,
        "PASS": status_counts["PASS"],
        "WARN": status_counts["WARN"],
        "FAIL": status_counts["FAIL"],
        "fail_items": "; ".join(check.item for check in failures),
        "warn_items": "; ".join(check.item for check in warnings),
        "fail_details": " || ".join(f"{check.item}: {check.detail}" for check in failures),
        "warn_details": " || ".join(f"{check.item}: {check.detail}" for check in warnings),
        "fail_suggestions": " || ".join(
            f"{check.item}: {check.suggestion}" for check in failures if check.suggestion
        ),
        "checks": [check.__dict__ for check in checks],
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fieldnames = [
        "test_id",
        "benchmark_id",
        "overall",
        "PASS",
        "WARN",
        "FAIL",
        "fail_items",
        "warn_items",
        "fail_details",
        "warn_details",
        "fail_suggestions",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def render_markdown(rows: list[dict[str, Any]], table: Path, package_dir: Path) -> str:
    overall_counts = Counter(row["overall"] for row in rows)
    fail_rows = [row for row in rows if row["FAIL"]]
    warn_rows = [row for row in rows if not row["FAIL"] and row["WARN"]]
    fail_by_item: dict[str, int] = defaultdict(int)
    warn_by_item: dict[str, int] = defaultdict(int)
    for row in rows:
        for check in row["checks"]:
            if check["status"] == "FAIL":
                fail_by_item[check["item"]] += 1
            elif check["status"] == "WARN":
                warn_by_item[check["item"]] += 1

    lines = [
        "# Batch h5ad QC Summary",
        "",
        f"- table: `{table}`",
        f"- package_dir: `{package_dir}`",
        f"- total test_id: {len(rows)}",
        f"- overall: PASS={overall_counts['PASS']}, WARN={overall_counts['WARN']}, FAIL={overall_counts['FAIL']}",
        "",
        "## Failure Item Counts",
        "",
    ]
    if fail_by_item:
        for item, count in sorted(fail_by_item.items(), key=lambda x: (-x[1], x[0])):
            lines.append(f"- {item}: {count}")
    else:
        lines.append("- none")

    lines.extend(["", "## Warning Item Counts", ""])
    if warn_by_item:
        for item, count in sorted(warn_by_item.items(), key=lambda x: (-x[1], x[0])):
            lines.append(f"- {item}: {count}")
    else:
        lines.append("- none")

    lines.extend(["", "## Failed Cases", ""])
    if fail_rows:
        lines.extend(["| test_id | benchmark_id | FAIL | fail_items | fail_details | suggestion |", "|---|---:|---:|---|---|---|"])
        for row in fail_rows:
            lines.append(
                f"| {row['test_id']} | {row['benchmark_id']} | {row['FAIL']} | "
                f"{escape(row['fail_items'])} | {escape(row['fail_details'])} | {escape(row['fail_suggestions'])} |"
            )
    else:
        lines.append("- none")

    lines.extend(["", "## Warning-Only Cases", ""])
    if warn_rows:
        lines.extend(["| test_id | benchmark_id | WARN | warn_items | warn_details |", "|---|---:|---:|---|---|"])
        for row in warn_rows:
            lines.append(
                f"| {row['test_id']} | {row['benchmark_id']} | {row['WARN']} | "
                f"{escape(row['warn_items'])} | {escape(row['warn_details'])} |"
            )
    else:
        lines.append("- none")

    return "\n".join(lines) + "\n"


def escape(value: Any) -> str:
    return str(value).replace("|", "\\|").replace("\n", "<br>")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", default=str(DEFAULT_TABLE), help="Path to result.xlsx or cell_line.xlsx")
    parser.add_argument("--package-dir", default="", help="Root directory for benchmark data; defaults to table parent and searches data/<benchmark_id>, h5ad/<benchmark_id>, and <benchmark_id>")
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR), help="Directory for CSV/JSON/Markdown reports")
    parser.add_argument("--limit", type=int, default=0, help="Only check the first N test_id rows")
    args = parser.parse_args()

    table = Path(args.table)
    package_dir = Path(args.package_dir) if args.package_dir else table.parent
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    checker = load_checker()
    test_ids = read_test_ids(table)
    if args.limit:
        test_ids = test_ids[: args.limit]

    rows = []
    for index, test_id in enumerate(test_ids, 1):
        print(f"[{index}/{len(test_ids)}] {test_id}", flush=True)
        rows.append(summarize_case(checker, table, package_dir, test_id))

    csv_path = out_dir / "summary.csv"
    json_path = out_dir / "summary.json"
    md_path = out_dir / "summary.md"
    write_csv(csv_path, rows)
    json_path.write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")
    md_path.write_text(render_markdown(rows, table, package_dir), encoding="utf-8")

    counts = Counter(row["overall"] for row in rows)
    print(f"Done. total={len(rows)} PASS={counts['PASS']} WARN={counts['WARN']} FAIL={counts['FAIL']}")
    print(f"CSV: {csv_path}")
    print(f"JSON: {json_path}")
    print(f"Markdown: {md_path}")
    return 1 if counts["FAIL"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
