---
name: validate-h5ad
description: "Run rule-based QC across Paper2Perturb h5ad files, test_case.json, and result.xlsx metadata. Use for a test_id or batch when checking orig.ident, perturb_var, dose_groups, time_groups, cell types, drug/control metadata, and concrete correction suggestions. Validate control.h5ad obs['orig.ident'] against the fixed label C."
---

# Validate H5AD

Use this skill to check one Paper2Perturb `test_id` against:

- `result.xlsx` (`Benchmark-level` and `Case-level` sheets), or legacy `cell_line.xlsx` (`Benchmark` and `Case` sheets)
- its benchmark directory's `test_case.json`
- `ground_truth.h5ad` and `control.h5ad`

Default table path:

```bash
result.xlsx
```

By default, benchmark files are discovered from `data/<benchmark_id>/`, `h5ad/<benchmark_id>/`, or a legacy package directory.

## Quick Start

Run the bundled checker from the repository root:

```bash
python3 skills/validate-h5ad/scripts/validate_h5ad.py \
  --test-id 32094658_01_01
```

Useful options:

```bash
python3 skills/validate-h5ad/scripts/validate_h5ad.py \
  --table result.xlsx \
  --test-id 32220329_01_02 \
  --format markdown
```

The script searches `data/<benchmark_id>/`, `h5ad/<benchmark_id>/`, and `<package-dir>/<benchmark_id>/`. Override the root only when needed:

```bash
python3 skills/validate-h5ad/scripts/validate_h5ad.py \
  --package-dir . \
  --table result.xlsx \
  --test-id 32094658_01_01
```

## Workflow

1. Locate the `benchmark_id` from `test_id`.
2. Read `result.xlsx` or `cell_line.xlsx` and find the matching benchmark and case rows.
3. Read the benchmark directory's `test_case.json` and find the matching case.
4. Read `ground_truth.h5ad` and `control.h5ad` obs metadata.
5. Apply rule-based checks and report `PASS` / `WARN` / `FAIL`.
6. If a rule fails, provide a concrete modification suggestion.

## Core Rules

Prefer deterministic checks over manual judgment. For example:

```python
if perturb_var == "dose":
    assert set(orig_ident_values) <= set(dose_groups)
if perturb_var == "time":
    assert set(orig_ident_values) <= set(time_groups)
```

Use exact string comparison first. Then report normalized comparisons separately for numeric/unit differences; never silently pass a naked numeric dose/time value.

For the full rule list and field mapping, read `references/rules.md` only when you need to adjust or extend the checker.

## Output Expectations

Report:

- input paths and resolved `benchmark_id`
- table-vs-JSON consistency
- h5ad-vs-JSON consistency
- `cell_type` existence and non-empty values in h5ad obs
- `orig.ident` existence and values
- `control.h5ad.obs['orig.ident']` fixed to `C`; do not compare it with the xlsx `control` value
- perturbation value set relationship (`O == E`, `E < O`, `O < E`, partial overlap, no overlap)
- concrete suggestions such as "set `obs['orig.ident']` to values from `dose_groups`" or "fix `test_case.json` from `cell_line.xlsx` first"

When `test_case.json` disagrees with `result.xlsx`, treat the table as the source of truth and ask to fix JSON before judging h5ad as correct. Apply the same rule to a legacy `cell_line.xlsx` input.
