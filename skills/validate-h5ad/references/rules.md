# h5ad orig.ident QC Rules

## Required Files

- Each benchmark directory must contain `test_case.json`, `ground_truth.h5ad`, and `control.h5ad`.
- `ground_truth.h5ad.obs` must contain non-empty `orig.ident`.

## JSON and Table Consistency

- `test_case.json` top-level `benchmark_id` must match `Benchmark.benchmark_id`.
- Check top-level fields against `Benchmark`: `tissue`, `source_type`, `perturbation_type`, `perturbation_name`, `smiles`, `description`.
- Each JSON case `test_id` must appear in `Case.test_id`.
- Each `Case` row for the benchmark should appear in JSON unless intentionally excluded.
- Check case fields against `Case`: `test_id`, `benchmark_id`, `drug`, `target_genes`, `perturb_var`, `control`, `dose_groups`, `time_groups`, `relation`.
- Parse JSON-like spreadsheet fields before comparing arrays: `target_genes`, `dose_groups`, `time_groups`, `paper_dose_index`, `paper_time_index`.
- If table and JSON disagree, use `result.xlsx` as the source of truth and recommend fixing JSON first. Accept `cell_line.xlsx` only as a legacy input.

## perturb_var

- Allowed values: `dose`, `time`.
- `dose`: `dose_groups` must be non-empty.
- `time`: `time_groups` must be non-empty.
- Other values such as `treatment` or `resistance_status` are unsupported unless the schema is intentionally expanded.

## h5ad Metadata

- `ground_truth.h5ad.obs` and `control.h5ad.obs` must contain a non-empty `cell_type` column.
- `control.h5ad.obs['orig.ident']` must exist, be non-empty, and contain only the fixed label `C`.
- Do not compare `control.h5ad.obs['orig.ident']` with the xlsx/JSON `control` value; that field is checked only for table-vs-JSON metadata consistency.
- Always check `orig.ident`.
- Also inspect available perturbation-like obs columns: `dose`, `dose_label`, `condition`, `time`, `treatment_time`, `sample`, `sample_id`, `Sample`, `batch`.
- For `perturb_var == "dose"`, dose-like h5ad values should map to `dose_groups`.
- For `perturb_var == "time"`, time-like h5ad values should map to `time_groups`.
- Values that are drug names, generic conditions (`treated`, `T`, `Pan`, `OSI`), or sample IDs are not valid replacements for dose/time groups without an explicit mapping.

## Set Relationships

For each case, compute:

- `O`: non-empty unique `orig.ident` values from `ground_truth.h5ad`
- `E`: expected values from `dose_groups` or `time_groups`

Classify:

- `O == E`: consistent
- `E < O`: h5ad includes extra groups; may mix other dose/time cells
- `O < E`: h5ad misses expected groups
- partial overlap: report extras and missing values
- no overlap: completely inconsistent

## Units and Normalization

- Dose values must include number and an accepted unit: `mM`, `uM`, `µM`, `nM`, `mg`, `mg/day`, or `mg/kg`. Optional trailing delivery descriptors such as `pellet` are allowed after the unit.
- Time values must include number and unit. Accepted examples: `min`, `minute`, `h`, `hr`, `hour`, `day`, `days`, `week`.
- Preserve original values in reports.
- Normalize only for diagnostics: whitespace, `µM` to `uM`, plural time units, `24 h` to `1 day` when safe.
- Naked numbers are not valid dose/time labels unless a separate source proves the unit.
