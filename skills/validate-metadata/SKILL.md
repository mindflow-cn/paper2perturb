---
name: validate-metadata
description: "Validate Paper2Perturb spreadsheet and JSON metadata with both deterministic schema rules and semantic paper-evidence review. Use for a PMID, result.xlsx, test_case.json, benchmark/case rows, or any request to verify drug perturbation fields and gene-direction claims against a source paper."
---

# Validate Metadata

## Overview

Run deterministic field validation, then verify whether metadata rows for one PMID are supported by the local paper. Report unsupported, inconsistent, or weakly supported fields.

## Expected Inputs

- User usually provides a PMID, such as `34591417`.
- Primary table to check: `result.xlsx`.
- Optional legacy or cross-check tables: `result_cell_line.xlsx`, `statement_verification_result.xlsx`, `statement_verification.xlsx`, and PMID-specific `test_case.json` files.
- Paper sources normally live under `papers/` or legacy `paper_md/` directories and are named by PMID, commonly:
  - `papers/<pmid>/<pmid>.md`
  - `paper_md/PMID_<pmid>*/PMID_<pmid>*.md`
  - `paper_md/<pmid>/<pmid>.md`
  - matching PDFs such as `paper_md/PMID_<pmid>*.pdf`

## Workflow

1. Read local task instructions if available:
   - Open `skill.md` in the current repository.
   - Follow its output requirement: report content that cannot be found in the paper.

2. Locate all rows for the PMID:
   - Inspect `result.xlsx` first, both `Benchmark-level` and `Case-level` sheets.
   - Match rows where `pmid` equals the user-provided PMID or `benchmark_id` starts with `<pmid>_`.
   - Also inspect `test_cases/test_case.json`, `data/<benchmark_id>/test_case.json`, and `h5ad/<benchmark_id>/test_case.json` if present.
   - Use `statement_verification_result.xlsx` only as a helper for text matching; do not treat it as sufficient evidence for biological direction or logic.

3. Locate paper evidence:
   - Search `paper_md/` for folders/files containing the PMID.
   - Prefer markdown files for text search.
   - If markdown evidence is incomplete, inspect the matching PDF or images if necessary.
   - Use `rg` for keywords from the table, including drug, cell line, accession, conditions, figure/table labels, genes, dose, and time.

4. **Drug perturbation scope check** (LLM-based, not a Python rule):
   - **Goal**: Keep only benchmarks that are same-cell drug treatment before/after (同一类细胞用药前后). Exclude three categories:
     - **Drug resistance testing** (耐药性检测): experiments comparing sensitive vs. resistant lines, drug-adapted cells, drug-tolerant persisters, dose-escalation-derived resistance, or survival under drug selection.
     - **Genetic perturbation** (基因扰动): CRISPR, shRNA, ORF overexpression, or similar non-drug interventions.
     - **Cross-cluster data** (不同聚类的用药数据): comparisons that span different cell clusters, subpopulations, or subtypes rather than a single cell type before and after drug treatment.
   - **How to check**: Read the benchmark's `perturbation_type`, `comparison_type`, `description`, and any `test_case.paper_reference_content` from the paper. Use the LLM to make a semantic judgment — do NOT rely on keyword matching alone. A description that says "tolerance" could still be a valid treatment time course; a description that says "A375 treated with dabrafenib" could still be a dose-escalation resistance study (read the paper to tell). The key question: *Is this the same cell population, before and after drug exposure, measuring the drug's effect (not resistance acquisition)?*
   - **Output**: If the benchmark falls into any excluded category, report it as `不适用` with the reason (e.g. "耐药性检测 — sensitive vs resistant comparison") and skip further checks for that benchmark. Only proceed to step 5 for benchmarks that pass this filter.

5. **Cell type singularity check** (LLM-based):
   - **Goal**: Each benchmark must involve exactly **one cell type** (一种细胞). A benchmark is defined as: one cell type × one drug → gene expression changes. If a row covers multiple distinct cell types, it must be split into separate benchmarks.
   - **What to check**: Read `cell_type_original`, `cell_type_standard`, `sample_system`, `cell_context`, and the paper text. Determine whether the benchmark uses a single homogeneous cell population or mixes multiple different cell types. Examples that violate this rule:
     - Multiple cell lines pooled together in one benchmark (e.g. "A375 and MCF7 treated with drug X")
     - A mixture of different primary cell types (e.g. "T cells and B cells from patient")
     - A co-culture system where two distinct cell types are present (e.g. "tumor cells co-cultured with fibroblasts")
     - "PBMC" without subsetting to a specific cell type
   - **What passes**: A single named cell line (e.g. "A375", "PC9", "MCF7"), a single primary cell type (e.g. "CD8+ T cells", "hepatocytes"), or a single PDX-derived population.
   - **How to check**: Use the LLM to read the paper and determine whether the experimental design involves one cell type or multiple. Do NOT rely on the column value alone — the paper may reveal additional cell types not recorded in the table. The key question: *Is this experiment performed on exactly one cell type, treated with one drug, measuring that cell type's gene response?*
   - **Output**: If multiple cell types are detected, report as `需拆分` with the list of cell types identified, and suggest splitting into separate benchmarks. The benchmark fails this check and should not proceed to detailed field verification until split.

6. Check benchmark-level fields:
   - Confirm `dataset_accession`, `secondary_accession`, `species`, `cell_type_original`, `cell_type_standard`, `tissue`, `source_type`, `cell_context`, `disease`, `platform`, `perturbation_type`, `perturbation_scope`, `perturbation_name`, `description`, `control_type`, `dose_design`, and `time_design` when present.
   - Mark a field as unsupported if it is only an inference and the task requires direct paper support. Example: `tissue = Skin` may be inferred from melanoma, but should be reported if the paper never directly states skin/tissue.

7. Check each case-level row:
   - Confirm `source_location` exists and points to a real figure, table, caption, or text passage.
   - Confirm `original_statement` is present verbatim or is a faithful paraphrase; quote only short snippets.
   - Confirm every `target_genes` entry appears in the cited table/figure/text or is justified by the cited gene set.
   - Confirm `drug`, `cell_type`, `control`, `dose_groups`, `time_groups`, and experimental design are directly stated or clearly derivable from stated values.
   - Confirm `relation` matches the paper logic for the stated comparison, not just a cluster label. If the paper says a cluster is enriched in a different treatment group than the row's comparison, report weak support or mismatch.

## Evidence Standards

- Directly supported: explicitly stated in paper text, tables, captions, or visible figure labels.
- Derivable: computed from explicit paper values, such as half of an IC50. Report as derivable if the exact value is not written.
- Weak support: evidence supports a related cluster/module/signature but not the row's exact comparison, dose, time, or drug condition.
- Unsupported: no matching paper evidence found in local markdown/PDF after targeted search.
- Inconsistent: paper evidence contradicts the spreadsheet/JSON.

## Rule-Based Validation (Python)

Before semantic checks, run deterministic field-format validation with `validate_metadata.py`. This script checks the following rules in pure Python:

### Benchmark-level rules

| Rule | Field | Check |
|------|-------|-------|
| R1 | `benchmark_id` | Must match pattern `<pmid>_<NN>` where NN is a 2-digit zero-padded index starting from 01, e.g. `31806696_01` |
| R2 | `source_type` | Must be one of: `cell_line`, `primary_culture`, `patient_sample`, `organoid`, `PDX` |
| R3 | `smiles` | Must not be empty; if empty, a reason must be recorded |

### Case-level rules

| Rule | Field | Check |
|------|-------|-------|
| R4 | `test_id` | Must follow pattern `<benchmark_id>_<index>`, e.g. `31806696_01_01` |
| R5 | `dose_groups` | Must be a non-empty list; each entry must contain a `mM`, `uM`, `nM`, `mg`, `mg/day`, or `mg/kg` unit |
| R6 | `time_groups` | If `perturb_var == time`, must be a non-empty list; each entry must contain `day` as unit |
| R7 | `perturb_var` | Must be one of: `dose`, `time` |
| R8 | `relation` | Must be one of: `UP`, `DOWN` |

### Run it

```bash
python3 skills/validate-metadata/scripts/validate_metadata.py result.xlsx
python3 skills/validate-metadata/scripts/validate_metadata.py result_fy.xlsx
```

Any violations are printed to stdout with row index, field name, current value, and reason. Fix these before proceeding to LLM-based semantic checks.

## Practical Commands

Use Python/pandas to inspect Excel rows:

```bash
python - <<'PY'
import pandas as pd
pmid = "34591417"
for path in ["result_cell_line.xlsx", "result.xlsx", "statement_verification_result.xlsx"]:
    try:
        xl = pd.ExcelFile(path)
    except FileNotFoundError:
        continue
    print("\\nFILE", path, xl.sheet_names)
    for sheet in xl.sheet_names:
        df = pd.read_excel(path, sheet_name=sheet)
        mask = df.astype(str).apply(lambda col: col.str.contains(pmid, na=False)).any(axis=1)
        if mask.any():
            print("SHEET", sheet)
            print(df.loc[mask].to_string(max_colwidth=220))
PY
```

Search local paper evidence:

```bash
find paper_md -maxdepth 3 -iname "*34591417*"
rg -n "A375|vemurafenib|GSE164897|Figure 3|Table 1|MITF|AXL" paper_md/PMID_34591417*
```

Use `nl -ba <paper.md> | sed -n '<start>,<end>p'` to capture line-numbered evidence for the final answer.

## Output

Respond in Chinese when the user asks in Chinese. Keep the answer concise and list only problems:

- `不适用`: benchmark falls outside the drug-perturbation-only scope (耐药性检测 / 基因扰动 / 不同聚类用药), skip remaining checks.
- `需拆分`: benchmark covers multiple cell types (多种细胞), must be split into separate single-cell-type benchmarks.
- `未直接找到`: exact field/value is not explicitly present in the paper.
- `证据不足`: related evidence exists but does not prove the row's exact claim.
- `不一致`: paper contradicts the row.

For each issue, include:

- row id, such as `benchmark_id` or `test_id`
- field name and current value
- reason
- short evidence reference with local file and line number when available

If everything checked is supported, say that no unsupported content was found and mention any derivable values separately.
