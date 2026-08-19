---
name: build-h5ad
description: "Download public scRNA-seq perturbation data from GEO, ArrayExpress, Zenodo, CELLxGENE, GitHub, and supported sources, then build standardized control.h5ad, ground_truth.h5ad, and test_case.json outputs from result.xlsx metadata. Use after perturbation metadata extraction."
---

# Build H5AD

Download scRNA-seq drug perturbation data from public repositories and convert to standardized h5ad format (control.h5ad + ground_truth.h5ad).

## Usage

```
/build-h5ad <xlsx_path> <ACCESSION>
```

Example:
```
/build-h5ad result.xlsx GSE139129
/build-h5ad result.xlsx E-MTAB-13502
/build-h5ad result.xlsx "Zenodo 7942968"
```

This will:
1. Download raw data to `raw_data/{ACCESSION}/` (auto-detects provider: GEO, ArrayExpress, etc.)
2. Read `result.xlsx` to find benchmarks matching the accession
3. Auto-detect if cell-type annotation is needed (based on `source_type`) and run it if so
4. Convert raw data to control.h5ad + ground_truth.h5ad per benchmark
5. Save output to `data/{benchmark_id}/`
6. Validate results against paper claims (effect size + direction check)

## Workflow

### Step 0: Generate test_case.json from xlsx

```bash
# Generate central JSON + per-benchmark files
python3 skills/build-h5ad/scripts/xlsx2json.py result.xlsx \
    --output test_cases/test_case.json \
    --per-benchmark-dir data/

# Or generate only per-benchmark files
python3 skills/build-h5ad/scripts/xlsx2json.py result.xlsx --per-benchmark-dir data/
```

Converts benchmark metadata and test case definitions from `result.xlsx` to JSON format.

- Reads `Benchmark-level` sheet for benchmark metadata
- Reads `Case-level` sheet for test case definitions
- Matches test cases to benchmarks by `benchmark_id`
- Parses JSON string fields (`target_genes`, `dose_groups`, `time_groups`)
- Writes central JSON to `test_cases/test_case.json` (with `--output`)
- Writes per-benchmark JSON to `data/{benchmark_id}/test_case.json` (with `--per-benchmark-dir`)

Optional flags:
- `--output <path>`: Write central JSON to file (default: stdout)
- `--per-benchmark-dir <dir>`, `-d`: Write per-benchmark files to `{dir}/{benchmark_id}/test_case.json`
- `--gse <GSE_ID>`: Filter benchmarks by `dataset_accession`
- `--benchmark-ids <ids>`: Comma-separated benchmark IDs to include
- `--list-benchmarks`: List all benchmark IDs in xlsx and exit

JSON output fields per benchmark:
- `benchmark_id`, `sample_system`, `tissue`, `source_type`
- `perturbation_type`, `perturbation_name`, `smiles` (optional)
- `default_cell_subset`, `description`
- `test_cases[]`: each with `test_id`, `target_genes`, `perturb_var`, `control`, `dose_groups`, `time_groups`, `relation`

Case granularity: one row per (cell line + drug + conditions), NOT per gene. Multiple genes sharing the same `perturb_var`, `control`, `dose_groups`, `time_groups`, and `relation` are merged into a single case with `target_genes` as a JSON array of all genes. The `test_id` reflects the merged genes (e.g. `MB0001BREAST001_STAT1_CEBPB_IRF7_UP`).

### Step 1: Download raw data

```bash
python3 skills/build-h5ad/scripts/download.py <ACCESSION>
```

Auto-detects the data provider from the accession prefix and downloads processed data files. Supports:

| Provider | Accession pattern | Source |
|----------|-------------------|--------|
| **NCBI GEO** | `GSE\d+` | FTP listing + HTTP download of supplementary files |
| **EBI ArrayExpress / BioStudies** | `E-MTAB-\d+`, `E-GEOD-\d+`, `E-MEXP-\d+` | BioStudies REST API — prioritises processed count matrices |
| **Zenodo** | `Zenodo \d+`, `zenodo.\d+` (DOI suffix) | Zenodo REST API — downloads all files from the record |
| **CZ CELLxGENE Discover** | `cellxgene:<collection_uuid>` | CELLxGENE curation API — downloads h5ad assets from collection |
| **GitHub** | `GitHub: owner/repo` | GitHub API — downloads data files from public repositories (supports Git LFS) |

- Downloads all files to `raw_data/{ACCESSION}/`
- ArrayExpress: skips raw FASTQ/BAM/tar archives, downloads only processed .mtx.gz + .tsv.gz files
- Extracts .tar archives and gunzips .gz files (when converters can't read compressed natively)
- Skips already-downloaded files

Optional flags:
- `--output-dir <dir>`: Custom output directory (default: `raw_data/{ACCESSION}/`)
- `--skip-download`: Only extract existing archives

Backward compatibility: `geo_download.py` still works — it delegates to `download.py`.

### Step 2: Convert to h5ad (with automatic cell-type annotation)

```bash
python3 skills/build-h5ad/scripts/convert.py <xlsx_path> <ACCESSION>
```

Reads the `Benchmark-level` sheet from the xlsx, finds all benchmarks whose `dataset_accession` matches the accession, then for each benchmark:

1. **Check if annotation is needed** — auto-detected from `source_type` AND actual data:
   - `cell_line` → skip annotation (cell lines are homogeneous)
   - If raw data already contains cell type labels (e.g. `cell_type` column in h5ad or metadata CSV) → skip annotation regardless of `source_type`
   - Other `source_type` values without pre-existing annotations → run annotation
   - Missing/unknown `source_type` without pre-existing annotations → conservatively run annotation
2. **If annotation is needed**: runs `annotate-cell-types` to identify cell types and select target cells matching the Case-level `cell_type` values. Produces `selected.h5ad` with only the target cell population.
3. **Split selected cells** into control and treated using Case-level split guidance (`perturb_var`, `control`, `dose_groups`, `time_groups`)
4. **Normalize** — CPM normalization (target_sum=1e4) + log1p transformation
5. **Save** — to `data/{benchmark_id}/control.h5ad` and `data/{benchmark_id}/ground_truth.h5ad`

Optional flags:
- `--raw-dir <dir>`: Custom raw data directory
- `--output-root <dir>`: Custom output root (default: `data/`)
- `--converter <name>`: Force a specific converter
- `--list-converters`: List available converters
- `--annotate`: Force annotation on ALL benchmarks, even `cell_line` ones
- `--no-annotate`: Skip annotation entirely (override auto-detection)
- `--annotation-mode <quick|balanced|deep>`: Annotation thoroughness (default: `balanced`)
- `--case-table <path>`: Path to Case-level CSV (auto-derived from xlsx if omitted)
- `--metadata-xlsx <path>`: Metadata xlsx with Case-level sheet (auto-derived from xlsx if omitted)
- `--allow-annotation-fallback`: Fall back to direct conversion if annotation fails (UNSAFE — skips target cell filtering)
- `--allow-heuristic-fallback`: Allow heuristic control/treat split when Case-level guidance fails

#### Annotation modes

| Mode | `cell_line_or_in_vitro` | Other source types |
|------|------------------------|-------------------|
| `quick` | Metadata label | Marker heuristics only |
| `balanced` (default) | Metadata label | CellTypist + marker heuristics |
| `deep` | Metadata label + sanity | CellTypist + marker heuristics + R adapters |

The strategy category is auto-classified from the Case-level `cell_type` values. The annotation pipeline is deterministic — no LLM is used.

### Step 3: Generate test_case.json from xlsx

```bash
python3 skills/build-h5ad/scripts/xlsx2json.py result.xlsx \
    --output test_cases/test_case.json \
    --per-benchmark-dir data/
```

Converts benchmark + test case metadata from xlsx to JSON format. Writes both the central aggregation (`test_cases/test_case.json`) and per-benchmark files (`data/{benchmark_id}/test_case.json`). See Step 0 for full usage.

### Step 4: Validate against paper claims

```bash
# Default: only validate NEW benchmarks (not already in check_true.csv), append results
python3 skills/build-h5ad/scripts/validate.py data/ --data-root data/ --sample-col orig.ident

# Validate all benchmarks (including already-validated), replace old rows for re-validated ones
python3 skills/build-h5ad/scripts/validate.py data/ --data-root data/ --sample-col orig.ident --all

# Validate specific benchmarks (always runs, replaces old rows)
python3 skills/build-h5ad/scripts/validate.py data/ --benchmark-ids MB0001BREAST001 --data-root data/ --sample-col orig.ident
```

For each benchmark and test case, verifies that the actual gene expression changes in the h5ad files match the direction (UP/DOWN) claimed in the paper.

**Methodology:**
1. Loads `control.h5ad` and `ground_truth.h5ad` per benchmark (backed mode for files >500 MB)
2. Parses `perturb_var` and filters control cells to the specified control label
3. Iterates over each treatment group (dose/time) separately — one row per group per gene
4. Creates pseudobulk replicates by random partitioning
5. Computes Hedges' g effect size with bootstrap 95% CI
6. Checks if the observed direction matches the expected relation from test_case.json

**Output:** `validation_results/check_true.csv`

| Column | Description |
|--------|-------------|
| `benchmark_id` | Benchmark identifier |
| `test_id` | Test case identifier |
| `gene` | Target gene symbol |
| `treatment_group` | Single treatment group value (e.g. "1.0 uM", "3 day") |
| `is_confident` | Whether bootstrap 95% CI excludes zero |
| `expected_direction` | Expected direction from paper (UP/DOWN) |
| `delta_true` | Hedges' g effect size |
| `delta_true_ci_lower` | Bootstrap 95% CI lower bound |
| `delta_true_ci_upper` | Bootstrap 95% CI upper bound |
| `n_pseudobulk` | Number of pseudobulk replicates |
| `cell_count` | Number of cells in this treatment group |
| `is_ok` | Whether observed direction matches expected (True if not confident or direction matches) |

Also saves `validation_results/validation_summary.json` with aggregate statistics per benchmark.

Optional flags:
- `--output-dir <dir>`: Save validation CSV + summary JSON (default: `validation_results/`)
- `--benchmark-ids <ids>`: Validate only specific benchmarks (always runs, replaces old rows for those benchmarks)
- `--all`: Re-validate ALL benchmarks (even those already in check_true.csv); replaces old rows for re-validated benchmarks
- `--k-pseudo <n>`: Number of pseudobulk replicates (default: 20)
- `--B-bootstrap <n>`: Bootstrap iterations (default: 2000)
- `--sample-col <col>`, `-s`: Column in adata.obs to filter by (default: auto-detect; set to `orig.ident` for Seurat-formatted h5ad)

**Default behavior**: Only benchmarks NOT already in `check_true.csv` are validated. Results are **appended** (never overwrite). Re-validating the same benchmark replaces its old rows in the CSV.

### Step 5: Quick verification

Check that:
- `data/{benchmark_id}/control.h5ad` exists (DMSO/vehicle cells)
- `data/{benchmark_id}/ground_truth.h5ad` exists (drug-treated cells)
- `data/{benchmark_id}/test_case.json` exists (benchmark metadata + test cases)
- Both h5ad files contain the same genes
- `adata.obs` contains condition/time/dose metadata
- Raw counts preserved in `adata.layers['counts']`

## Core Normalization

All converters apply the same normalization:

```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
```

Raw UMI counts are preserved in `adata.layers['counts']`.

## Supported Formats

| Format | Example | Description |
|--------|-------------|-------------|
| gse139129 | GSE139129 | Tab-separated UMI count matrices (.txt.gz), one per sample, metadata in filenames |
| gse274905 | GSE274905 | 10x Chromium HDF5/MTX + Cell Hashing metadata (skeleton) |
| gse275330 | GSE275330 | 10x Chromium MTX, per-sample matrix/barcodes/features |
| gse229617 | GSE229617 | Processed count matrices |
| emtab13502 | E-MTAB-13502 | 10x MTX per-condition (_raw_matrix.mtx.gz + barcodes + features), ArrayExpress |

## Adding a New Converter

1. Create `skills/build-h5ad/scripts/converters/gseXXXXX.py`
2. Subclass `BaseConverter` and implement:
   - `detect(raw_dir) -> bool` (staticmethod)
   - `convert(self, raw_dir) -> tuple[AnnData, AnnData]`

That's it — the skill auto-discovers all `BaseConverter` subclasses in `skills/build-h5ad/scripts/converters/`, no registration needed.

Example skeleton:

```python
from pathlib import Path
from .base import BaseConverter

class GSEXXXXXConverter(BaseConverter):
    @staticmethod
    def detect(raw_dir: Path) -> bool:
        # Check for characteristic files
        return bool(list(raw_dir.glob("*.tsv.gz")))

    def convert(self, raw_dir: Path) -> tuple:
        # 1. Read expression matrix
        # 2. Read metadata
        # 3. Align cells
        # 4. Build AnnData
        adata = self.build_adata(...)
        adata = self.normalize(adata)
        control, treated = self.split_by_condition(
            adata, "condition", "control", "treated"
        )
        return control, treated
```

## Output Structure

```
data/
└── MB0001BREAST001/
    ├── control.h5ad         # DMSO control cells
    ├── ground_truth.h5ad    # paclitaxel-treated cells
    └── test_case.json       # benchmark metadata + test cases
```

## Migration: existing data to per-benchmark layout

If you already have data directories with `control.h5ad` + `ground_truth.h5ad` but no `test_case.json`:

```bash
python3 skills/build-h5ad/scripts/migrate_to_per_benchmark.py
```

This reads `test_cases/test_case.json` and writes each benchmark to `data/{benchmark_id}/test_case.json`. Idempotent — safe to re-run. Use `--dry-run` to preview.

## Dependencies

- Python 3.9+
- scanpy, anndata, pandas, openpyxl, scipy
- wget or curl (for downloading)
