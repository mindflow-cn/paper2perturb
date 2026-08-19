---
name: annotate-cell-types
description: "Annotate scRNA-seq observations with standardized cell_type labels and select benchmark-relevant cell populations. Use after raw data download, either directly or as part of build-h5ad, when source annotations are missing, incomplete, or inconsistent."
---

# Annotate Cell Types

Annotate cells and extract benchmark target cell types. Run after data download and before final h5ad generation by `build-h5ad`.

## Usage

```bash
python3 skills/annotate-cell-types/scripts/annotate.py \
  --input-dir raw_data/GSE228421 \
  --dataset-id GSE228421 \
  --benchmark-id GSE228421_psoriasis_lesional_skin_risankizumab_scRNA \
  --case-table "path/to/Case-level.csv" \
  --output-root prepared \
  --mode balanced
```

## Preflight check

```bash
python3 skills/annotate-cell-types/scripts/preflight.py
```

## Outputs

- `prepared/datasets/{dataset_id}/annotated.h5ad` — dataset-level, all cells with `cell_type` + `annotation_method`
- `prepared/benchmarks/{benchmark_id}/selected.h5ad` — benchmark-level, only target cell types
- `prepared/benchmarks/{benchmark_id}/annotation_manifest.json` — strategy, methods, counts, warnings

## Dependencies

Required (default Python-only): scanpy, anndata, pandas, numpy, scipy, celltypist

Optional (R adapters, `--mode deep --enable-r-adapters` only): SingleR, Azimuth, scATOMIC.
**R adapters are currently stubs** — they will be skipped even in deep mode with a message in the manifest.
To implement R adapter inline execution, see `scripts/optional_r/*.R`.

## Supported input formats

- 10x MEX (direct or prefixed multi-sample): `matrix.mtx`, `features.tsv`, `barcodes.tsv`
- Cell Ranger output: `outs/filtered_feature_bc_matrix` or `outs/raw_feature_bc_matrix`
- 10x H5: `filtered_feature_bc_matrix.h5`, `raw_feature_bc_matrix.h5`
- h5ad
- CSV/TSV/TXT
- FASTQ (requires external Cell Ranger preprocessing or `--run-cellranger`)

## Fail-fast behavior

- Missing required Python packages: preflight check fails immediately
- benchmark_id not found in Case-level table: raises ValueError
- No non-empty cell_type values for benchmark: raises ValueError
- FASTQ without `--run-cellranger`: fails with clear instructions
- Unknown input format: fails with expected formats listed
- Target cell selection yields 0 cells: fails unless `--allow-empty-selection`
- build-h5ad `--annotate` failure: fails unless `--allow-annotation-fallback`

## Tests

```bash
python3 skills/annotate-cell-types/scripts/test_annotate.py
```
