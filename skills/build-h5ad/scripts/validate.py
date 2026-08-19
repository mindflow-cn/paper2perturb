#!/usr/bin/env python3
"""Validate that gene expression changes in h5ad files match paper claims.

Reads test_case.json, loads control.h5ad + ground_truth.h5ad per benchmark,
and checks whether the observed direction of change for each target gene
matches the expected relation (UP/DOWN) from the paper.

Methodology: pseudobulk aggregation + Hedges' g effect size with bootstrap CI.
Each treatment group (dose/time) is validated separately, producing one row
per gene per group.

Output: check_true.csv with columns:
  benchmark_id, test_id, gene, treatment_group, is_confident,
  expected_direction, delta_true, delta_true_ci_lower, delta_true_ci_upper,
  n_pseudobulk, cell_count, is_ok

Usage:
    python3 validate.py data/ --data-root data/
    python3 validate.py test_cases/test_case.json --data-root data/
    python3 validate.py data/MB0001BREAST001/test_case.json --data-root data/
"""

import argparse
import gc
import json
import sys
from pathlib import Path

import numpy as np
import anndata as ad


# ---- Statistical functions ----

def create_pseudobulk(expr, k_pseudo=20):
    """Create pseudobulk samples by randomly partitioning cells into k groups."""
    n_cells = len(expr)
    k = min(max(k_pseudo, 2), n_cells // 2)
    if k < 2:
        return expr
    indices = np.random.permutation(n_cells)
    pb_list = []
    cpg = n_cells // k
    for i in range(k):
        start = i * cpg
        end = start + cpg if i < k - 1 else n_cells
        pb_list.append(expr[indices[start:end]].mean(axis=0))
    return np.array(pb_list)


def calculate_hedges_g(x1, x2):
    """Hedges' g effect size (x1=treatment, x2=control)."""
    n1, n2 = len(x1), len(x2)
    if n1 < 2 or n2 < 2:
        return np.nan
    s1 = np.std(x1, ddof=1)
    s2 = np.std(x2, ddof=1)
    sp = np.sqrt(((n1 - 1) * s1 ** 2 + (n2 - 1) * s2 ** 2) / (n1 + n2 - 2))
    if sp == 0:
        return 0.0
    d = (np.mean(x1) - np.mean(x2)) / sp
    j = 1 - (3 / (4 * (n1 + n2) - 9))
    return d * j


def bootstrap_hedges_g(x1, x2, B=2000):
    """Bootstrap Hedges' g → (g_obs, ci_lower, ci_upper)."""
    g_obs = calculate_hedges_g(x1, x2)
    n1, n2 = len(x1), len(x2)
    g_vals = []
    for _ in range(B):
        x1b = np.random.choice(x1, size=n1, replace=True)
        x2b = np.random.choice(x2, size=n2, replace=True)
        gb = calculate_hedges_g(x1b, x2b)
        if not np.isnan(gb):
            g_vals.append(gb)
    if not g_vals:
        return g_obs, np.nan, np.nan
    return g_obs, np.percentile(g_vals, 2.5), np.percentile(g_vals, 97.5)


# ---- Obs column helpers ----

def find_col_from_obs(obs, candidates):
    """Find the first matching column in obs."""
    for c in candidates:
        if c in obs.columns:
            return c
    return None


def filter_obs_by_col(obs, col, value):
    """Case-insensitive filter of obs by column value. Returns (mask, n_cells)."""
    mask = obs[col].astype(str).str.strip().str.lower() == value.strip().lower()
    n = int(mask.sum())
    return (mask, n) if n > 0 else (None, 0)


def get_treatment_label(obs, perturb_var):
    """Get human-readable treatment label from ground_truth obs."""
    candidates = {
        'treatment': ['Treatment_group', 'treatment', 'condition', 'group'],
        'resistance_status': ['resistance_status', 'condition'],
        'dose': ['dose', 'dose_label', 'dose_level'],
        'time': ['time', 'timepoint', 'time_point'],
    }.get(perturb_var, ['condition'])
    for c in candidates:
        if c in obs.columns:
            vals = obs[c].unique()
            return '; '.join(sorted(str(v) for v in vals))
    return 'treated'


def resolve_control_mask(obs, perturb_var, control_label, sample_col=None):
    """Return (mask_or_slice, n_cells) for control filtering.

    When sample_col is specified and exists in obs, uses it directly
    (like check_ground_truth.py). Otherwise performs column-aware search.
    """
    # Direct mode: filter by the specified sample column (e.g. orig.ident)
    if sample_col and sample_col in obs.columns:
        mask, n = filter_obs_by_col(obs, sample_col, control_label)
        if n > 0:
            return mask, n
        return slice(None), len(obs)  # fallback: all cells

    if perturb_var in ('treatment', 'resistance_status'):
        return slice(None), len(obs)

    col_map = {
        'time': ['time', 'time_point', 'timepoint'],
        'dose': ['dose', 'dose_label', 'dose_level', 'concentration'],
    }
    cols = col_map.get(perturb_var, [])

    # Strategy 1: direct full-label match
    for col in cols:
        if col in obs.columns:
            mask, n = filter_obs_by_col(obs, col, control_label)
            if n > 0:
                return mask, n

    # Strategy 2: try each word of the control label against metadata columns
    for word in control_label.strip().split():
        for col in cols:
            if col in obs.columns:
                mask, n = filter_obs_by_col(obs, col, word)
                if n > 0:
                    return mask, n

    return slice(None), len(obs)  # fallback: all cells


def resolve_treatment_mask(obs, perturb_var, group_value, sample_col=None):
    """Return (mask_or_slice, n_cells) for treatment group filtering.

    When sample_col is specified and exists in obs, uses it directly.
    """
    # Direct mode: filter by the specified sample column (e.g. orig.ident)
    if sample_col and sample_col in obs.columns:
        mask, n = filter_obs_by_col(obs, sample_col, group_value)
        if n > 0:
            return mask, n
        return slice(None), len(obs)  # fallback: all cells

    if perturb_var == 'treatment':
        return slice(None), len(obs)

    col_map = {
        'time': ['time', 'time_point', 'timepoint'],
        'dose': ['dose', 'dose_label', 'dose_level', 'concentration'],
        'resistance_status': ['resistance_status', 'condition'],
    }
    for col in col_map.get(perturb_var, []):
        if col in obs.columns:
            mask, n = filter_obs_by_col(obs, col, group_value)
            if n > 0:
                return mask, n

    return slice(None), len(obs)  # fallback: all cells


# ---- Validation ----

def validate_benchmark(benchmark, data_root, k_pseudo=20, B_bootstrap=2000, sample_col=None):
    """Validate all test cases for one benchmark. Returns list of record dicts.

    Each record represents one gene × one treatment group combination.
    """
    benchmark_id = benchmark.get("benchmark_id", "")
    test_cases = benchmark.get("test_cases", [])

    data_dir = data_root / benchmark_id
    ctrl_path = data_dir / "control.h5ad"
    gt_path = data_dir / "ground_truth.h5ad"

    if not ctrl_path.exists():
        print(f"  SKIP: {ctrl_path} not found")
        return []
    if not gt_path.exists():
        print(f"  SKIP: {gt_path} not found")
        return []

    total_mb = (ctrl_path.stat().st_size + gt_path.stat().st_size) / 1e6
    print(f"\n{'=' * 60}")
    print(f"Processing {benchmark_id} ({total_mb:.0f} MB)")

    use_backed = total_mb > 500

    if use_backed:
        ctrl = ad.read_h5ad(ctrl_path, backed='r')
        gt = ad.read_h5ad(gt_path, backed='r')
        print(f"  Backed mode. Ctrl: {ctrl.n_obs} cells, GT: {gt.n_obs} cells")
    else:
        ctrl = ad.read_h5ad(ctrl_path)
        gt = ad.read_h5ad(gt_path)
        print(f"  Ctrl: {ctrl.n_obs} cells, GT: {gt.n_obs} cells")

    ctrl_obs = ctrl.obs
    gt_obs = gt.obs
    ctrl_var_names = list(ctrl.var_names)
    gt_var_names = list(gt.var_names)

    # Build gene name -> positional index maps (handles non-unique var_names)
    # Also strip _NCBI_ID suffix (e.g. AKR1C1_123 -> AKR1C1) for matching
    import re
    ctrl_name_to_idx = {}
    for i, name in enumerate(ctrl_var_names):
        if name not in ctrl_name_to_idx:
            ctrl_name_to_idx[name] = i
        base = re.sub(r'_\d+$', '', name)
        if base and base not in ctrl_name_to_idx:
            ctrl_name_to_idx[base] = i
    gt_name_to_idx = {}
    for i, name in enumerate(gt_var_names):
        if name not in gt_name_to_idx:
            gt_name_to_idx[name] = i
        base = re.sub(r'_\d+$', '', name)
        if base and base not in gt_name_to_idx:
            gt_name_to_idx[base] = i

    common_genes = set(ctrl_name_to_idx) & set(gt_name_to_idx)

    records = []

    for tc in test_cases:
        test_id = tc.get('test_id', 'unknown')
        target_genes = tc.get('target_genes', [])
        expected_relation = str(tc.get('relation', '')).upper()
        perturb_var = tc.get('perturb_var', 'treatment')
        control_label = str(tc.get('control', ''))
        dose_groups = tc.get('dose_groups', [])
        time_groups = tc.get('time_groups', [])
        cell_type = tc.get('cell_type', None)

        # Determine treatment groups to iterate
        if perturb_var == 'dose':
            groups = dose_groups
        elif perturb_var == 'time':
            groups = time_groups
        elif perturb_var == 'resistance_status':
            groups = dose_groups
        else:
            groups = []

        if not groups:
            groups = [get_treatment_label(gt_obs, perturb_var)]

        print(f"  {test_id}: {len(target_genes)} genes, groups={groups}")

        # Resolve cell_type mask (per test case, for multi-cell-line benchmarks)
        ctrl_ct_mask = None
        gt_ct_mask = None
        if cell_type and "cell_type" in ctrl_obs.columns and "cell_type" in gt_obs.columns:
            ct_norm = cell_type.strip().lower().replace("_", " ").replace("-", " ")
            ctrl_ct_mask = ctrl_obs["cell_type"].astype(str).str.strip().str.lower() \
                .str.replace("_", " ").str.replace("-", " ") == ct_norm
            gt_ct_mask = gt_obs["cell_type"].astype(str).str.strip().str.lower() \
                .str.replace("_", " ").str.replace("-", " ") == ct_norm
            if ctrl_ct_mask.sum() < 2 or gt_ct_mask.sum() < 2:
                print(f"    cell_type='{cell_type}': insufficient cells "
                      f"(ctrl={int(ctrl_ct_mask.sum())}, gt={int(gt_ct_mask.sum())}), skip")
                continue

        # Resolve control mask
        ctrl_mask, n_ctrl = resolve_control_mask(ctrl_obs, perturb_var, control_label, sample_col)

        # Intersect with cell_type mask if applicable
        if ctrl_ct_mask is not None:
            if isinstance(ctrl_mask, slice):
                ctrl_mask = ctrl_ct_mask.values
            else:
                ctrl_mask = ctrl_mask & ctrl_ct_mask.values
            n_ctrl = int(ctrl_mask.sum())

        # Collect valid genes (present in both ctrl and gt)
        gene_names = [g for g in target_genes if g in common_genes]
        missing = [g for g in target_genes if g not in common_genes]
        if missing:
            print(f"    missing genes: {missing}")
        if not gene_names:
            continue

        ctrl_gene_idx = [ctrl_name_to_idx[g] for g in gene_names]
        gt_gene_idx = [gt_name_to_idx[g] for g in gene_names]

        # Extract control expression once for all groups
        if use_backed:
            ctrl_expr_all = ctrl[ctrl_mask, ctrl_gene_idx].to_memory().X
        else:
            ctrl_expr_all = ctrl[ctrl_mask, ctrl_gene_idx].X
        if hasattr(ctrl_expr_all, 'toarray'):
            ctrl_expr_all = ctrl_expr_all.toarray()

        # Process each treatment group separately
        for grp in groups:
            grp_str = str(grp).strip()

            gt_mask, n_gt = resolve_treatment_mask(gt_obs, perturb_var, grp_str, sample_col)

            # Intersect with cell_type mask if applicable
            if gt_ct_mask is not None:
                if isinstance(gt_mask, slice):
                    gt_mask = gt_ct_mask.values
                else:
                    gt_mask = gt_mask & gt_ct_mask.values
                n_gt = int(gt_mask.sum())

            if n_ctrl < 2 or n_gt < 2:
                print(f"    [{grp_str}] insufficient: ctrl={n_ctrl}, gt={n_gt}, skip")
                continue

            # Extract treatment expression for this group
            if use_backed:
                gt_expr = gt[gt_mask, gt_gene_idx].to_memory().X
            else:
                gt_expr = gt[gt_mask, gt_gene_idx].X
            if hasattr(gt_expr, 'toarray'):
                gt_expr = gt_expr.toarray()

            # Pseudobulk
            pb_ctrl = create_pseudobulk(ctrl_expr_all, k_pseudo)
            pb_gt = create_pseudobulk(gt_expr, k_pseudo)

            for i, gene in enumerate(gene_names):
                delta, ci_lo, ci_hi = bootstrap_hedges_g(
                    pb_gt[:, i], pb_ctrl[:, i], B=B_bootstrap
                )

                confident = not (np.isnan(ci_lo) or np.isnan(ci_hi) or ci_lo <= 0 <= ci_hi)

                if np.isnan(delta):
                    is_ok = False
                elif not confident:
                    is_ok = True
                elif expected_relation == 'UP' and delta > 0:
                    is_ok = True
                elif expected_relation == 'DOWN' and delta < 0:
                    is_ok = True
                else:
                    is_ok = False

                records.append({
                    'benchmark_id': benchmark_id,
                    'test_id': test_id,
                    'gene': gene,
                    'treatment_group': grp_str,
                    'is_confident': confident,
                    'expected_direction': expected_relation,
                    'delta_true': round(delta, 6) if not np.isnan(delta) else np.nan,
                    'delta_true_ci_lower': round(ci_lo, 6) if not np.isnan(ci_lo) else np.nan,
                    'delta_true_ci_upper': round(ci_hi, 6) if not np.isnan(ci_hi) else np.nan,
                    'n_pseudobulk': len(pb_gt),
                    'cell_count': n_gt,
                    'is_ok': is_ok,
                })

            print(f"    [{grp_str}] ctrl={n_ctrl}, gt={n_gt}, genes={len(gene_names)}")
            del gt_expr, pb_gt

    if use_backed:
        del ctrl, gt
    else:
        del ctrl, gt
    gc.collect()

    return records


def load_benchmarks(path: Path) -> list[dict]:
    """Load benchmark definitions from a file or directory.

    Accepts:
    - Central JSON file: {"description": "...", "benchmarks": [...]}
    - Single per-benchmark file: {"benchmark_id": "...", "test_cases": [...]}
    - Directory: scans for */test_case.json per-benchmark files
    """
    if path.is_dir():
        benchmarks = []
        for tc_file in sorted(path.glob("*/test_case.json")):
            with open(tc_file, encoding="utf-8") as f:
                benchmarks.append(json.load(f))
        return benchmarks

    with open(path, encoding="utf-8") as f:
        data = json.load(f)

    if "benchmarks" in data:
        return data["benchmarks"]

    if "benchmark_id" in data:
        return [data]

    return []


def main():
    parser = argparse.ArgumentParser(
        description="Validate h5ad gene expression changes against paper claims"
    )
    parser.add_argument(
        "test_case_json", nargs="?", default="data/",
        help="Path to test_case.json, per-benchmark file, or data/ directory "
             "(default: data/)",
    )
    parser.add_argument(
        "--data-root", "-d", default="data",
        help="Root directory containing {benchmark_id}/control.h5ad + ground_truth.h5ad",
    )
    parser.add_argument(
        "--output-dir", "-o", default="validation_results",
        help="Output directory for check_true.csv (default: validation_results/)",
    )
    parser.add_argument(
        "--k-pseudo", type=int, default=20,
        help="Number of pseudobulk replicates (default: 20)",
    )
    parser.add_argument(
        "--B-bootstrap", type=int, default=2000,
        help="Bootstrap iterations (default: 2000)",
    )
    parser.add_argument(
        "--benchmark-ids", "-b", default=None,
        help="Comma-separated benchmark IDs to validate (default: only new ones)",
    )
    parser.add_argument(
        "--all", action="store_true",
        help="Validate all benchmarks, including those already in check_true.csv",
    )
    parser.add_argument(
        "--sample-col", "-s", default=None,
        help="Column name in adata.obs to filter by (default: auto-detect; "
             "set to 'orig.ident' for Seurat-formatted h5ad)",
    )
    args = parser.parse_args()

    json_path = Path(args.test_case_json)
    if not json_path.exists():
        print(f"ERROR: {json_path} not found.", file=sys.stderr)
        print("Run xlsx2json.py first to generate test_case.json.", file=sys.stderr)
        sys.exit(1)

    data_root = Path(args.data_root)

    all_benchmarks = load_benchmarks(json_path)
    print(f"Loaded {len(all_benchmarks)} benchmark(s) from {json_path}")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "check_true.csv"

    # Determine which benchmarks are already validated
    existing_ids = set()
    if not args.all and csv_path.exists():
        import pandas as pd
        try:
            existing_df = pd.read_csv(csv_path)
        except Exception as e:
            # Back up corrupted file instead of crashing
            import time
            backup = csv_path.with_suffix(f".csv.bak.{int(time.time())}")
            csv_path.rename(backup)
            print(f"WARNING: {csv_path} is corrupted ({e}), backed up to {backup}",
                  file=sys.stderr)
            existing_df = pd.DataFrame()
        existing_ids = set(existing_df["benchmark_id"].unique()) if "benchmark_id" in existing_df.columns else set()
        if existing_ids:
            print(f"Skipping {len(existing_ids)} benchmark(s) already in {csv_path}: "
                  f"{sorted(existing_ids)}")

    # Filter benchmarks
    if args.benchmark_ids:
        id_set = set(x.strip() for x in args.benchmark_ids.split(","))
        benchmarks = [b for b in all_benchmarks
                      if b.get("benchmark_id", "") in id_set]
        print(f"Filtered to {len(benchmarks)} benchmark(s) by --benchmark-ids")
    elif not args.all:
        benchmarks = [b for b in all_benchmarks
                      if b.get("benchmark_id", "") not in existing_ids]
        skipped = len(all_benchmarks) - len(benchmarks)
        if skipped:
            print(f"Skipping {skipped} already-validated benchmark(s). "
                  f"Use --all to re-validate all.")
    else:
        benchmarks = all_benchmarks
        print("--all: re-validating all benchmarks")

    if not benchmarks:
        print("No benchmarks to validate.")
        sys.exit(0)

    # Validate each benchmark
    all_records = []
    for b in benchmarks:
        bid = b.get("benchmark_id", "unknown")
        n_cases = len(b.get("test_cases", []))
        print(f"\n--- {bid} ({n_cases} test case(s)) ---")

        records = validate_benchmark(
            b, data_root,
            k_pseudo=args.k_pseudo,
            B_bootstrap=args.B_bootstrap,
            sample_col=args.sample_col,
        )
        all_records.extend(records)

        if records:
            n_ok = sum(1 for r in records if r["is_ok"])
            n_conf = sum(1 for r in records if r["is_confident"])
            print(f"    Result: {n_ok}/{len(records)} ok, "
                  f"{n_conf}/{len(records)} confident")

    if not all_records:
        print("\nNo records generated.")
        sys.exit(0)

    # Build new dataframe
    import pandas as pd
    df_new = pd.DataFrame(all_records)

    # Rename test_ids to match xlsx naming convention
    test_id_rename = {
        'MB0015HEART001_等6genes_UP': 'MB0015HEART001_EOMES_UP',
        'MB0006PROSTATE001_等7genes_DOWN': 'MB0006PROSTATE001_CD44_DOWN',
    }
    df_new['test_id'] = df_new['test_id'].replace(test_id_rename)

    cols = ['benchmark_id', 'test_id', 'gene', 'treatment_group', 'is_confident',
            'expected_direction', 'delta_true', 'delta_true_ci_lower', 'delta_true_ci_upper',
            'n_pseudobulk', 'cell_count', 'is_ok']
    df_new = df_new[cols].sort_values(['benchmark_id', 'test_id', 'gene', 'treatment_group']).reset_index(drop=True)

    # Append to existing CSV (deduplicate by benchmark_id)
    if csv_path.exists():
        try:
            df_old = pd.read_csv(csv_path)
        except Exception as e:
            import time
            backup = csv_path.with_suffix(f".csv.bak.{int(time.time())}")
            csv_path.rename(backup)
            print(f"WARNING: {csv_path} is corrupted ({e}), backed up to {backup}",
                  file=sys.stderr)
            df_old = pd.DataFrame()
        # Remove old rows for benchmarks we just re-validated
        new_bids = set(df_new['benchmark_id'].unique())
        if "benchmark_id" in df_old.columns:
            df_old = df_old[~df_old['benchmark_id'].isin(new_bids)]
        df = pd.concat([df_old, df_new], ignore_index=True)
        print(f"Appended {len(df_new)} new rows (removed {len(new_bids)} old benchmark(s) if present)")
    else:
        df = df_new

    df = df.sort_values(['benchmark_id', 'test_id', 'gene', 'treatment_group']).reset_index(drop=True)
    df.to_csv(csv_path, index=False)

    # Summary
    print(f"\n{'=' * 60}")
    print("SUMMARY")
    print(f"{'=' * 60}")
    print(f"Total rows: {len(df)}")
    print(f"is_ok rate: {df['is_ok'].sum()}/{len(df)} = {df['is_ok'].sum() / len(df) * 100:.1f}%")
    print(f"is_confident: {df['is_confident'].sum()}/{len(df)} = {df['is_confident'].sum() / len(df) * 100:.1f}%")

    for bid in sorted(df['benchmark_id'].unique()):
        bdf = df[df['benchmark_id'] == bid]
        n_ok = bdf['is_ok'].sum()
        n_conf = bdf['is_confident'].sum()
        groups = sorted(bdf['treatment_group'].unique())
        print(f"  {bid}: {len(bdf)} rows, ok={n_ok}/{len(bdf)}, conf={n_conf}/{len(bdf)}, groups={groups}")

    bad = df[~df['is_ok']]
    if len(bad):
        print(f"\nMISMATCHES ({len(bad)}):")
        for _, r in bad.iterrows():
            print(f"  {r['test_id']}/{r['gene']}/{r['treatment_group']}: "
                  f"delta={r['delta_true']}, conf={r['is_confident']}")

    print(f"\nResults saved to {csv_path}")

    # Save summary JSON
    bench_summary = {}
    for bid in sorted(df['benchmark_id'].unique()):
        bdf = df[df['benchmark_id'] == bid]
        bench_summary[bid] = {
            "total": len(bdf),
            "n_ok": int(bdf["is_ok"].sum()),
            "n_confident": int(bdf["is_confident"].sum()),
            "n_mismatch": int((bdf["is_confident"] & ~bdf["is_ok"]).sum()),
        }

    summary = {
        "total": len(df),
        "n_ok": int(df["is_ok"].sum()),
        "n_confident": int(df["is_confident"].sum()),
        "n_mismatch": int((df["is_confident"] & ~df["is_ok"]).sum()),
        "ok_rate": round(df["is_ok"].sum() / len(df), 4) if len(df) else 0,
        "confidence_rate": round(df["is_confident"].sum() / len(df), 4) if len(df) else 0,
        "by_benchmark": bench_summary,
    }
    summary_path = out_dir / "validation_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"Summary saved to {summary_path}")


if __name__ == "__main__":
    main()
