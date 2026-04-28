"""
06_summary_leaderboard.py
Builds the full method leaderboard across all TM variants.

Rows: one per (method, metric, TM_variant) combination
  - 4 PLMs x 5 metrics = 20 PLM rows
  - BLAST x 4 metrics   =  4 BLAST rows
  - 3 TM variants
  Total: 24 methods x 3 TM variants = 72 rows

Columns: method, metric, target, spearman, pearson, n

Sorted by target=TMmin, then spearman descending.

Output:
  results/leaderboard_20445.csv
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import pandas as pd
from scipy import stats

import config


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _corr(x, y):
    mask = np.isfinite(x) & np.isfinite(y)
    n = mask.sum()
    if n < 3:
        return np.nan, np.nan, n
    xs, ys = x[mask], y[mask]
    sp = stats.spearmanr(xs, ys)
    pe = stats.pearsonr(xs, ys)
    rho = sp.statistic if hasattr(sp, "statistic") else sp.correlation
    r   = pe.statistic if hasattr(pe, "statistic") else pe[0]
    return rho, r, n


# ---------------------------------------------------------------------------
# Method definitions
# ---------------------------------------------------------------------------

# PLM methods: (display_name, column_name)
PLM_METHODS = []
for model_key in config.MODEL_LABELS:
    for metric in config.METRIC_SUFFIXES:
        col  = f"{metric}_{model_key}"
        name = f"{config.MODEL_LABELS[model_key]}_{metric}"
        PLM_METHODS.append((name, col, model_key, metric))

# BLAST methods
BLAST_METHODS = [
    ("BLAST_log_evalue",    "blast_log_evalue"),
    ("BLAST_bitscore",      "blast_bitscore"),
    ("BLAST_identity_pct",  "blast_identity_pct"),
    ("BLAST_qcovs",         "blast_qcovs"),
]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 60)
    print("06_summary_leaderboard.py  —  starting")
    print("=" * 60)

    os.makedirs(config.RESULTS_DIR, exist_ok=True)

    print(f"Loading data from:\n  {config.DATA_PATH}")
    df = pd.read_csv(config.DATA_PATH)
    print(f"  Loaded {len(df):,} rows")

    rows = []

    # PLM methods
    print("\nComputing PLM correlations ...")
    for display_name, col, model_key, metric in PLM_METHODS:
        for tm in config.TM_VARIANTS:
            rho, r, n = _corr(df[tm].values.astype(float),
                               df[col].values.astype(float))
            rows.append({
                "method":   display_name,
                "model":    model_key,
                "metric":   metric,
                "target":   tm,
                "spearman": rho,
                "pearson":  r,
                "n":        n,
            })

    # BLAST methods
    print("Computing BLAST correlations ...")
    for display_name, col in BLAST_METHODS:
        for tm in config.TM_VARIANTS:
            rho, r, n = _corr(df[tm].values.astype(float),
                               df[col].values.astype(float))
            rows.append({
                "method":   display_name,
                "model":    "BLAST",
                "metric":   col,
                "target":   tm,
                "spearman": rho,
                "pearson":  r,
                "n":        n,
            })

    lb = pd.DataFrame(rows)

    # Sort: TMmin rows first by spearman desc, then other targets
    tmmin = lb[lb["target"] == "TMmin"].sort_values("spearman", ascending=False)
    rest  = lb[lb["target"] != "TMmin"]
    lb_sorted = pd.concat([tmmin, rest], ignore_index=True)

    out_path = os.path.join(config.RESULTS_DIR, "leaderboard_20445.csv")
    lb_sorted.to_csv(out_path, index=False)
    print(f"\nSaved {out_path}  ({len(lb_sorted)} rows)")

    print(f"\n--- Top 10 by TMmin Spearman ---")
    top10 = tmmin.head(10)[["method", "metric", "spearman", "pearson", "n"]]
    print(top10.to_string(index=False, float_format="{:.4f}".format))

    print(f"\nTotal leaderboard rows: {len(lb_sorted)}")
    print("\n06_summary_leaderboard.py  —  DONE")
