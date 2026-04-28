"""
03_length_ratio_correlations.py
Correlations stratified by min/max length ratio bins.

Column: len_ratio_min_vs_max
Bins: [0.0, 0.5)  and  [0.5, 1.0]

For each bin: each PLM model x each metric -> Spearman/Pearson vs TMmin.

Output CSVs:
  results/spearman/length_ratio_spearman_tmmin.csv
  results/pearson/length_ratio_pearson_tmmin.csv
  results/spearman/length_ratio_n_tmmin.csv

Rows = models, columns = <bin_label>_<metric>  (e.g.  0.0_0.5_cosine)
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


BINS = [
    ("0.0_0.5",  0.0, 0.5, True,  False),   # lo <= x < hi
    ("0.5_1.0",  0.5, 1.0, True,  True),    # lo <= x <= hi
]

MODELS  = list(config.MODEL_LABELS.keys())
METRICS = config.METRIC_SUFFIXES
TARGET  = "TMmin"


def build_lenratio_table(df, target):
    col_names = [f"{bl}_{m}" for bl, *_ in BINS for m in METRICS]

    sp_rows, pe_rows, n_rows = [], [], []
    for model in MODELS:
        sp_row, pe_row, n_row = [], [], []
        for bin_label, lo, hi, lo_inc, hi_inc in BINS:
            lo_mask = (df["len_ratio_min_vs_max"] >= lo) if lo_inc else (df["len_ratio_min_vs_max"] > lo)
            hi_mask = (df["len_ratio_min_vs_max"] <= hi) if hi_inc else (df["len_ratio_min_vs_max"] < hi)
            sub = df[lo_mask & hi_mask]
            for metric in METRICS:
                col = f"{metric}_{model}"
                rho, r, n = _corr(sub[target].values.astype(float),
                                   sub[col].values.astype(float))
                sp_row.append(rho); pe_row.append(r); n_row.append(n)
        sp_rows.append(sp_row); pe_rows.append(pe_row); n_rows.append(n_row)

    sp_df = pd.DataFrame(sp_rows, index=MODELS, columns=col_names)
    pe_df = pd.DataFrame(pe_rows, index=MODELS, columns=col_names)
    n_df  = pd.DataFrame(n_rows,  index=MODELS, columns=col_names)
    return sp_df, pe_df, n_df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 60)
    print("03_length_ratio_correlations.py  —  starting")
    print("=" * 60)

    os.makedirs(config.SPEARMAN_DIR, exist_ok=True)
    os.makedirs(config.PEARSON_DIR,  exist_ok=True)

    print(f"Loading data from:\n  {config.DATA_PATH}")
    df = pd.read_csv(config.DATA_PATH)
    print(f"  Loaded {len(df):,} rows")

    print("\nLength-ratio bin boundaries and sizes:")
    for bin_label, lo, hi, lo_inc, hi_inc in BINS:
        lo_mask = (df["len_ratio_min_vs_max"] >= lo) if lo_inc else (df["len_ratio_min_vs_max"] > lo)
        hi_mask = (df["len_ratio_min_vs_max"] <= hi) if hi_inc else (df["len_ratio_min_vs_max"] < hi)
        n = (lo_mask & hi_mask).sum()
        bracket_lo = "[" if lo_inc else "("
        bracket_hi = "]" if hi_inc else ")"
        print(f"  {bracket_lo}{lo}, {hi}{bracket_hi}  n = {n:,}")

    print(f"\nBuilding length-ratio table for target = {TARGET} ...")
    sp_df, pe_df, n_df = build_lenratio_table(df, TARGET)

    sp_path = os.path.join(config.SPEARMAN_DIR, "length_ratio_spearman_tmmin.csv")
    pe_path = os.path.join(config.PEARSON_DIR,  "length_ratio_pearson_tmmin.csv")
    n_sp    = os.path.join(config.SPEARMAN_DIR, "length_ratio_n_tmmin.csv")
    n_pe    = os.path.join(config.PEARSON_DIR,  "length_ratio_n_tmmin.csv")

    sp_df.to_csv(sp_path); print(f"  Saved {sp_path}")
    pe_df.to_csv(pe_path); print(f"  Saved {pe_path}")
    n_df.to_csv(n_sp);     n_df.to_csv(n_pe)

    print("\n  Spearman (TMmin) cosine columns:")
    cosine_cols = [c for c in sp_df.columns if c.endswith("cosine")]
    print(sp_df[cosine_cols].to_string())

    print("\n03_length_ratio_correlations.py  —  DONE")
