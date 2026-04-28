"""
02_identity_bin_correlations.py
Correlations stratified by sequence identity bins.

Bins (on PercentID column): [0,10)  [11,20)  [21,30]
For each bin:  each PLM model x each metric -> Spearman/Pearson vs TMmin and TMmax.

Output CSVs (in results/spearman/ and results/pearson/):
  identity_bin_spearman_tmmin.csv
  identity_bin_pearson_tmmin.csv
  identity_bin_spearman_tmmax.csv
  identity_bin_pearson_tmmax.csv

Rows = models, columns = <bin_label>_<metric>  (e.g.  0_10_cosine)
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
    ("0_10",   0,  10),
    ("11_20", 11,  20),
    ("21_30", 21,  30),
]

MODELS  = list(config.MODEL_LABELS.keys())
METRICS = config.METRIC_SUFFIXES


def build_bin_table(df, target):
    """
    Returns (spearman_df, pearson_df, n_df).
    Rows = models, columns = bin_metric combinations.
    """
    col_names = [f"{bl}_{m}" for bl, _, _ in BINS for m in METRICS]

    sp_rows, pe_rows, n_rows = [], [], []
    for model in MODELS:
        sp_row, pe_row, n_row = [], [], []
        for bin_label, lo, hi in BINS:
            mask = (df["PercentID"] >= lo) & (df["PercentID"] <= hi)
            sub  = df[mask]
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
    print("02_identity_bin_correlations.py  —  starting")
    print("=" * 60)

    os.makedirs(config.SPEARMAN_DIR, exist_ok=True)
    os.makedirs(config.PEARSON_DIR,  exist_ok=True)

    print(f"Loading data from:\n  {config.DATA_PATH}")
    df = pd.read_csv(config.DATA_PATH)
    print(f"  Loaded {len(df):,} rows")

    # Print bin sizes
    print("\nIdentity bin sizes:")
    for label, lo, hi in BINS:
        mask = (df["PercentID"] >= lo) & (df["PercentID"] <= hi)
        print(f"  [{lo:>2},{hi:>2}]  n = {mask.sum():,}")

    for target in ["TMmin", "TMmax"]:
        print(f"\nBuilding identity-bin table for target = {target} ...")
        sp_df, pe_df, n_df = build_bin_table(df, target)

        sp_path = os.path.join(config.SPEARMAN_DIR, f"identity_bin_spearman_{target.lower()}.csv")
        pe_path = os.path.join(config.PEARSON_DIR,  f"identity_bin_pearson_{target.lower()}.csv")
        n_sp    = os.path.join(config.SPEARMAN_DIR, f"identity_bin_n_{target.lower()}.csv")
        n_pe    = os.path.join(config.PEARSON_DIR,  f"identity_bin_n_{target.lower()}.csv")

        sp_df.to_csv(sp_path);  print(f"  Saved {sp_path}")
        pe_df.to_csv(pe_path);  print(f"  Saved {pe_path}")
        n_df.to_csv(n_sp);      n_df.to_csv(n_pe)

        print(f"\n  Spearman ({target}) cosine columns:")
        cosine_cols = [c for c in sp_df.columns if c.endswith("cosine")]
        print(sp_df[cosine_cols].to_string())

    print("\n02_identity_bin_correlations.py  —  DONE")
