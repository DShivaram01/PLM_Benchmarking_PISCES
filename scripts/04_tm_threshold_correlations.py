"""
04_tm_threshold_correlations.py
Correlations above and below the TMmin = 0.5 structural homology threshold.

Splits:
  below: TMmin < 0.5
  above: TMmin >= 0.5

For each split: Spearman and Pearson of cosine similarity (all 4 PLMs)
and blast_log_evalue vs TMmin.

Output CSVs:
  results/spearman/tm_threshold_spearman.csv
  results/pearson/tm_threshold_pearson.csv

Rows  = method (esm1b, esm2, prostt5, prott5, blast)
Cols  = TMmin_lt0.5 | TMmin_ge0.5 | delta (ge - lt)
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


METHODS = [
    ("esm1b",   "cosine_esm1b"),
    ("esm2",    "cosine_esm2"),
    ("prostt5", "cosine_prostt5"),
    ("prott5",  "cosine_prott5"),
    ("blast",   "blast_log_evalue"),
]

TARGET = "TMmin"
THRESHOLD = config.TM_THRESHOLD


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 60)
    print("04_tm_threshold_correlations.py  —  starting")
    print("=" * 60)

    os.makedirs(config.SPEARMAN_DIR, exist_ok=True)
    os.makedirs(config.PEARSON_DIR,  exist_ok=True)

    print(f"Loading data from:\n  {config.DATA_PATH}")
    df = pd.read_csv(config.DATA_PATH)
    print(f"  Loaded {len(df):,} rows")

    below = df[df[TARGET] < THRESHOLD]
    above = df[df[TARGET] >= THRESHOLD]
    print(f"\nSplit at TMmin = {THRESHOLD}:")
    print(f"  TMmin < {THRESHOLD}:   n = {len(below):,}")
    print(f"  TMmin >= {THRESHOLD}: n = {len(above):,}")

    sp_rows, pe_rows = [], []
    index_labels = []

    for label, col in METHODS:
        index_labels.append(label)

        rho_lt, r_lt, n_lt = _corr(below[TARGET].values.astype(float),
                                    below[col].values.astype(float))
        rho_ge, r_ge, n_ge = _corr(above[TARGET].values.astype(float),
                                    above[col].values.astype(float))

        sp_rows.append({
            "TMmin_lt0.5": rho_lt,
            "TMmin_ge0.5": rho_ge,
            "delta":       (rho_ge - rho_lt) if (not np.isnan(rho_lt) and not np.isnan(rho_ge)) else np.nan,
            "n_lt0.5":     n_lt,
            "n_ge0.5":     n_ge,
        })
        pe_rows.append({
            "TMmin_lt0.5": r_lt,
            "TMmin_ge0.5": r_ge,
            "delta":       (r_ge - r_lt) if (not np.isnan(r_lt) and not np.isnan(r_ge)) else np.nan,
            "n_lt0.5":     n_lt,
            "n_ge0.5":     n_ge,
        })

    sp_df = pd.DataFrame(sp_rows, index=index_labels)
    pe_df = pd.DataFrame(pe_rows, index=index_labels)

    sp_path = os.path.join(config.SPEARMAN_DIR, "tm_threshold_spearman.csv")
    pe_path = os.path.join(config.PEARSON_DIR,  "tm_threshold_pearson.csv")

    sp_df.to_csv(sp_path); print(f"\nSaved {sp_path}")
    pe_df.to_csv(pe_path); print(f"Saved {pe_path}")

    print("\n--- TM threshold Spearman table ---")
    print(sp_df.to_string(float_format="{:.4f}".format))

    print("\n04_tm_threshold_correlations.py  —  DONE")
