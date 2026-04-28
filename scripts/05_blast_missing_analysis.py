"""
05_blast_missing_analysis.py
Characterizes pairs where BLAST returns no hit (NaN blast_log_evalue).

Analyses performed:
  1. Count and fraction of pairs with missing BLAST scores
  2. TM-score distribution: pairs with vs without BLAST hits
  3. Identity distribution for missing BLAST pairs
  4. Spearman correlation for non-missing BLAST pairs only

Output CSVs:
  results/blast_missing_summary.csv
  results/blast_nonmissing_spearman.csv
  results/blast_nonmissing_pearson.csv
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
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 60)
    print("05_blast_missing_analysis.py  —  starting")
    print("=" * 60)

    os.makedirs(config.RESULTS_DIR,   exist_ok=True)
    os.makedirs(config.SPEARMAN_DIR,  exist_ok=True)
    os.makedirs(config.PEARSON_DIR,   exist_ok=True)

    print(f"Loading data from:\n  {config.DATA_PATH}")
    df = pd.read_csv(config.DATA_PATH)
    print(f"  Loaded {len(df):,} rows")

    blast_col   = "blast_log_evalue"
    has_blast   = df[blast_col].notna()
    n_total     = len(df)
    n_missing   = (~has_blast).sum()
    n_present   = has_blast.sum()

    print(f"\nBLAST hit coverage:")
    print(f"  Total pairs:           {n_total:>7,}")
    print(f"  With BLAST hit:        {n_present:>7,}  ({100*n_present/n_total:.1f}%)")
    print(f"  Without BLAST hit:     {n_missing:>7,}  ({100*n_missing/n_total:.1f}%)")

    df_miss    = df[~has_blast]
    df_present = df[has_blast]

    # ---- TM-score stats ----
    print("\nTMmin stats (pairs WITHOUT BLAST hit):")
    print(df_miss["TMmin"].describe().to_string())
    print("\nTMmin stats (pairs WITH BLAST hit):")
    print(df_present["TMmin"].describe().to_string())

    print("\nIdentity stats (pairs WITHOUT BLAST hit):")
    print(df_miss["PercentID"].describe().to_string())

    # ---- Summary CSV ----
    summary = pd.DataFrame({
        "group": ["with_blast", "without_blast"],
        "n":     [n_present, n_missing],
        "frac":  [n_present/n_total, n_missing/n_total],
        "TMmin_mean": [df_present["TMmin"].mean(), df_miss["TMmin"].mean()],
        "TMmin_median": [df_present["TMmin"].median(), df_miss["TMmin"].median()],
        "PercentID_mean": [df_present["PercentID"].mean(), df_miss["PercentID"].mean()],
    })
    sum_path = os.path.join(config.RESULTS_DIR, "blast_missing_summary.csv")
    summary.to_csv(sum_path, index=False)
    print(f"\nSaved {sum_path}")

    # ---- Correlation on non-missing BLAST pairs ----
    print("\nComputing Spearman/Pearson on non-missing BLAST pairs ...")
    all_cols = config.PLM_METRICS + config.BLAST_METRICS
    sp_rows, pe_rows = [], []
    for tm in config.TM_VARIANTS:
        sp_row, pe_row = {}, {}
        for col in all_cols:
            rho, r, n = _corr(df_present[tm].values.astype(float),
                               df_present[col].values.astype(float))
            sp_row[col] = rho
            pe_row[col] = r
        sp_rows.append(sp_row); pe_rows.append(pe_row)

    sp_df = pd.DataFrame(sp_rows, index=config.TM_VARIANTS)
    pe_df = pd.DataFrame(pe_rows, index=config.TM_VARIANTS)

    sp_path = os.path.join(config.SPEARMAN_DIR, "blast_nonmissing_spearman.csv")
    pe_path = os.path.join(config.PEARSON_DIR,  "blast_nonmissing_pearson.csv")
    sp_df.to_csv(sp_path); print(f"Saved {sp_path}")
    pe_df.to_csv(pe_path); print(f"Saved {pe_path}")

    print("\nSpearman (non-missing BLAST pairs) — BLAST metrics vs TMmin:")
    print(sp_df[config.BLAST_METRICS].to_string(float_format="{:.4f}".format))

    print("\n05_blast_missing_analysis.py  —  DONE")
