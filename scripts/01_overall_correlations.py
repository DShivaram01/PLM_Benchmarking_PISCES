"""
01_overall_correlations.py
Computes overall Spearman and Pearson correlations between TM variants
and PLM / BLAST metrics across the full 20,445-pair dataset.

Tables produced
---------------
table1_all_tm_vs_all_metrics      – TM variants x (PLM + BLAST metrics)
table2_all_tm_vs_plm_metrics      – TM variants x PLM metrics only
table3_plms_x_metrics_x_tmmin     – rows=models, cols=metrics, target=TMmin
table4_tmmin_vs_blast             – TMmin x BLAST metrics
table5_tmmin_vs_cosine_all_plms   – TMmin x cosine per PLM
table6_tmmax_vs_cosine_all_plms   – TMmax x cosine per PLM

Each table is saved as <name>_rho.csv + <name>_n.csv in
results/spearman/ and results/pearson/.
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
    """Return (spearman_rho, pearson_r, n) for two 1-D arrays (drop NaN pairs)."""
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


def compute_table(df, row_cols, col_cols):
    """
    Build (rho_df, r_df, n_df) where rows index row_cols and cols index col_cols.
    """
    rho_rows, r_rows, n_rows = [], [], []
    for rc in row_cols:
        rho_row, r_row, n_row = [], [], []
        for cc in col_cols:
            rho, r, n = _corr(df[rc].values.astype(float),
                               df[cc].values.astype(float))
            rho_row.append(rho); r_row.append(r); n_row.append(n)
        rho_rows.append(rho_row); r_rows.append(r_row); n_rows.append(n_row)
    rho_df = pd.DataFrame(rho_rows, index=row_cols, columns=col_cols)
    r_df   = pd.DataFrame(r_rows,   index=row_cols, columns=col_cols)
    n_df   = pd.DataFrame(n_rows,   index=row_cols, columns=col_cols)
    return rho_df, r_df, n_df


def save_table(rho_df, r_df, n_df, name):
    """Save rho/r/n CSVs to both spearman/ and pearson/ directories."""
    for d, df in [(config.SPEARMAN_DIR, rho_df), (config.PEARSON_DIR, r_df)]:
        suffix = "_rho" if d == config.SPEARMAN_DIR else "_r"
        df.to_csv(os.path.join(d, f"{name}{suffix}.csv"))
    # n is the same for both; save once in spearman/
    n_df.to_csv(os.path.join(config.SPEARMAN_DIR, f"{name}_n.csv"))
    n_df.to_csv(os.path.join(config.PEARSON_DIR,  f"{name}_n.csv"))


# ---------------------------------------------------------------------------
# Per-model x per-metric table (Table 3 style)
# ---------------------------------------------------------------------------

def compute_plm_x_metric_table(df, target):
    """
    Rows = models (esm1b, esm2, prostt5, prott5)
    Cols = metric suffixes (cosine, euclid_sim, rbf_sim, manhattan_sim, dot)
    """
    models  = list(config.MODEL_LABELS.keys())
    metrics = config.METRIC_SUFFIXES
    rho_rows, r_rows, n_rows = [], [], []
    for model in models:
        rho_row, r_row, n_row = [], [], []
        for metric in metrics:
            col = f"{metric}_{model}"
            rho, r, n = _corr(df[target].values.astype(float),
                               df[col].values.astype(float))
            rho_row.append(rho); r_row.append(r); n_row.append(n)
        rho_rows.append(rho_row); r_rows.append(r_row); n_rows.append(n_row)
    rho_df = pd.DataFrame(rho_rows, index=models, columns=metrics)
    r_df   = pd.DataFrame(r_rows,   index=models, columns=metrics)
    n_df   = pd.DataFrame(n_rows,   index=models, columns=metrics)
    return rho_df, r_df, n_df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 60)
    print("01_overall_correlations.py  —  starting")
    print("=" * 60)

    os.makedirs(config.SPEARMAN_DIR, exist_ok=True)
    os.makedirs(config.PEARSON_DIR,  exist_ok=True)

    print(f"Loading data from:\n  {config.DATA_PATH}")
    df = pd.read_csv(config.DATA_PATH)
    print(f"  Loaded {len(df):,} rows x {len(df.columns)} columns")

    # ---- Table 1: all TM variants vs all metrics (PLM + BLAST) ----
    print("\nComputing table1_all_tm_vs_all_metrics ...")
    rho, r, n = compute_table(df, config.TM_VARIANTS, config.ALL_METRICS)
    save_table(rho, r, n, "table1_all_tm_vs_all_metrics")
    print(f"  Saved  (shape {rho.shape})")

    # ---- Table 2: all TM variants vs PLM metrics only ----
    print("Computing table2_all_tm_vs_plm_metrics ...")
    rho, r, n = compute_table(df, config.TM_VARIANTS, config.PLM_METRICS)
    save_table(rho, r, n, "table2_all_tm_vs_plm_metrics")
    print(f"  Saved  (shape {rho.shape})")

    # ---- Table 3: PLMs x metrics x TMmin ----
    print("Computing table3_plms_x_metrics_x_tmmin ...")
    rho, r, n = compute_plm_x_metric_table(df, "TMmin")
    save_table(rho, r, n, "table3_plms_x_metrics_x_tmmin")
    print(f"  Saved  (shape {rho.shape})")

    # ---- Table 4: TMmin vs BLAST metrics ----
    print("Computing table4_tmmin_vs_blast ...")
    rho, r, n = compute_table(df, ["TMmin"], config.BLAST_METRICS)
    save_table(rho, r, n, "table4_tmmin_vs_blast")
    print(f"  Saved  (shape {rho.shape})")

    # ---- Table 5: TMmin vs cosine (all 4 PLMs) ----
    print("Computing table5_tmmin_vs_cosine_all_plms ...")
    rho, r, n = compute_table(df, ["TMmin"], config.COSINE_COLS)
    save_table(rho, r, n, "table5_tmmin_vs_cosine_all_plms")
    print(f"  Saved  (shape {rho.shape})")

    # ---- Table 6: TMmax vs cosine (all 4 PLMs) ----
    print("Computing table6_tmmax_vs_cosine_all_plms ...")
    rho, r, n = compute_table(df, ["TMmax"], config.COSINE_COLS)
    save_table(rho, r, n, "table6_tmmax_vs_cosine_all_plms")
    print(f"  Saved  (shape {rho.shape})")

    # ---- Console summary ----
    print("\n--- Overall TMmin Spearman (cosine) ---")
    sp_rho, _, _ = compute_table(df, ["TMmin"], config.COSINE_COLS)
    print(sp_rho.to_string())

    print("\n01_overall_correlations.py  —  DONE")
    print(f"Output dirs:\n  {config.SPEARMAN_DIR}\n  {config.PEARSON_DIR}")
