"""
Central configuration for PLM Benchmarking PISCES repository.
All paths and constants are defined here — no hardcoded values elsewhere.
"""

import os

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
DATA_PATH = (
    "/home/aubsxd002/sbhattacharya_research/shivaram/datasets/"
    "pisces_benchmark/results/final_tm_blast_plm_metrics_20445.csv"
)

REPO_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(REPO_DIR, "results")
SPEARMAN_DIR = os.path.join(RESULTS_DIR, "spearman")
PEARSON_DIR = os.path.join(RESULTS_DIR, "pearson")

# ---------------------------------------------------------------------------
# Column groups
# ---------------------------------------------------------------------------
TM_VARIANTS = ["TMmin", "TMmax", "TMavg"]

PLM_METRICS = [
    "cosine_esm1b", "euclid_sim_esm1b", "rbf_sim_esm1b",
    "manhattan_sim_esm1b", "dot_esm1b",
    "cosine_esm2", "euclid_sim_esm2", "rbf_sim_esm2",
    "manhattan_sim_esm2", "dot_esm2",
    "cosine_prostt5", "euclid_sim_prostt5", "rbf_sim_prostt5",
    "manhattan_sim_prostt5", "dot_prostt5",
    "cosine_prott5", "euclid_sim_prott5", "rbf_sim_prott5",
    "manhattan_sim_prott5", "dot_prott5",
]

BLAST_METRICS = [
    "blast_log_evalue",
    "blast_bitscore",
    "blast_identity_pct",
    "blast_qcovs",
]

ALL_METRICS = PLM_METRICS + BLAST_METRICS

COSINE_COLS = [
    "cosine_esm1b",
    "cosine_esm2",
    "cosine_prostt5",
    "cosine_prott5",
]

# ---------------------------------------------------------------------------
# Model labels: internal column prefix -> display name
# ---------------------------------------------------------------------------
MODEL_LABELS = {
    "esm1b":   "ESM-1b",
    "esm2":    "ESM-2",
    "prostt5": "ProstT5",
    "prott5":  "ProtT5",
}

# Metric suffixes (column suffix after the underscore after model name)
METRIC_SUFFIXES = ["cosine", "euclid_sim", "rbf_sim", "manhattan_sim", "dot"]

# ---------------------------------------------------------------------------
# Binning parameters
# ---------------------------------------------------------------------------
IDENTITY_BIN_EDGES  = [0, 10, 20, 30]          # PercentID bins
LENGTH_RATIO_EDGES  = [0.0, 0.5, 1.0]          # len_ratio_min_vs_max bins
TM_THRESHOLD        = 0.5                       # split on TMmin
