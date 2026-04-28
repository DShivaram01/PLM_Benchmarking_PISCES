# PLM Benchmarking on PISCES

Reproducible analysis code for the paper:
**"Benchmarking Protein Language Model Embeddings Against Structural Similarity on the PISCES Dataset"**

This repository contains every script needed to regenerate all result tables in the
paper from a single CSV file. One command runs the full pipeline.

---

## Dataset

The benchmark uses **20,445 non-redundant protein chain pairs** drawn from the PISCES
culled sequence list. All pairs satisfy:

- Pairwise sequence identity ≤ 30%
- No duplicate pairs
- TM-score computed with TM-align (TMmin, TMmax, TMavg reported)
- Four PLM embedding similarities: ESM-1b, ESM-2, ProstT5, ProtT5
- Five similarity metrics per PLM: cosine, Euclidean, RBF, Manhattan, dot product
- BLAST baseline: log e-value, bit score, identity %, query coverage

The raw data lives on ASC at:

```
/home/aubsxd002/sbhattacharya_research/shivaram/datasets/pisces_benchmark/results/final_tm_blast_plm_metrics_20445.csv
```

See `data/README.md` for the full column reference.

---

## How to reproduce

### 1. Clone the repository

```bash
git clone <repo-url>
cd PLM_Benchmarking_PISCES
```

### 2. Set the data path (if running outside ASC)

Edit `config.py` and update `DATA_PATH` to point to the CSV on your system.

### 3. Set up the environment (ASC)

```bash
module load python
conda activate plm          # or whichever env has pandas/scipy installed
pip install -r requirements.txt   # if needed
```

### 4. Run everything

```bash
bash reproduce_all.sh
```

All output CSVs are written to the `results/` directory.

---

## Output files

| File | Location | Paper table / figure |
|---|---|---|
| `table1_all_tm_vs_all_metrics_rho.csv` | `results/spearman/` | Table 1 (all metrics) |
| `table2_all_tm_vs_plm_metrics_rho.csv` | `results/spearman/` | Table 2 (PLM only) |
| `table3_plms_x_metrics_x_tmmin_rho.csv` | `results/spearman/` | Table 3 (PLM x metric grid) |
| `table4_tmmin_vs_blast_rho.csv` | `results/spearman/` | Table 4 (BLAST baseline) |
| `table5_tmmin_vs_cosine_all_plms_rho.csv` | `results/spearman/` | Table 5 (cosine, TMmin) |
| `table6_tmmax_vs_cosine_all_plms_rho.csv` | `results/spearman/` | Table 6 (cosine, TMmax) |
| `identity_bin_spearman_tmmin.csv` | `results/spearman/` | Identity-bin analysis |
| `identity_bin_spearman_tmmax.csv` | `results/spearman/` | Identity-bin analysis |
| `length_ratio_spearman_tmmin.csv` | `results/spearman/` | Length-ratio analysis |
| `tm_threshold_spearman.csv` | `results/spearman/` | TM-threshold analysis |
| `blast_nonmissing_spearman.csv` | `results/spearman/` | BLAST coverage analysis |
| `blast_missing_summary.csv` | `results/` | BLAST missing-hit summary |
| `leaderboard_20445.csv` | `results/` | Full method leaderboard |

Every `_rho.csv` (Spearman) has a matching `_r.csv` (Pearson) and `_n.csv`
(pair count) in the same directory. Pearson variants live under `results/pearson/`.

---

## Repository structure

```
PLM_Benchmarking_PISCES/
├── README.md
├── requirements.txt
├── config.py                         # all paths and constants
├── reproduce_all.sh                  # single entry point
│
├── data/
│   └── README.md                     # dataset description and column reference
│
├── results/
│   ├── spearman/                     # Spearman correlation CSVs
│   └── pearson/                      # Pearson correlation CSVs
│
└── scripts/
    ├── 01_overall_correlations.py
    ├── 02_identity_bin_correlations.py
    ├── 03_length_ratio_correlations.py
    ├── 04_tm_threshold_correlations.py
    ├── 05_blast_missing_analysis.py
    └── 06_summary_leaderboard.py
```

---

## ASC-specific notes

- Load Python before activating conda: `module load python`
- The `plm` conda environment contains all required packages.
- The data file is on the shared filesystem; no download needed on ASC.
- Scripts do not require a GPU.

---

## Dependencies

```
pandas>=1.5.0
numpy>=1.23.0
scipy>=1.9.0
```
