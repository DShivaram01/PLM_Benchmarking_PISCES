# Data

## Source file

The single source of truth for all analyses is:

```
/home/aubsxd002/sbhattacharya_research/shivaram/datasets/pisces_benchmark/results/final_tm_blast_plm_metrics_20445.csv
```

This path is set in `config.py` (`DATA_PATH`). If you are running on a different
machine, update `DATA_PATH` to point to your local copy.

## Dataset description

| Property | Value |
|---|---|
| Pairs | 20,445 |
| Max sequence identity | 30% |
| Duplicates | None |
| Structural similarity | TM-score (TMmin, TMmax, TMavg) |
| Sequence methods | ESM-1b, ESM-2, ProstT5, ProtT5 |
| Sequence baseline | BLAST (log e-value, bitscore, identity %, query coverage) |

## Column reference

| Column | Description |
|---|---|
| `pdbchain1`, `pdbchain2` | PDB chain identifiers |
| `PercentID` | Sequence identity (0–30%) |
| `len1`, `len2` | Chain lengths |
| `len_ratio_min_vs_max` | min(len1,len2) / max(len1,len2) |
| `TMmin`, `TMmax`, `TMavg` | TM-score variants (ground truth) |
| `cosine_<model>` | Cosine similarity of PLM embeddings |
| `euclid_sim_<model>` | Euclidean similarity (transformed distance) |
| `rbf_sim_<model>` | RBF kernel similarity |
| `manhattan_sim_<model>` | Manhattan similarity |
| `dot_<model>` | Dot product of PLM embeddings |
| `blast_log_evalue` | log10(E-value); NaN when no BLAST hit |
| `blast_bitscore` | BLAST bit score |
| `blast_identity_pct` | BLAST alignment identity % |
| `blast_qcovs` | BLAST query coverage % |

## How to obtain

The file is stored on ASC at the path listed above. It is not included in this
repository due to size. Contact the data curator if you need access.
