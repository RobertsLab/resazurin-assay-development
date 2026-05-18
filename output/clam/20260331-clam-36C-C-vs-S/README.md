`output/clam/20260331-clam-36C-C-vs-S`

---

### OVERVIEW

Output files produced by
`scripts/clam/00.00-resazurin-20260331-clam-36C-C-vs-S.Rmd`.
All tabular files are comma-separated values (CSV).

See `data/clam/20260331-clam-36C-C-vs-S/README.md` for experimental context.

---

### FILES

#### `auc_all_metrics.csv`

Per-individual area under the curve (AUC), one row per individual per metric.

| Column | Description |
|--------|-------------|
| `trace_id` | Individual identifier (sample ID when available, otherwise `plate_id_well_id`) |
| `family_id_group` | Family group |
| `treatment_group` | Treatment group |
| `AUC` | Trapezoid-rule AUC integrated across all timepoints |
| `n_timepoints` | Number of finite timepoints used to compute AUC |
| `metric` | Measurement metric (e.g. `metabolism_per_area_cm2_measurement`) |

---

#### `auc_summary.csv`

Group-level descriptive statistics of AUC values.

| Column | Description |
|--------|-------------|
| `metric` | Measurement metric |
| `family_id_group` | Family group |
| `treatment_group` | Treatment group |
| `n` | Number of individuals |
| `mean` | Mean AUC |
| `sd` | Standard deviation |
| `se` | Standard error of the mean |
| `median` | Median AUC |

---

#### `metabolism.csv`

Full per-well per-timepoint data frame. One row per well per timepoint.
Includes raw values, blank correction steps, and size-normalised metabolism.

| Column | Description |
|--------|-------------|
| `plate_id` | Plate identifier |
| `well_id` | Well identifier |
| `time_hr` | Timepoint (hours from start of assay) |
| `value` | Raw fluorescence (RFU) |
| `family_id_group` | Family group |
| `sample_id_group` | Sample identifier |
| `treatment_group` | Treatment group |
| `is_blank` | `TRUE` if the well is a blank control |
| `exclude_from_analysis` | `TRUE` if the well is flagged for exclusion |
| `trace_id` | Individual identifier used in plots and AUC calculations |
| `fold_change` | Fluorescence relative to each well's own T0 reading |
| `mean_blank_fc` | Mean blank fold-change for this plate × timepoint |
| `corrected_fc` | Blank-corrected fold-change (`fold_change − mean_blank_fc`) |
| `metabolism_per_*` | `corrected_fc` divided by each active size-measurement column |

---

#### `pairwise_stats.csv`

Tukey-adjusted pairwise comparisons from AUC linear models.
One row per contrast per metric.

| Column | Description |
|--------|-------------|
| `contrast` | Group contrast label (e.g. `blue control - tweed control`) |
| `estimate` | Estimated difference in mean AUC |
| `SE` | Standard error of the estimate |
| `df` | Degrees of freedom |
| `t.ratio` | t-statistic |
| `p.value` | Tukey-adjusted p-value |
| `comparison` | Comparison type: `family:treatment`, `family`, or `treatment` |
| `metric` | Measurement metric |
