# Python API reference

This page documents the main public entry points after the path / evaluation refactors. Names re-exported from facades (`eval_framework`, `framework_functions`, `process_mocks`) behave the same as the canonical modules listed here.

**Suggested imports**

```python
from tax_credit.mock_evaluation import evaluate_results, compute_mock_results
from tax_credit.novel_evaluation import (
    novel_taxa_classification_evaluation,
    extract_per_level_accuracy,
)
from tax_credit.eval_framework import (
    seek_results,
    get_expected_tables_lookup,
    mount_observations,
    parameter_comparisons,
    merge_expected_and_observed_tables,
)
from tax_credit import paths
from tax_credit import simulation_names
```

---

## `tax_credit.mock_evaluation`

Orchestrates **mock-style** evaluation: observed BIOM tables vs expected BIOM tables, multiple taxonomic levels, optional per-sequence P/R/F. Uses `eval_framework.mount_observations`, `compute_taxon_accuracy`, and optionally `per_sequence_precision`. Collapsed tables and raw feature tables can be **cached** in-process (see `enable_biom_cache`).

### `evaluate_results(...)`

High-level driver: discovers result tables under `results_dirs`, loads or computes metrics, reads/writes the summary TSV at `results_fp`.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `results_dirs` | `list[str]` | (required) | Roots to search; each tree must contain `dataset/reference/method/params/table.biom` (see `paths`). |
| `expected_results_dir` | `str` | (required) | Root containing `dataset/reference/expected/` expected BIOMs. |
| `results_fp` | `str` | (required) | Output (and optional input) path for the **tab-separated** summary table. |
| `mock_dir` | `str` | (required) | Parent of per-dataset folders; each dataset folder must contain `feature_table.biom` for per-sequence logic. |
| `taxonomy_level_range` | iterable of int | `range(2, 7)` | 0-based taxonomy depths passed to `mount_observations` (Greengenes-style levels). |
| `min_count` | `int` | `0` | Minimum count on **observed** tables before collapse (via `mount_observations`). |
| `taxa_to_keep` | `list` or `None` | `None` | If set, restricts observations during filtering (prefix / metadata rules in `filter_table`). |
| `md_key` | `str` | `'taxonomy'` | Observation metadata key for taxonomy on observed tables. |
| `dataset_ids` | `list` or `None` | `None` | Restrict to these dataset IDs. |
| `reference_ids` | `list` or `None` | `None` | Restrict to these reference IDs. |
| `method_ids` | `list` or `None` | `None` | Restrict to these method names. |
| `parameter_ids` | `list` or `None` | `None` | Restrict to these parameter folder names. |
| `subsample` | `bool` | `False` | If `True`, shuffle and take first `size` result rows (debug / smoke tests). |
| `filename_pattern` | `str` | `paths.DEFAULT_EXPECTED_TABLE_PATTERN` | `str.format(level)` for expected BIOM filename; default `table.L{0}-taxa.biom` with level from `get_expected_tables_lookup` (`level=6`). |
| `size` | `int` | `10` | Subsample size when `subsample=True`. |
| `per_seq_precision` | `bool` | `False` | If `True` and `trueish-taxonomies.tsv` exists beside expected BIOM, compute per-sequence P/R/F. |
| `exclude` | `list` | `['other']` | Taxonomy labels excluded from **per-sequence** scoring (`compute_prf` with `test_type='mock'`). |
| `backup` | `bool` | `True` | Before overwrite, copy existing `results_fp` to `results_fp + '.bk'` (when write path is used). |
| `force` | `bool` | `False` | If `True`, recompute even when `results_fp` exists. |
| `append` | `bool` | `False` | Merge behavior with existing file; see docstring matrix (`force` × `append` × filters). |
| `enable_biom_cache` | `bool` | `True` | Reuse mounted BIOMs / feature tables across rows sharing paths. |
| `biom_cache_max_entries` | `int` or `None` | `None` | Optional LRU cap on cache size for one `compute_mock_results` run. |

**Returns:** `pandas.DataFrame` — same rows as written to `results_fp` (tab-separated, first column index).

**`force` / `append`:** The docstring in source spells out four combinations (overwrite vs load-only vs append new result directories vs filter loaded frame). When `force=False` and the file exists, existing numeric results are not recomputed unless `append` adds missing `(Dataset, Reference, Method, Parameters)` tuples.

---

### `compute_mock_results(...)`

Lower-level: given an explicit list of result tuples and an expected-path lookup, returns the metrics `DataFrame` **without** the `evaluate_results` file I/O branches.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `result_tables` | `list[tuple]` | (required) | Each tuple: `(dataset_id, reference_id, method_id, parameter_id, actual_table_fp)`. |
| `expected_table_lookup` | `dict` | (required) | `lookup[dataset_id][reference_id] -> expected_biom_fp`. |
| `results_fp` | `str` | (required) | Passed through for API compatibility; **not** used to read/write inside this function. |
| `mock_dir` | `str` | (required) | Root for `join(mock_dir, dataset_id, FEATURE_TABLE_BIOM)`. |
| `taxonomy_level_range` | iterable | `range(2, 7)` | Levels to evaluate. |
| `min_count` | `int` | `0` | Observed-table filtering threshold. |
| `taxa_to_keep` | `list` or `None` | `None` | Passed to `mount_observations` / `filter_table`. |
| `md_key` | `str` | `'taxonomy'` | Observed taxonomy metadata key. |
| `per_seq_precision` | `bool` | `False` | Enable per-sequence branch. |
| `exclude` | `list` or `None` | `None` | Per-sequence exclude list (defaults handled in caller). |
| `enable_biom_cache` | `bool` | `True` | Use `BiomTableCache` vs `NO_CACHE`. |
| `biom_cache_max_entries` | `int` or `None` | `None` | LRU limit. |

**Output columns**

| Column | Meaning |
|--------|---------|
| `Dataset`, `Level`, `SampleID`, `Reference`, `Method`, `Parameters` | Keys for the evaluation row. |
| `Precision`, `Recall`, `F-measure` | From `per_sequence_precision` when enabled and sidecar files exist; else `-1.0`. |
| `Taxon Accuracy Rate`, `Taxon Detection Rate` | From `compute_taxon_accuracy` (two presence/absence style rates); `-1.0` if `ZeroDivisionError` in that path. |

**Side files for per-sequence metrics:** beside expected BIOM: `trueish-taxonomies.tsv`; beside observed BIOM: `rep_seqs_tax_assignments.txt` or `taxonomy.tsv` (see `paths`).

---

## `tax_credit.novel_evaluation`

Text-based assignment evaluation for **novel-taxa**, **cross-validated**, **cross-validated-trad**, and **self-validated** layouts (not BIOM composition tables). Uses `framework_functions.load_prf`, `compute_prf`, `evaluate_classification`, `find_last_common_ancestor`, and `paths.QUERY_*` filenames.

### `novel_taxa_classification_evaluation(results_dirs, expected_results_dir, summary_fp, test_type='novel-taxa')`

| Parameter | Type | Description |
|-----------|------|-------------|
| `results_dirs` | iterable of `str` | Each path must end with `dataset_id/method_id/params_id` (`parse_assignment_results_dir`). Must contain `query_tax_assignments.txt`. Build with `paths.list_assignment_result_dirs(results_root)` when outputs use the standard four-level sweep tree. |
| `expected_results_dir` | `str` | Must contain `join(expected_results_dir, dataset_id, query_taxa.tsv)` for each dataset. |
| `summary_fp` | `str` | Where to write the summary CSV (pandas default comma separator). |
| `test_type` | `str` | `'novel-taxa'`, `'cross-validated'`, `'cross-validated-trad'`, or `'self-validated'`. Selects how `dataset_id` is parsed (`parse_novel_dataset_id`, `parse_cv_dataset_id`, or `parse_self_validated_dataset_id`). Other values raise `ValueError`. |

**Per directory:** writes `classification_accuracy_log.tsv` under that results dir (`CLASSIFICATION_ACCURACY_LOG_TSV`), appends one summary row, returns the full `DataFrame`.

**Output columns:** `Dataset`, `level`, `iteration`, `Method`, `Parameters`, `match_ratio`, `overclassification_ratio`, `underclassification_ratio`, `misclassification_ratio`, `mismatch_level_list`, `Precision`, `Recall`, `F-measure`.

**Loading saved summaries:** `pd.read_csv(summary_fp, index_col=0)`.

---

### `extract_per_level_accuracy(df, columns=[...])`

Expands summary rows into **per-level** rows for plotting (levels `1..6` in the implementation).

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `df` | `DataFrame` | (required) | Typically the return value of `novel_taxa_classification_evaluation`. |
| `columns` | `list` | `['Precision','Recall','F-measure','mismatch_level_list']` | Which columns to expand; `mismatch_level_list` is converted from string form if needed. |

**Returns:** `DataFrame` with columns `Dataset`, `level`, `iteration`, `Method`, `Parameters`, plus derived metric columns. When `mismatch_level_list` is listed, `match_ratio` is added and equals `Recall` at that level.

### `extract_per_level_classification_ratios_by_fold(results_dirs)` / `extract_per_level_classification_ratios(results_dirs)`

Recompute match / over- / under- / misclassification ratios at levels `1..6` (phylum..species) from each results directory's `classification_accuracy_log.tsv`. The `_by_fold` form returns one row per directory and level with `Dataset`, `novel_level` (novel-taxa simulation level, `<NA>` otherwise), `iteration`, `Method`, `Parameters`, `level` and the four ratios. The other form averages across iterations, grouped by `Dataset`, `novel_level`, `Method`, `Parameters` and `level`.

### `select_best_runs(df, metrics, group_cols=("Dataset",), run_cols=("Method", "Parameters"), tolerance=1e-9)`

For each group and metric, returns the run (method + parameters) with the best mean score across its rows. Metrics in `LOWER_IS_BETTER_METRICS` (mis-, over- and underclassification ratio) are minimised; all others are maximised. Ties within `tolerance` go to the first run sorted by `run_cols`. Output columns: the group columns, `metric`, `direction`, the run columns, `value`, `n_folds`, `n_tied`. Filter `df` to a single level first.

---

## `tax_credit.eval_framework`

### Discovery and path processing

| Function | Purpose |
|----------|---------|
| `find_and_process_result_tables(start_dir, biom_processor=abspath, filename_pattern=DEFAULT_MOCK_RESULT_TABLE_PATTERN)` | Glob observed tables; returns `(dataset_id, reference_id, method_id, parameter_id, processed_path_or_table)`. |
| `find_and_process_expected_tables(start_dir, biom_processor=abspath, filename_pattern=..., level=6)` | Glob expected tables at one collapsed level. |
| `get_expected_tables_lookup(start_dir, biom_processor=abspath, filename_pattern=..., level=6)` | Nested dict `dataset_id -> reference_id -> path` (or processed object). |
| `seek_results(results_dirs, dataset_ids=None, reference_ids=None, method_ids=None, parameter_ids=None)` | Union of `find_and_process_result_tables` over dirs, then filter. Asserts each `results_dir` exists. |

### Table operations

| Function | Key parameters | Notes |
|----------|----------------|-------|
| `mount_observations(table_fp, min_count=0, taxonomy_level=6, taxa_to_keep=None, md_key='taxonomy', normalize=True, clean_obs_ids=True, filter_obs=True)` | Loads BIOM, optional `filter_table`, collapses to `taxonomy_level`, optional `norm`. | Core primitive for mock metrics. |
| `filter_table(table, min_count=0, taxonomy_level=None, taxa_to_keep=None, md_key='taxonomy')` | Observation filter callback for BIOM. | Used inside `mount_observations` when counts/taxa filters apply. |
| `compute_taxon_accuracy(actual_table, expected_table, actual_sample_id=None, expected_sample_id=None)` | Sample-wise presence/absence overlap. | Returns two floats `(p, r)`-style rates. |
| `per_sequence_precision(expected_table_fp, actual_table_fp, feature_table, sample_id, taxonomy_level, exclude=None)` | Per-rep-seq P/R/F for one sample. | Returns `(-1,-1,-1)` if no `trueish-taxonomies.tsv`. |

### Summaries and comparisons

| Function | Key parameters | Returns |
|----------|----------------|---------|
| `get_sample_to_top_params(df, metric, sample_col='SampleID', method_col='Method', dataset_col='Dataset', ascending=False)` | Uses mean absolute deviation from max/min to collect “near-best” parameter sets per method. | Wide `DataFrame` indexed by `(Dataset, SampleID)`. |
| `parameter_comparisons(df, method, metrics=[...], sample_col=..., method_col=..., dataset_col=..., ascending=None)` | Counts how often each parameter set is “top” per metric. | `DataFrame` indexed by parameter id. |
| `filter_df(df_in, column_name=None, values=None, exclude=False)` | Row filter helper. | Filtered frame. |
| `method_by_dataset(df, dataset, sort_field, display_fields, group_by='Dataset', test_field='Method')` | First row per method after sort. | Subframe with `display_fields`. |
| `method_by_dataset_a1` | `functools.partial` of `method_by_dataset` with `sort_field="F-measure"` and fixed display tuple. | Convenience for notebooks. |
| `method_by_reference_comparison(df, group_by='Reference', dataset='Dataset', level_range=range(4,7), ...)` | Nested loops over dataset / level / reference calling `method_by_dataset`. | Concatenated summary. |

### Merging BIOMs (notebooks)

`merge_expected_and_observed_tables(expected_results_dir, results_dirs, md_key='taxonomy', min_count=0, taxonomy_level=6, taxa_to_keep=None, biom_fp=MERGED_TABLE_BIOM, filename_pattern=DEFAULT_EXPECTED_TABLE_PATTERN, dataset_ids=None, reference_ids=None, method_ids=None, parameter_ids=None, force=False)`

| Parameter | Notes |
|-----------|--------|
| `biom_fp` | Output filename under each `dataset/reference/` (default `merged_table.biom`). |
| `force` | If **`False`**, the function calls **`exit()`** with a message (intended to stop accidental “Run all” merges). Set **`force=True`** to generate or overwrite merged tables. |

---

## `tax_credit.framework_functions` (selected)

Large module: simulation generation, parameter sweeps, PRF utilities, QIIME helpers, runtime benchmarking. Functions below are the ones **novel evaluation** depends on.

| Function | Signature highlights | Role |
|----------|---------------------|------|
| `load_prf(obs_fp, exp_fp, level=slice(0,7), sort=True)` | Paths or list-like inputs accepted via `load_taxa`. | Align expected/observed taxon lists. |
| `compute_prf(exp, obs, test_type='cross-validated', l_range=range(1,7), sample_weight=None, exclude=None)` | `test_type` in `mock`, `novel-taxa`, `cross-validated`, `cross-validated-trad`. | Micro-averaged P/R/F; CV modes fill length-7 vectors by level. |
| `precision_recall_fscore(exp, obs, sample_weight=None, exclude=None)` | Internal to `compute_prf`. | Match / underclassification / misclassification logic. |
| `evaluate_classification(obs_taxon, exp_taxon)` | String taxonomies. | `'match'`, `'underclassification'`, `'overclassification'`, `'misclassification'`. |
| `find_last_common_ancestor(obs, exp)` | Taxonomy strings. | Index of shallowest mismatch (used for mismatch histograms). |

Simulation and sweep entry points (e.g. `generate_simulated_datasets`, `parameter_sweep`, `recall_simulated_taxa_dirs`) remain here; see below and source docstrings for full parameter lists.

### `recall_simulated_taxa_dirs(data_dir, databases, iterations, ref_seqs=..., ref_taxa=..., max_level=6, min_level=0, multilevel=True)`

Builds ``(dataset_reference_combinations, reference_dbs)`` for taxonomy-assignment parameter sweeps: fold directory names under ``data_dir`` (novel-taxa ``<db>-L<level>-iter<n>`` when ``multilevel=True``, or cross-validated ``<db>-iter<n>`` when ``multilevel=False``) mapped to ``ref_seqs`` / ``ref_taxa`` paths. Pass the appropriate root (e.g. ``novel_taxa_simulations_root(data_dir)`` or ``cross_validated_root(data_dir)``) as the first argument so paths resolve to the simulated tree you generated.

### `recall_self_validated_dirs(data_dir, databases, ref_seqs=..., ref_taxa=...)`

Same return shape as ``recall_simulated_taxa_dirs``, but for **self-validated** datasets: one ``(database, database)`` pair per reference database under ``self_validated_root(data_dir)`` (no CV fold iterations).

### `generate_self_validated_datasets(dataframe, data_dir, ...)`

Builds one self-validation dataset per reference database: every sequence is classified against the full database (no held-out folds, no ID removal from the reference). Writes under ``self-validated/<database>/``.

### `trad_cv_shared_reference_qzas(project_data_dir, reference_id)`

Returns ``(ref_seqs_qza, ref_taxa_qza)`` paths under ``ref_dbs/<reference_id>/`` for **cross-validated-trad** shared training artifacts (``_trad_cv_shared_ref_seqs.qza``, ``_trad_cv_shared_ref_taxa.qza``).

### `trad_cv_naive_bayes_commands_single_classifier(trad_sim_data_dir, project_data_dir, results_dir, database_names, iterations, method_parameters_combinations, ...)`

For **cross-validated-trad**, builds two shell-command lists: **fit** naive Bayes once per (database, **fit** parameter combo) into ``results_dir/<db>/<db>/<method>/<fit-params>/classifier.qza``, then **classify** each fold’s ``query.qza`` into ``results_dir/<fold-id>/<fold-id>/<method>/<run-id>/`` (same depth as ``parameter_sweep`` with ``multilevel=False``). Run all fit commands before classify. See ``tax-credit_example.ipynb`` (cross-validated assignment section).

**``method_parameters_combinations``:** per method, either a **flat** dict (all keys go to ``fit-classifier-naive-bayes``; classify uses only ``confidence`` and ``classify_n_jobs``), or ``{'fit': {...}, 'classify': {...}}`` where each inner dict maps QIIME flag stems to lists (Cartesian product). Classify flags (e.g. ``p-confidence``, ``p-n-jobs``, ``p-reads-per-batch``) are passed to ``classify-sklearn``; omitted ``p-confidence`` / ``p-n-jobs`` default from the function kwargs. When both sides sweep, result dirs use ``<fit-id>__cls__<classify-id>``.

`generate_simulated_datasets(..., simulation_method=...)` accepts a single value or any combination of:
- `cross-validated-taxa` (original taxonomy-aware CV output under `cross-validated/`)
- `cross-validated-trad` (random KFold CV under `cross-validated-trad/`; each fold keeps only test sequences in `query.fasta` / `query_taxa.tsv`, while `ref_seqs.fasta` and `ref_taxa.tsv` are **symlinks** to the full simulated-reads FASTA and cleaned taxonomy TSV so the reference still includes every sequence. Per-fold `ref_seqs.qza` and `ref_taxa.qza` symlink to shared QIIME artifacts in the ref database directory (`_trad_cv_shared_ref_seqs.qza` and `_trad_cv_shared_ref_taxa.qza`) to save disk space.)
- `novel-taxa` (novel-taxa output under `novel-taxa-simulations/`)

Default behavior generates **all three** simulation types.  
Backward compatibility: `cross-validated` is treated as an alias of `cross-validated-taxa`.

---

## `tax_credit.biom_cache`

Used internally by `compute_mock_results`; exposed for tests or custom tooling.

| Name | Description |
|------|-------------|
| `mount_observations_cache_key(table_fp, min_count, taxonomy_level, taxa_to_keep, md_key, filter_obs)` | Stable tuple key including `realpath` of `table_fp`. |
| `feature_table_cache_key(table_fp)` | Key for raw `load_table` cache. |
| `BiomTableCache(max_entries=None)` | `get_or_put(key, factory)` with optional LRU eviction. |
| `NO_CACHE` | `_NoBiomCache` singleton; always runs `factory`. |

---

## `tax_credit.paths` and `tax_credit.simulation_names`

- **`paths`:** filename constants (`FEATURE_TABLE_BIOM`, `QUERY_TAXA_TSV`, …), default glob patterns, and `parse_*` helpers for path segments. Prefer these over hard-coded strings or `split(sep)[-5]` indexing.

### `assignment_result_leaf_glob(results_root)`

Returns the glob string `join(results_root, '*', '*', '*', '*')` for sweep leaves (`dataset/reference/method/parameters`).

### `list_assignment_result_dirs(results_root, assignments_filename=QUERY_TAX_ASSIGNMENTS_TXT, sort=True)`

Lists directories under *results_root* at that depth that contain the assignments file. Use this instead of a raw `glob` when **cross-validated-trad** (or similar) also has classifier-only directories at the same depth.
- **`simulation_names`:** constants such as `DIR_CROSS_VALIDATED`, `DIR_CROSS_VALIDATED_TRAD`, `DIR_NOVEL_TAXA_SIMULATIONS`, `DIR_REF_DBS`; helpers `cross_validated_root`, `cross_validated_trad_root`, `novel_taxa_simulations_root`, `ref_dbs_root`; `format_*` / `parse_*` for fold IDs (`parse_cv_dataset_id` applies to both CV trees; novel IDs support hyphenated DB names).

See [directory-layout.md](directory-layout.md).

---

## Mock community QIIME 2 pipeline (Phase 5)

| Module | Responsibility |
|--------|----------------|
| `tax_credit.mockrobiota_extract` | Mockrobiota metadata, downloads, expected TSV → BIOM, `amend_biom_taxonomy_ids`. |
| `tax_credit.mock_denoise` | Demux, DADA2, feature table export, optional tree. |
| `tax_credit.mock_transport` | Copy artifacts into repo `data/` layout. |
| `tax_credit.process_mocks` | Re-exports all public functions from the three modules above (`__all__` in source). |

---

## `tax_credit.plot_theme`

Shared plot style. `apply_tax_credit_theme()` sets fonts (Arial, falling back to DejaVu Sans), editable TrueType PDF text, `constrained_layout`, and trimmed axes; every plotting function below calls it. Also provides `METHOD_COLORS` / `method_palette(methods, override=None)` (a fixed Okabe-Ito colour per classify method), `CLASSIFICATION_RATIO_COLORS` and `RATIO_STACK_ORDER`, `metric_cmap(metric)` (`mako_r` for scores, `rocket_r` for error ratios), `metric_limits(values)` (0-1 unless every value is within 0.25 of 0 or 1, then zoomed; returns `(low, high, zoomed)`), and `metric_label` / `eval_method_label` / `ratio_label` for readable names.

## `tax_credit.plotting_functions`

Seaborn/matplotlib helpers (boxplots, heatmaps, PCoA, etc.); dependencies match the QIIME amplicon environment described in [installation.md](installation.md). Every plotting function returns a matplotlib `Figure` and never calls `plt.show()`.

Evaluation metric plots used by the Tourmaline tax-credit step:

| Function | Draws |
|----------|-------|
| `pointplot_from_data_frame(df, x, metric, hue="Method", col="Dataset", x_order=None, col_order=None, palette=None, x_label=None, title=None)` | Mean metric per `x`, one line per `hue`, one panel per `col`; error bars span min-max. |
| `faceted_boxplot_from_data_frame(df, x, metric, hue="Method", col=None, col_order=None, palette=None, title=None)` | Boxplots with each row drawn as a point, one panel per `col`. |
| `heatmap_from_data_frame(df, metric, rows=("Method", "Parameters"), cols=("Dataset",), cmap=None, vmin=None, vmax=None, annotate=None, title=None)` | Mean metric per row/column group; colour map and limits default from the metric; values printed for 120 cells or fewer. |
| `stacked_classification_barplot_from_data_frame(df, run_cols=("Method", "Parameters"), col="Dataset", level_col="level", level_labels=None, level_axis_label=..., title=None)` | Classification ratios by level; one row per run, one column per dataset. |
| `stacked_classification_panels_from_data_frames(panels, ncols, ..., row_labels=None, col_titles=None, panel_size=(2.6, 2.6), title=None)` | Grid of `(title, df)` panels, each one run's ratios by level. |

`tax_credit.log_plotting.method_parameter_sensitivity_heatmap_from_data_frame(pivot_df, title=None, value_label=..., annotate_max_cells=120)` draws one heatmap panel per dataset from a `(dataset, expected_taxonomy)`-indexed pivot, keeping its row order; hatched cells have no data. Rank rows first with `log_analysis.select_top_sensitivity_taxa(pivot_df, top_n)`, which keeps the `top_n` worst taxa per dataset.

**API changes (notebooks in `ipynb/` not yet updated):**

- `show=` removed from every plotting function; figures are returned instead. Call `plt.show()` or display the figure in notebooks.
- `pointplot_from_data_frame` takes one `metric` (was a `y_vars` list), uses `x` / `hue` / `col` (were `x_axis` / `color_by` / `group_by`) and returns a `Figure` (was a dict of `FacetGrid`s).
- `heatmap_from_data_frame`, `boxplot_from_data_frame` and `method_parameter_sensitivity_heatmap_from_data_frame` return a `Figure` (were `Axes`). `boxplot_from_data_frame` no longer fixes the y axis to 0-1 by default.
- `faceted_boxplot_from_data_frame` returns a `Figure` (was a `FacetGrid`) and takes `palette` (was `color_palette`).
- `stacked_classification_barplot_from_data_frame` draws a grid of runs x datasets (was one axis of nested clusters) and returns a `Figure`.
- `select_top_sensitivity_taxa` keeps `top_n` taxa per dataset (was across all datasets) in ranked order.
- Removed: `lmplot_from_data_frame` (did not run on seaborn 0.12) and `DEFAULT_CLASSIFICATION_RATIO_COLORS` (use `plot_theme.CLASSIFICATION_RATIO_COLORS`).

---

## Related

- [Overview](overview.md) — scientific modes (mock / CV / novel).
- [Directory layout](directory-layout.md) — on-disk contracts.
- [Notebooks](notebooks.md) — typical import patterns in `ipynb/`.
