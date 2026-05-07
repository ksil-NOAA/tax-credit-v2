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

Text-based assignment evaluation for **novel-taxa** and **cross-validated** layouts (not BIOM composition tables). Uses `framework_functions.load_prf`, `compute_prf`, `evaluate_classification`, `find_last_common_ancestor`, and `paths.QUERY_*` filenames.

### `novel_taxa_classification_evaluation(results_dirs, expected_results_dir, summary_fp, test_type='novel-taxa')`

| Parameter | Type | Description |
|-----------|------|-------------|
| `results_dirs` | iterable of `str` | Each path must end with `dataset_id/method_id/params_id` (`parse_assignment_results_dir`). Must contain `query_tax_assignments.txt`. |
| `expected_results_dir` | `str` | Must contain `join(expected_results_dir, dataset_id, query_taxa.tsv)` for each dataset. |
| `summary_fp` | `str` | Where to write the summary CSV (pandas default comma separator). |
| `test_type` | `str` | `'novel-taxa'` or `'cross-validated'`. Selects how `dataset_id` is parsed (`parse_novel_dataset_id` vs `parse_cv_dataset_id`). Other values raise `ValueError`. |

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

**Returns:** `DataFrame` with columns `Dataset`, `level`, `iteration`, `Method`, `Parameters`, plus derived metric columns (e.g. `match_ratio` when `mismatch_level_list` is processed).

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
| `compute_prf(exp, obs, test_type='cross-validated', l_range=range(1,7), sample_weight=None, exclude=None)` | `test_type` in `mock`, `novel-taxa`, `cross-validated`. | Micro-averaged P/R/F; CV mode fills length-7 vectors by level. |
| `precision_recall_fscore(exp, obs, sample_weight=None, exclude=None)` | Internal to `compute_prf`. | Match / underclassification / misclassification logic. |
| `evaluate_classification(obs_taxon, exp_taxon)` | String taxonomies. | `'match'`, `'underclassification'`, `'overclassification'`, `'misclassification'`. |
| `find_last_common_ancestor(obs, exp)` | Taxonomy strings. | Index of shallowest mismatch (used for mismatch histograms). |

Simulation and sweep entry points (e.g. `generate_simulated_datasets`, `parameter_sweep`, `recall_novel_taxa_dirs`) remain here; see source docstrings for full parameter lists.

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
- **`simulation_names`:** directory names for CV vs novel trees, `format_*` / `parse_*` for fold IDs (supports hyphenated DB names in novel IDs).

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

## `tax_credit.plotting_functions`

Seaborn/matplotlib helpers for notebooks (boxplots, heatmaps, PCoA, etc.). Import only in analysis contexts; dependencies match the QIIME amplicon environment described in [installation.md](installation.md).

---

## Related

- [Overview](overview.md) — scientific modes (mock / CV / novel).
- [Directory layout](directory-layout.md) — on-disk contracts.
- [Notebooks](notebooks.md) — typical import patterns in `ipynb/`.
