# Directory layout and path helpers

tax-credit assumes **fixed-depth** directory trees for discovering results. Central definitions live in `tax_credit/paths.py` and `tax_credit/simulation_names.py` so notebooks and library code share one contract.

## Mock-style results (observed)

Observed mock (or mock-like) BIOM tables are expected under:

```text
<root>/<dataset_id>/<reference_id>/<method_id>/<parameter_id>/<table files>
```

Default observed table glob pattern: `table*biom` (see `DEFAULT_MOCK_RESULT_TABLE_PATTERN`).

## Expected composition

Expected tables live under:

```text
<root>/<dataset_id>/<reference_id>/expected/<file>
```

Default expected filename pattern: `table.L{level}-taxa.biom` (see `DEFAULT_EXPECTED_TABLE_PATTERN`).

## Per-sequence and sidecar files

Common filenames (constants in `paths.py`) include:

| Constant | Typical use |
|----------|-------------|
| `FEATURE_TABLE_BIOM` | Mock feature table without taxonomy (`feature_table.biom`) |
| `TRUEISH_TAXONOMIES_TSV` | Per-sequence expected labels |
| `REP_SEQS_TAX_ASSIGNMENTS_TXT` / `TAXONOMY_TSV` | Observed taxonomy for rep seqs |
| `QUERY_TAX_ASSIGNMENTS_TXT` / `QUERY_TAXA_TSV` | Novel/CV query vs expected lists |
| `CLASSIFICATION_ACCURACY_LOG_TSV` | Per-read classification log |

## Novel and cross-validated simulations

Under the data directory, simulation roots use names exposed in `simulation_names.py`:

- `DIR_CROSS_VALIDATED` — `cross-validated`
- `DIR_NOVEL_TAXA_SIMULATIONS` — `novel-taxa-simulations`
- Reference databases often under `ref_dbs/`

Fold directory names are generated and parsed with helpers such as:

- `format_cv_fold_dirname` / `parse_cv_dataset_id` — e.g. `<db>-iter<n>`
- `format_novel_fold_dirname` / `parse_novel_dataset_id` — e.g. `<db>-L<level>-iter<n>` (supports hyphenated database IDs like `B1-REF`)

## Assignment result directories (novel / CV)

Paths ending in:

```text
.../<dataset_id>/<method_id>/<params_id>
```

are parsed by `parse_assignment_results_dir`.

## Using the API

Prefer **named helpers** over manual `split`/`join`:

- `mock_observed_tables_glob`, `expected_tables_glob`
- `parse_mock_result_table_path`, `parse_expected_table_path`
- `parse_result_leaf_dir_to_parts`, `parse_taxonomy_map_path_to_dataset_id`

This makes renames and validation errors easier to manage than hard-coded path indices.
