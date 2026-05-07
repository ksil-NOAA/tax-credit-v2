# Jupyter notebooks

Analysis and reproduction workflows live under **`ipynb/`**. Start from **`ipynb/Index.ipynb`** or run:

```bash
cd ipynb
jupyter notebook Index.ipynb
```

A narrative index also appears in [ipynb/README.md](../ipynb/README.md).

## Sections

| Folder | Focus |
|--------|--------|
| `mock-community/` | Mockrobiota retrieval, QIIME 2 preprocessing, taxonomy assignment templates, mock accuracy evaluation |
| `cross-validated/` | CV simulation outputs, classification evaluation |
| `novel-taxa/` | Novel-taxa simulation and evaluation |
| `simulated-community/` | Simulated communities, assignment, evaluation |
| `runtime/` | Runtime benchmarking |

## Landmark result notebooks

These match the “Quick Links” in the supplementary README:

- **Mock community performance** — [`mock-community/evaluate-classification-accuracy.ipynb`](../ipynb/mock-community/evaluate-classification-accuracy.ipynb): methods on known mixtures; precision, recall, F-measure vs expected BIOM.
- **Cross-validated performance** — [`cross-validated/evaluate-classification.ipynb`](../ipynb/cross-validated/evaluate-classification.ipynb): classic P/R/F at multiple ranks when queries are held out of the reference.
- **Novel-taxa performance** — [`novel-taxa/evaluate-classification.ipynb`](../ipynb/novel-taxa/evaluate-classification.ipynb): LCA-style correctness; match vs over/under/misclassification.
- **Simulated community performance** — [`simulated-community/evaluate-classification-accuracy.ipynb`](../ipynb/simulated-community/evaluate-classification-accuracy.ipynb): composition recovery on synthetic communities.
- **Runtime comparison** — [`runtime/analysis.ipynb`](../ipynb/runtime/analysis.ipynb).
- **Reference database comparison** — [`mock-community/evaluate-classification-database-comparison.ipynb`](../ipynb/mock-community/evaluate-classification-database-comparison.ipynb): same mock/method, different references.

## Import patterns (updated)

Notebooks have been aligned with the refactored package layout:

- **`evaluate_results`** — `from tax_credit.mock_evaluation import evaluate_results` (plus `eval_framework` for plotting/merge helpers).
- **Novel/CV classification** — `from tax_credit.novel_evaluation import novel_taxa_classification_evaluation` (and `extract_per_level_accuracy` where used).
- **Mockrobiota pipeline** — `mockrobiota_extract`, `mock_denoise`, `mock_transport` (or the combined `process_mocks` facade).

## Loading saved tables in pandas

Use **`pd.read_csv(..., index_col=0)`** instead of deprecated `DataFrame.from_csv`. Mock TSV summaries often need **`sep='\t'`**; novel/CV summary CSVs are usually comma-separated.

## Extending analyses with new methods

As described in the supplementary README: add precomputed results under the expected directory layout, run the same evaluation notebooks, and compare against published baselines. Contributing results back via pull request keeps the benchmark useful for others.

## Additional examples

- **`tax-credit_example.ipynb`** (repo root) — extended example workflow.
- **`tax-credit_mock_example.ipynb`** — mock-focused example.

Static viewing without a local install is possible via **nbviewer** or GitHub’s notebook preview (links may point at upstream repos depending on your fork).
