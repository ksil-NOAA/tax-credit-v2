# tax-credit documentation

**TAX CREdiT** (TAXonomic ClassifieR Evaluation Tool) is a standardized, extensible framework for comparing taxonomic classifiers on mock communities, cross-validated reference data, and novel-taxa simulations.

## Guides

| Document | Description |
|----------|-------------|
| [Installation](installation.md) | Recommended QIIME 2 amplicon conda setup, `pip install`, optional PyPI extras |
| [Overview](overview.md) | Evaluation modes, scientific goals, and how analyses fit together |
| [Directory layout](directory-layout.md) | On-disk conventions, `paths` and `simulation_names` helpers |
| [Python API](python-api.md) | Main packages and import patterns after the recent refactor |
| [Jupyter notebooks](notebooks.md) | Notebook index, workflows, and extending analyses |

## Quick start

1. Create a [QIIME 2 amplicon](https://docs.qiime2.org/) environment (e.g. 2024.10, Python 3.10).
2. From the repository root: `pip install -e .`
3. Open `ipynb/Index.ipynb` or a specific analysis notebook under `ipynb/`.

## Further reading

- Repository [README](../README.md) (clone URL, hardware notes, Jupyter launch)
- Supplementary notebook introduction [ipynb/README.md](../ipynb/README.md)
- License: [COPYING.txt](../COPYING.txt)
