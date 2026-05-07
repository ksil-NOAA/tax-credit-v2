# Installation

tax-credit is a Python library used alongside a scientific stack for BIOM tables, scikit-bio, pandas, and (for many workflows) **QIIME 2**.

## Recommended: QIIME 2 amplicon 2024.10

The package is designed to install **into** an existing QIIME 2 amplicon environment so duplicate PyPI dependencies are avoided.

1. **Create and activate** a conda environment from the official QIIME 2 amplicon distribution for your platform. Example files (adjust for OS/arch):

   - [Linux conda YAML](https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-linux-conda.yml)
   - [macOS Intel](https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-osx-conda.yml)
   - [macOS ARM64](https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-osx-arm64-conda.yml)

   ```bash
   conda env create -n taxCredit-q2-2024.10 --file <url-or-path-to-yaml>
   conda activate taxCredit-q2-2024.10
   ```

2. **Install tax-credit** in editable mode from the repository root:

   ```bash
   cd tax-credit
   pip install -e .
   ```

   `setup.py` declares **no** `install_requires` entries: numpy, pandas, scipy, biom-format, scikit-bio, matplotlib, seaborn, Jupyter, etc. are expected to come from the QIIME environment.

3. **Python version**: `python_requires` is pinned to **3.10.x** to match current QIIME 2 amplicon releases.

### Running tests

Use the same activated environment:

```bash
python -m unittest discover -s tax_credit/tests -v
```

Some modules (e.g. anything importing `biom`, `skbio`, or QIIME 2) need that full stack; lightweight tests for `paths` and `simulation_names` may run with a smaller environment.

## Optional: `full-pypi` extra

If you are **not** using QIIME 2 conda but want a reasonable PyPI stack for partial development or plotting:

```bash
pip install -e ".[full-pypi]"
```

This does **not** install QIIME 2 (`qiime2`, `q2-types`, plugins). Notebooks and functions that call QIIME 2 APIs will still fail without a proper QIIME environment.

## Legacy or alternate methods

Some taxonomy assignment workflows (e.g. QIIME 1-style commands) may require **separate** conda environments. The [taxonomy-assignment-qiime1](../ipynb/mock-community/taxonomy-assignment-qiime1.ipynb) notebook illustrates switching environments for those tools.

## Jupyter

Install `jupyter` or `notebook` via conda in the QIIME environment (typically already present), then:

```bash
cd ipynb
jupyter notebook Index.ipynb
```

See [notebooks.md](notebooks.md) for an overview of analysis notebooks.
