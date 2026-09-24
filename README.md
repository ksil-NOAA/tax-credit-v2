# tax-credit (Tourmaline fork)

**A modified version of [tax-credit](https://github.com/caporaso-lab/tax-credit) (TAXonomic ClassifieR Evaluation Tool), adapted to run as the reference-database benchmarking step of [Tourmaline 2](https://github.com/aomlomics/tourmaline).**

This is **not** the upstream tax-credit repository. It is a fork maintained for use with Tourmaline, and it is not intended to reproduce the analyses of the original tax-credit paper. If you want the original framework and its supplementary notebooks, use [caporaso-lab/tax-credit](https://github.com/caporaso-lab/tax-credit).

---

## What this fork is for

Upstream tax-credit is a notebook-driven framework: you open a Jupyter notebook per analysis, edit paths and parameter sweeps, and run the cells. This fork keeps the scientific core of that framework — the simulation, evaluation and plotting code — but the orchestration has moved into Tourmaline, which drives everything from a single YAML config.

In practice that means you do not call this package directly. Tourmaline's `scripts/run_tax_credit.py` imports it and runs the whole benchmark: it builds the simulated datasets, emits and executes the classification jobs, scores the assignments, and writes summary tables and plots.

## How it differs from upstream

- **Driven by Tourmaline, not notebooks.** The `ipynb/` tree of supplementary notebooks from the original paper has been removed, along with the precomputed `analyses/` results that backed them.
- **Trimmed to what Tourmaline uses.** The legacy notebook-only modules (`eval_framework`, `mock_evaluation`, `biom_cache`, `process_mocks`, `mock_denoise`, `mock_transport`, `mock_quality`, `mockrobiota_extract`, `simulated_communities`) have been removed. See [Package layout](#package-layout) for what remains.
- **New evaluation modes.** Traditional random K-fold cross-validation (`cross-validated-trad`) has been added alongside the original taxonomy-stratified folds, with shared reference artifacts to avoid duplicating a database per fold.
- **New mock-community implementation.** `tax_credit.mock_community` replaces the old `mock_evaluation` / `eval_framework` scoring path, and works from feature tables, ASV sequences and expected composition or per-ASV "trueish" taxonomies.
- **Log analysis and plot theming.** `tax_credit.log_analysis`, `tax_credit.log_plotting` and `tax_credit.plot_theme` were added to summarize per-taxon classifier behaviour and to give the generated figures a consistent look.
- **Modernized environment.** Targets QIIME 2 amplicon 2024.10 (Python 3.10). The upstream Travis CI config (Python 3.5) has been removed.

## Installation

This package expects the scientific stack from a QIIME 2 amplicon environment; it declares no PyPI dependencies of its own.

```bash
conda env create -n qiime2-amplicon-2024.10 \
  --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-osx-conda.yml
conda activate qiime2-amplicon-2024.10

git clone https://github.com/ksil-NOAA/tax-credit-v2.git tax-credit
cd tax-credit
pip install -e .
```

Clone it as a sibling of your Tourmaline directory — Tourmaline's `tax_credit_package_dir` defaults to `../tax-credit`.

See [docs/installation.md](docs/installation.md) for platform-specific environment files and the non-conda fallback.

## Running a benchmark

Benchmarks are run from Tourmaline, not from this repository:

```bash
conda activate snakemake-tour2
cd /path/to/tourmaline
./tourmaline.sh --step tax-credit --configfile config_04_tax_credit.yaml --cores 8
```

Set `tax_credit_package_dir` in `config_04_tax_credit.yaml` to point at this clone. Tourmaline's [docs/steps/tax_credit.md](https://github.com/aomlomics/tourmaline/blob/main/docs/steps/tax_credit.md) and `docs/configuration.md` are the reference for every available option.

Five evaluation modes are supported. The first four are simulated from the reference database itself; `mock-community` requires real sequencing data you supply.

| Mode | What it does |
|---|---|
| `cross-validated` | Taxonomy-aware K-fold splits; classify held-out sequences. |
| `cross-validated-trad` | Traditional random K-fold splits. |
| `novel-taxa` | Hold out whole taxa, so a query's own taxon is absent from the reference. |
| `self-validated` | Classify the full database against itself; a best-case ceiling. |
| `mock-community` | Classify real reads from communities of known composition. |

Benchmark runs are large — a full matrix of databases × methods × parameter sets × folds is hundreds of assignment jobs. Start with a reduced matrix.

## Package layout

| Module | Role |
|---|---|
| `framework_functions` | Dataset simulation, parameter sweeps, classification command generation |
| `novel_evaluation` | Scoring for novel-taxa and cross-validated assignments; best-run selection |
| `mock_community` | Mock-community scoring (precision/recall, taxon detection, Bray-Curtis) |
| `log_analysis`, `log_plotting` | Per-taxon classifier behaviour summaries and their figures |
| `plotting_functions`, `plot_theme` | Metric plots and shared figure styling |
| `paths`, `simulation_names` | On-disk layout conventions and simulation naming |
| `taxa_manipulator` | Taxonomy string parsing, normalization and filtering |

Full API documentation is in [docs/python-api.md](docs/python-api.md); on-disk conventions are in [docs/directory-layout.md](docs/directory-layout.md).

## Examples

The notebooks in [`examples/`](examples/) are reference material, not part of the Tourmaline workflow:

- [`cross-validated-and-novel-taxa.ipynb`](examples/cross-validated-and-novel-taxa.ipynb) — comparing rCRUX and CRABS reference databases interactively, using the current API.
- [`mock-community-legacy.ipynb`](examples/mock-community-legacy.ipynb) — **historical reference only.** Its evaluation cells import modules that were removed from this fork and will raise `ModuleNotFoundError`. Use Tourmaline's mock-community mode instead.

## Documentation

| Document | Description |
|---|---|
| [docs/index.md](docs/index.md) | Documentation index |
| [docs/installation.md](docs/installation.md) | Environment setup |
| [docs/overview.md](docs/overview.md) | Evaluation modes and how they fit together |
| [docs/directory-layout.md](docs/directory-layout.md) | On-disk conventions |
| [docs/python-api.md](docs/python-api.md) | Python API reference |

## Citation

This fork is derived from tax-credit, the evaluation framework published with QIIME 2's `q2-feature-classifier`. **If you use this software, please cite the original paper:**

> Bokulich NA, Kaehler BD, Rideout JR, Dillon M, Bolyen E, Knight R, Huttley GA, Caporaso JG. Optimizing taxonomic classification of marker-gene amplicon sequences with QIIME 2's q2-feature-classifier plugin. *Microbiome*. 2018;6(1):90. doi:[10.1186/s40168-018-0470-z](https://doi.org/10.1186/s40168-018-0470-z)

```bibtex
@article{bokulich2018optimizing,
  title   = {Optimizing taxonomic classification of marker-gene amplicon
             sequences with {QIIME} 2's q2-feature-classifier plugin},
  author  = {Bokulich, Nicholas A. and Kaehler, Benjamin D. and
             Rideout, Jai Ram and Dillon, Matthew and Bolyen, Evan and
             Knight, Rob and Huttley, Gavin A. and Caporaso, J. Gregory},
  journal = {Microbiome},
  volume  = {6},
  number  = {1},
  pages   = {90},
  year    = {2018},
  doi     = {10.1186/s40168-018-0470-z}
}
```

## License

BSD 3-Clause, unchanged from upstream — see [COPYING.txt](COPYING.txt). Copyright (c) 2014--, tax-credit development team.
