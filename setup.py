#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2014--, tax-credit development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""Distribution metadata for tax-credit.

QIIME 2 amplicon conda environments supply the scientific Python stack; this
package declares no extra PyPI dependencies so ``pip install -e .`` only
registers tax-credit in the environment.
"""

from setuptools import find_packages, setup

# Reference: qiime2-amplicon 2024.10 (released) — see
# https://raw.githubusercontent.com/qiime2/distributions/dev/2024.10/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml
# Plotting uses matplotlib / seaborn (already in that environment); no bokeh.

INSTALL_REQUIRES = []

# Optional: recreate a PyPI scientific stack when *not* using the QIIME 2 amplicon
# conda env (QIIME 2 APIs will still be missing).
EXTRAS_FULL_PYPI = [
    "numpy>=1.26,<2",
    "pandas>=2.2",
    "scipy>=1.13",
    "biom-format>=2.1.14,<3",
    "scikit-bio>=0.5.0",
    "scikit-learn>=1.4",
    "matplotlib>=3.8",
    "seaborn>=0.12",
    "statsmodels>=0.14",
    "ipython>=8.12",
    "jupyter",
    "notebook",
    "ipywidgets",
    "pytest>=8",
]

LONG_DESCRIPTION = """# tax-credit

Systematic benchmarking of taxonomic classification methods.

## Recommended install (QIIME 2 amplicon 2024.10 + tax-credit)

1. Create and activate a conda environment from the official **QIIME 2 amplicon**
   distribution for your platform (Python 3.10). The released dependency snapshot is
   published here (Linux “latest” example):

   https://raw.githubusercontent.com/qiime2/distributions/dev/2024.10/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml

   Platform-specific install files are also linked from the QIIME 2 docs, e.g.:

   - https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-linux-conda.yml
   - https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-osx-conda.yml
   - https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-osx-arm64-conda.yml

   ```bash
   conda env create -n taxCredit-q2-2024.10 --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-osx-conda.yml
   conda activate taxCredit-q2-2024.10
   ```

2. From the `tax-credit` repository root:

   ```bash
   pip install -e .
   ```

   **No additional PyPI packages are required:** numpy, pandas, scipy, biom-format,
   scikit-bio, matplotlib, seaborn, statsmodels, ipython, pytest, Jupyter, etc.
   come from the QIIME 2 amplicon environment.

**Note:** Without the QIIME 2 conda stack, `qiime2` / `q2-*` imports will not work.
For a non-conda venv you can try ``pip install -e ".[full-pypi]"`` (still no QIIME 2).
"""

setup(
    name="tax-credit",
    version="0.0.0-dev",
    license="BSD-3-Clause",
    python_requires=">=3.10,<3.11",
    packages=find_packages(exclude=("*.tests", "*.tests.*", "tests.*", "tests")),
    install_requires=INSTALL_REQUIRES,
    extras_require={
        "full-pypi": EXTRAS_FULL_PYPI,
    },
    author="Nicholas Bokulich",
    author_email="nbokulich@gmail.com",
    description="Systematic benchmarking of taxonomic classification methods",
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    url="https://github.com/ksilnoaa/tax-credit",
    project_urls={
        "QIIME 2": "https://docs.qiime2.org/",
        "QIIME 2 install": "https://docs.qiime2.org/2024.10/install/native/#install-qiime-2-within-a-miniconda-or-anaconda-distribution",
        "QIIME 2 amplicon 2024.10 (conda deps)": "https://raw.githubusercontent.com/qiime2/distributions/dev/2024.10/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml",
    },
)
