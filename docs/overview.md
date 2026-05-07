# Overview

tax-credit supports **systematic benchmarking** of marker-gene taxonomic classifiers. It standardizes where results live on disk, how metrics are computed, and how plots summarize comparisons across methods and parameters.

## Evaluation modes

The framework stresses three complementary designs (described in more detail in the supplementary notebooks):

### Mock communities

**Known composition.** Sequences come from mixtures of organisms with known taxonomy and relative abundance (e.g. [Mockrobiota](http://caporasolab.us/mockrobiota/)). Observed BIOM tables are compared to **expected** composition tables. This measures performance under realistic sequencing error and community complexity.

Typical outputs: precision, recall, F-measure, taxon detection, and optional per-sequence metrics when `trueish-taxonomies.tsv` and per-sequence assignment files are available.

### Cross-validated reference classification

**Train / test splits on a reference database.** Query sequences are held out from the reference used for classification; labels are still known, so classic **precision and recall** apply at multiple taxonomic levels.

### Novel-taxa simulations

**Queries with no exact match in the reference.** Reference sequences that share the query taxonomy are removed; correct behavior is often assignment to the **last common ancestor** (LCA). Metrics emphasize match vs. overclassification vs. underclassification vs. misclassification, summarized in classification logs and aggregate tables.

## Data flow (high level)

1. **Simulate or acquire** communities and reference data (see `framework_functions`, `simulated_communities`, and mock-community notebooks).
2. **Run classifiers** (QIIME 2, BLAST, legacy QIIME 1, etc.) and place outputs under the expected directory depth (see [directory-layout.md](directory-layout.md)).
3. **Evaluate** with `mock_evaluation.evaluate_results` (mock-style BIOM) or `novel_evaluation.novel_taxa_classification_evaluation` (novel/CV text assignments).
4. **Visualize** with `plotting_functions` and pandas/seaborn in notebooks.

## BIOM and QIIME 2

Mock-style evaluation is built around **BIOM** tables and observation metadata (taxonomy). Many pipelines produce QIIME 2 artifacts; tax-credit often consumes exported `.biom` files or uses QIIME 2 in **mock community preprocessing** notebooks (demux, DADA2, etc.).

## Hardware expectations

Moderate laptops can run many steps; large classifier sweeps or millions of reads may need **cluster** resources or long runtimes. See the repository [README](../README.md) for original hardware notes.
