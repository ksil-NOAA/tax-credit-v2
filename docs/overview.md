# Overview

tax-credit supports **systematic benchmarking** of marker-gene taxonomic classifiers. It standardizes where results live on disk, how metrics are computed, and how plots summarize comparisons across methods and parameters.

## Evaluation modes

The framework stresses three complementary designs:

### Mock communities

**Known composition.** Sequences come from mixtures of organisms with known taxonomy and relative abundance (e.g. [Mockrobiota](http://caporasolab.us/mockrobiota/)). Observed BIOM tables are compared to **expected** composition tables. This measures performance under realistic sequencing error and community complexity.

Typical outputs: precision, recall, F-measure, taxon detection, and optional per-sequence metrics when `trueish-taxonomies.tsv` and per-sequence assignment files are available.

### Cross-validated reference classification

Simulations split a reference database into **query** (test) and **reference** (training) folds. Labels remain known, so classic **precision and recall** apply at multiple taxonomic levels. `framework_functions.generate_simulated_datasets` can produce one or more layouts via **`simulation_method`** (default: all of the following):

- **`cross-validated-taxa`** (maps to on-disk `cross-validated/`): stratified folds by taxonomic strata; query labels may be trimmed so each expected taxonomy prefix appears somewhere in the training taxonomies. Queries are **not** present in the per-fold reference FASTA / taxonomy tables.

- **`cross-validated-trad`** (`cross-validated-trad/`): random **KFold** splits by sequence ID; query FASTA and `query_taxa.tsv` contain only the test fold. **`ref_seqs.fasta`** and **`ref_taxa.tsv`** are **symbolic links** to the full simulated-reads FASTA and cleaned taxonomy for that database, so the classifier sees **every** sequence (including queries). Evaluation still compares assignments to the held-out query labels. Shared **`ref_seqs.qza`** / **`ref_taxa.qza`** artifacts under `ref_dbs/` reduce duplication across folds (see [directory-layout.md](directory-layout.md)).

The legacy alias **`cross-validated`** means **`cross-validated-taxa`**.

### Novel-taxa simulations

**Queries with no exact match in the reference.** Reference sequences that share the query taxonomy are removed; correct behavior is often assignment to the **last common ancestor** (LCA). Metrics emphasize match vs. overclassification vs. underclassification vs. misclassification, summarized in classification logs and aggregate tables.

## Data flow (high level)

1. **Simulate or acquire** communities and reference data (see `framework_functions.generate_simulated_datasets`).
2. **Run classifiers** (QIIME 2 naive-bayes, BLAST, VSEARCH consensus, BLCA, etc.) and place outputs under the expected directory depth (see [directory-layout.md](directory-layout.md)).
3. **Evaluate** with `mock_community.evaluate_mock_samples` (mock communities) or `novel_evaluation.novel_taxa_classification_evaluation` (novel / CV / CV-trad / self-validated text assignments; set `test_type` accordingly).
4. **Visualize** with `plotting_functions` and `log_plotting`.

In this fork, Tourmaline's `scripts/run_tax_credit.py` performs all four steps from `config_04_tax_credit.yaml`.

## BIOM and QIIME 2

Mock-community evaluation consumes feature tables (BIOM or TSV), ASV sequences, and expected composition and/or per-ASV taxonomies. Simulated modes consume QIIME 2 artifacts (`.qza`) and exported FASTA/TSV produced by `framework_functions`.

## Hardware expectations

Moderate laptops can run small benchmarks, but a full matrix of databases x methods x parameter sets x folds is hundreds of assignment jobs and generally wants **cluster** resources. Start with a reduced matrix; Tourmaline ships an sbatch wrapper (`scripts/sbatch_tourmaline2_step4.sh`) for SLURM.
