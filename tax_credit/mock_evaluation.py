#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2014--, tax-credit development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""Mock-community BIOM evaluation: seek result tables, compute P/R/F vs expected.

Orchestration for ``evaluate_results`` / ``compute_mock_results`` lives here;
table loading, taxonomy collapse, and per-sequence helpers remain in
``eval_framework`` to avoid import cycles.
"""

from os.path import dirname, exists, join
from random import shuffle
from shutil import copy

from biom import load_table
import pandas as pd

from tax_credit.paths import (
    DEFAULT_EXPECTED_TABLE_PATTERN,
    FEATURE_TABLE_BIOM,
    TRUEISH_TAXONOMIES_TSV,
)


def _write_mock_results(mock_results, results_fp, backup=True):
    if backup:
        copy(results_fp, ''.join([results_fp, '.bk']))
    mock_results.to_csv(results_fp, sep='\t')


def _filter_mock_results(mock_results, dataset_ids, reference_ids, method_ids,
                         parameter_ids):
    """Filter mock results on dataset / reference / method / parameter ids."""
    from tax_credit.eval_framework import filter_df

    if dataset_ids:
        mock_results = filter_df(mock_results, 'Dataset', dataset_ids)
    if reference_ids:
        mock_results = filter_df(
            mock_results, 'Reference', reference_ids)
    if method_ids:
        mock_results = filter_df(mock_results, 'Method', method_ids)
    if parameter_ids:
        mock_results = filter_df(
            mock_results, 'Parameters', parameter_ids)
    return mock_results


def _result_row_exists(mock_results, r):
    """Return True if tuple r (dataset, ref, method, params) is in mock_results."""
    return ((mock_results['Dataset'] == r[0]) &
            (mock_results['Reference'] == r[1]) &
            (mock_results['Method'] == r[2]) &
            (mock_results['Parameters'] == r[3])).any()


def _pending_results_for_append(mock_results, results):
    """Rows from *results* not yet represented in *mock_results*."""
    return [r for r in results if not _result_row_exists(mock_results, r)]


def evaluate_results(results_dirs, expected_results_dir, results_fp, mock_dir,
                     taxonomy_level_range=range(2, 7), min_count=0,
                     taxa_to_keep=None, md_key='taxonomy',
                     dataset_ids=None, reference_ids=None,
                     method_ids=None, parameter_ids=None, subsample=False,
                     filename_pattern=DEFAULT_EXPECTED_TABLE_PATTERN, size=10,
                     per_seq_precision=False, exclude=['other'], backup=True,
                     force=False, append=False,
                     enable_biom_cache=True, biom_cache_max_entries=None):
    '''Load observed and expected observations from tax-credit, compute
        precision, recall, F-measure, and correlations, and return results
        as dataframe.

        results_dirs: list of directories containing precomputed taxonomy
            assignment results to evaluate. Must be in format:
                results_dirs/<dataset name>/
                    <reference name>/<method>/<parameters>/
        expected_results_dir: directory containing expected composition data in
            the structure:
            expected_results_dir/<dataset name>/<reference name>/expected/
        results_fp: path to output file containing evaluation results summary.
        mock_dir: path
            Directory of mock community directiories containing mock feature
            tables without taxonomy.
        taxonomy_level_range: RANGE of taxonomic levels to evaluate.
        min_count: int
            Minimum abundance threshold for filtering taxonomic features.
        taxa_to_keep: list of taxonomies to retain, others are removed before
            evaluation.
        md_key: metadata key containing taxonomy metadata in observed taxonomy
            biom tables.
        dataset_ids: list
            dataset ids (mock community study ID) to process. Defaults to None
            (process all).
        reference_ids: list
            reference database data to process. Defaults to None (process all).
        method_ids: list
            methods to process. Defaults to None (process all).
        parameter_ids: list
            parameters to process. Defaults to None (process all).
        subsample: bool
            Randomly subsample results for test runs.
        size: int
            Size of subsample to take.
        exclude: list
            taxonomies to explicitly exclude from precision scoring.
        backup: bool
            Backup pre-existing results before overwriting? Will overwrite
            previous backups, and will only backup if force or append ==True.
        force: bool
            Overwrite pre-existing results_fp?
        append: bool
            Append new data to results_fp? Behavior of force and append will
            depend on whether the data in results_dirs have already been
            calculated in results_fp, and have interacting effects:

            if force=   append= Action
                True	True	Append new to results_fp; pre-existing results
                                are overwritten if they are requested by the
                                "results params": dataset_ids, reference_ids,
                                method_ids, parameter_ids. If these should be
                                excluded and results_fp should only include
                                results specifically requested, use force==True
                                and append==False.
                True	False	Overwrite results_fp with results requested by
                                "results params".
                False	True	Load results_fp and append new to results_fp;
                                pre-existing results are not overwritten even
                                if requested by "results params".
                False	False	Load results_fp. If "results params" are set,
                                the dataframe returned by this function is
                                automatically filtered to include only those
                                results.
        enable_biom_cache: bool
            If True (default), reuse mounted BIOM tables and feature tables
            across methods/parameters that share the same files (Phase 3).
        biom_cache_max_entries: int or None
            Optional LRU cap for the in-memory cache; ``None`` means unbounded
            for the duration of ``compute_mock_results``.

    '''
    from tax_credit.eval_framework import (
        get_expected_tables_lookup,
        seek_results,
    )

    results = seek_results(
        results_dirs, dataset_ids, reference_ids, method_ids, parameter_ids)

    if subsample is True:
        shuffle(results)
        results = results[:size]

    expected_tables = get_expected_tables_lookup(
        expected_results_dir, filename_pattern=filename_pattern)

    if not exists(results_fp) or force:
        if exists(results_fp) and append and (
                dataset_ids or reference_ids or method_ids or parameter_ids):
            old_results = pd.read_csv(results_fp, sep='\t', index_col=0)
            old_results = _filter_mock_results(
                old_results, dataset_ids, reference_ids, method_ids,
                parameter_ids)
        else:
            old_results = None

        mock_results = compute_mock_results(
            results, expected_tables, results_fp, mock_dir,
            taxonomy_level_range, min_count=min_count,
            taxa_to_keep=taxa_to_keep, md_key=md_key,
            per_seq_precision=per_seq_precision, exclude=exclude,
            enable_biom_cache=enable_biom_cache,
            biom_cache_max_entries=biom_cache_max_entries)

        if old_results is not None:
            mock_results = pd.concat([mock_results, old_results])
        _write_mock_results(mock_results, results_fp, backup)

    else:
        print("{0} already exists.".format(results_fp))
        print("Reading in pre-computed evaluation results.")
        print("To overwrite, set force=True")
        mock_results = pd.read_csv(results_fp, sep='\t', index_col=0)

        if append:
            results = _pending_results_for_append(mock_results, results)
            print("append==True and force==False")
            print(len(results), "new results have been appended to results.")
            if len(results) >= 1:
                new_mock_results = compute_mock_results(
                    results, expected_tables, results_fp, mock_dir,
                    taxonomy_level_range, min_count=min_count,
                    taxa_to_keep=taxa_to_keep, md_key=md_key,
                    per_seq_precision=per_seq_precision, exclude=exclude,
                    enable_biom_cache=enable_biom_cache,
                    biom_cache_max_entries=biom_cache_max_entries)
                mock_results = pd.concat([mock_results, new_mock_results])
                _write_mock_results(mock_results, results_fp, backup)

        elif dataset_ids or reference_ids or method_ids or parameter_ids:
            print("Results have been filtered to only include datasets or "
                  "reference databases or methods or parameters that are "
                  "explicitly set by results params. To disable this "
                  "function and load all results, set dataset_ids and "
                  "reference_ids and method_ids and parameter_ids to None.")
            mock_results = _filter_mock_results(
                mock_results, dataset_ids, reference_ids, method_ids,
                parameter_ids)

    return mock_results


def compute_mock_results(result_tables, expected_table_lookup, results_fp,
                         mock_dir, taxonomy_level_range=range(2, 7),
                         min_count=0,
                         taxa_to_keep=None, md_key='taxonomy',
                         per_seq_precision=False, exclude=None,
                         enable_biom_cache=True, biom_cache_max_entries=None):
    """Compute precision, recall, and f-measure for result_tables at each
    taxonomy level in *taxonomy_level_range*.

        result_tables: 2d list of tables to be compared to expected tables,
         where the data in the inner list is:
          [dataset_id, reference_database_id, method_id,
           parameter_combination_id, table_fp]
        expected_table_lookup: 2d dict of dataset_id, reference_db_id to BIOM
         table filepath, for the expected result tables
        taxonomy_level_range: range of levels to compute results
        results_fp: path to output file containing evaluation results summary
        mock_dir: path
            Directory of mock community directories that contain feature tables
            without taxonomy.
        per_seq_precision: bool
            Compute per-sequence precision/recall scores from expected
            taxonomy assignments?
        exclude: list
            taxonomies to explicitly exclude from precision scoring.
        enable_biom_cache: bool
            Reuse collapsed BIOM tables and raw feature tables when the same
            paths and parameters appear multiple times (default True).
        biom_cache_max_entries: int or None
            Optional LRU limit; ``None`` keeps all entries for this run.
    """
    from tax_credit.biom_cache import (
        BiomTableCache,
        NO_CACHE,
        feature_table_cache_key,
        mount_observations_cache_key,
    )
    from tax_credit.eval_framework import (
        compute_taxon_accuracy,
        mount_observations,
        per_sequence_precision,
    )

    if enable_biom_cache:
        cache = BiomTableCache(max_entries=biom_cache_max_entries)
    else:
        cache = NO_CACHE

    results = []
    for dataset_id, ref_id, method, params, actual_table_fp in result_tables:

        try:
            expected_table_fp = expected_table_lookup[dataset_id][ref_id]
        except KeyError:
            raise KeyError("Can't find expected table for \
                            ({0}, {1}).".format(dataset_id, ref_id))

        feature_table_fp = join(mock_dir, dataset_id, FEATURE_TABLE_BIOM)

        def _load_feature_table():
            try:
                return load_table(feature_table_fp)
            except ValueError:
                raise ValueError(
                    "Couldn't parse BIOM table: {0}".format(feature_table_fp))

        feature_table = cache.get_or_put(
            feature_table_cache_key(feature_table_fp), _load_feature_table)

        for taxonomy_level in taxonomy_level_range:
            exp_key = mount_observations_cache_key(
                expected_table_fp, 0, taxonomy_level, taxa_to_keep,
                'taxonomy', False)

            def _mount_expected(tl=taxonomy_level):
                return mount_observations(
                    expected_table_fp, min_count=0, taxonomy_level=tl,
                    taxa_to_keep=taxa_to_keep, filter_obs=False)

            expected_table = cache.get_or_put(exp_key, _mount_expected)

            act_key = mount_observations_cache_key(
                actual_table_fp, min_count, taxonomy_level, taxa_to_keep,
                md_key, True)

            def _mount_actual(tl=taxonomy_level):
                return mount_observations(
                    actual_table_fp, min_count=min_count,
                    taxonomy_level=tl, taxa_to_keep=taxa_to_keep,
                    md_key=md_key)

            actual_table = cache.get_or_put(act_key, _mount_actual)

            for sample_id in actual_table.ids(axis="sample"):
                try:
                    accuracy, detection = compute_taxon_accuracy(
                        actual_table, expected_table,
                        actual_sample_id=sample_id,
                        expected_sample_id=sample_id)
                except ZeroDivisionError:
                    accuracy, detection = -1., -1.

                if per_seq_precision and exists(join(
                        dirname(expected_table_fp), TRUEISH_TAXONOMIES_TSV)):
                    p, r, f = per_sequence_precision(
                        expected_table_fp, actual_table_fp, feature_table,
                        sample_id, taxonomy_level, exclude=exclude)
                else:
                    p, r, f = -1., -1., -1.

                results.append((dataset_id, taxonomy_level, sample_id,
                                ref_id, method, params, p, r, f, accuracy,
                                detection))

    result = pd.DataFrame(results, columns=["Dataset", "Level", "SampleID",
                                            "Reference", "Method",
                                            "Parameters", "Precision",
                                            "Recall", "F-measure",
                                            "Taxon Accuracy Rate",
                                            "Taxon Detection Rate"])
    return result
