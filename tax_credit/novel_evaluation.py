#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2016--, tax-credit development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""Novel-taxa and cross-validated assignment evaluation from query tax files.

Consumes paths laid out under ``tax_credit.paths`` and dataset ids from
``tax_credit.simulation_names``; PRF and taxonomy helpers stay in
``framework_functions``.
"""

from collections import Counter
from os.path import exists, join

import pandas as pd

from tax_credit.paths import (
    CLASSIFICATION_ACCURACY_LOG_TSV,
    QUERY_TAX_ASSIGNMENTS_TXT,
    QUERY_TAXA_TSV,
    parse_assignment_results_dir,
)
from tax_credit.simulation_names import (
    parse_cv_dataset_id,
    parse_novel_dataset_id,
    parse_self_validated_dataset_id,
)
from tax_credit.taxa_manipulator import export_list_to_file, normalize_taxon


def novel_taxa_classification_evaluation(results_dirs, expected_results_dir,
                                         summary_fp, test_type='novel-taxa'):
    '''Input glob of novel taxa results, receive a summary of accuracy results.
    results_dirs = list or glob of novel taxa observed results in format:
                    precomputed_results_dir/dataset_id/method_id/params_id/
                    For the usual four-level sweep tree under a single root, use
                    ``tax_credit.paths.list_assignment_result_dirs`` instead of
                    hand-written globs.
    expected_results_dir = directory containing expected novel-taxa results in
                    format:
                    expected_results_dir/dataset_id/method_id/params_id/
    summary_fp = filepath to contain summary of results
    test_type = one of 'novel-taxa', 'cross-validated',
        'cross-validated-trad', or 'self-validated'

    Returns results as df, in addition to printing summary_fp
    '''
    from tax_credit.framework_functions import (
        compute_prf,
        count_records,
        evaluate_classification,
        find_last_common_ancestor,
        load_prf,
    )

    results = []

    for results_dir in results_dirs:
        res_parts = parse_assignment_results_dir(results_dir)
        dataset_id, method_id, params_id = res_parts

        if test_type == 'novel-taxa':
            novel_parts = parse_novel_dataset_id(dataset_id)
            index, level, iteration = (
                novel_parts.database,
                novel_parts.level,
                novel_parts.iteration,
            )
        elif test_type in ('cross-validated', 'cross-validated-trad'):
            cv_parts = parse_cv_dataset_id(dataset_id)
            index, iteration = cv_parts.database, cv_parts.iteration
            level = 6
        elif test_type == 'self-validated':
            sv_parts = parse_self_validated_dataset_id(dataset_id)
            index, iteration = sv_parts.database, sv_parts.iteration
            level = 6
        else:
            raise ValueError(
                'test_type must be "novel-taxa", "cross-validated", '
                '"cross-validated-trad", or "self-validated"')

        obs_fp = join(results_dir, QUERY_TAX_ASSIGNMENTS_TXT)
        exp_fp = join(expected_results_dir, dataset_id, QUERY_TAXA_TSV)
        exp_taxa, obs_taxa = load_prf(obs_fp, exp_fp)

        p, r, f = compute_prf(exp_taxa, obs_taxa, test_type=test_type)

        mismatch_level_list = [0] * 8
        log = ['dataset\tlevel\titeration\tmethod\tparameters\
               \tobserved_taxonomy\texpected_taxonomy\tresult\tmismatch_level\
               \tPrecision\tRecall\tF-measure']

        record_counter = Counter()
        for obs, exp in zip(obs_taxa, exp_taxa):
            mismatch_level = find_last_common_ancestor(obs, exp)
            mismatch_level_list[mismatch_level] += 1

            result = evaluate_classification(obs, exp)

            record_counter.update({'line_count': 1})
            record_counter.update({result: 1})
            log.append('\t'.join(map(str, [index, level, iteration,
                                           method_id, params_id,
                                           obs, exp, result,
                                           mismatch_level, p, r, f])))

        log_fp = join(results_dir, CLASSIFICATION_ACCURACY_LOG_TSV)
        export_list_to_file(log, log_fp)

        match_ratio = count_records(record_counter, 'match', 'line_count')
        overclass = count_records(record_counter, 'overclassification',
                                  'line_count')
        underclass = count_records(record_counter, 'underclassification',
                                   'line_count')
        misclass = count_records(record_counter, 'misclassification',
                                 'line_count')

        results.append((index, level, iteration, method_id, params_id,
                        match_ratio, overclass, underclass, misclass,
                        mismatch_level_list, p, r, f))

    result = pd.DataFrame(results, columns=["Dataset", "level", "iteration",
                                            "Method", "Parameters",
                                            "match_ratio",
                                            "overclassification_ratio",
                                            "underclassification_ratio",
                                            "misclassification_ratio",
                                            "mismatch_level_list", "Precision",
                                            "Recall", "F-measure"])
    result.to_csv(summary_fp)
    return result


def extract_per_level_accuracy(df, columns=['Precision', 'Recall', 'F-measure',
                                            'mismatch_level_list']):
    '''Generate new pandas dataframe, containing match ratios for taxonomic
    assignments at each taxonomic level. Extracts mismatch_level_list from a
    dataframe and splits this into separate df entries for plotting.

    df: dataframe
        pandas dataframe
    column: list
        column names containing mismatch_level_list or other lists to be
        separated into multiple dataframe entries for plotting.

        mismatch_level_list reports mismatches at each level of taxonomic
        assignment (8 levels).

        Currently levels  are hardcoded, but could be adjusted
        below in lines:
            for level in range(1, 7):
    '''
    results = []

    for index, data in df.iterrows():
        for level in range(1, 7):
            level_results = []
            col_names = []
            for column in columns:
                if isinstance(data[column], str):
                    col = list(map(float, data[column].strip('[]').split(',')))
                else:
                    col = data[column]
                if column == 'mismatch_level_list':
                    # Match ratio at a level is recall at that level. Summing
                    # mismatch_level_list counted matches against references
                    # with trailing NA ranks as mismatches.
                    col_names.append("match_ratio")
                    recall = data['Recall']
                    if isinstance(recall, str):
                        recall = list(map(float, recall.strip('[]').split(',')))
                    score = recall[level]
                else:
                    score = col[level]
                    col_names.append(column)

                level_results.append(score)

            results.append((data['Dataset'], level, data['iteration'],
                            data['Method'], data['Parameters'],
                            *[s for s in level_results]))

    result = pd.DataFrame(results, columns=["Dataset", "level", "iteration",
                                            "Method", "Parameters",
                                            *[s for s in col_names]])
    return result


CLASSIFICATION_RATIO_COLS = [
    "match_ratio",
    "misclassification_ratio",
    "overclassification_ratio",
    "underclassification_ratio",
]


def _truncate_taxonomy_at_level(taxon, level):
    """Return taxonomy truncated through rank ``level`` (1–6, matching PRF indices).

    Internal NA ranks keep their position; trailing NA ranks are dropped.
    """
    return normalize_taxon(";".join(taxon.split(";")[:level + 1]))


def _dataset_parts_from_dataset_id(dataset_id):
    """Return ``(database, novel_level, iteration)`` for an assignment dataset id.

    ``novel_level`` is ``None`` for cross-validated and self-validated ids.
    """
    try:
        parts = parse_novel_dataset_id(dataset_id)
        return parts.database, parts.level, str(parts.iteration)
    except (ValueError, TypeError):
        pass
    for parser in (parse_cv_dataset_id, parse_self_validated_dataset_id):
        try:
            parts = parser(dataset_id)
            return parts.database, None, str(parts.iteration)
        except (ValueError, TypeError):
            continue
    return dataset_id.split("-", 1)[0], None, "0"


def per_level_classification_ratios_from_log(log_fp):
    """Compute classification ratios at ranks 1–6 from a per-read accuracy log."""
    from tax_credit.framework_functions import count_records, evaluate_classification

    counters = {level: Counter() for level in range(1, 7)}
    meta = {}

    with open(log_fp) as f:
        header = [col.strip() for col in f.readline().strip().split("\t")]
        col_idx = {name: idx for idx, name in enumerate(header)}
        obs_idx = col_idx["observed_taxonomy"]
        exp_idx = col_idx["expected_taxonomy"]
        for key in ("dataset", "method", "parameters"):
            if key in col_idx:
                meta[key] = None

        for line in f:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            obs = fields[obs_idx]
            exp = fields[exp_idx]
            if meta:
                for key in meta:
                    meta[key] = fields[col_idx[key]]
            for level in range(1, 7):
                obs_t = _truncate_taxonomy_at_level(obs, level)
                exp_t = _truncate_taxonomy_at_level(exp, level)
                result = evaluate_classification(obs_t, exp_t)
                counters[level].update({"line_count": 1, result: 1})

    rows = []
    for level in range(1, 7):
        counter = counters[level]
        rows.append({
            "level": level,
            "match_ratio": count_records(counter, "match", "line_count"),
            "overclassification_ratio": count_records(
                counter, "overclassification", "line_count"
            ),
            "underclassification_ratio": count_records(
                counter, "underclassification", "line_count"
            ),
            "misclassification_ratio": count_records(
                counter, "misclassification", "line_count"
            ),
        })
    return rows, meta


def extract_per_level_classification_ratios_by_fold(results_dirs):
    """Per-level classification ratios, one row per result directory and level.

    Reads ``classification_accuracy_log.tsv`` from each assignment result
    directory. Columns: ``Dataset``, ``novel_level`` (novel-taxa simulation
    level, ``<NA>`` otherwise), ``iteration``, ``Method``, ``Parameters``,
    ``level`` and ``CLASSIFICATION_RATIO_COLS``.
    """
    rows = []
    for results_dir in results_dirs:
        log_fp = join(results_dir, CLASSIFICATION_ACCURACY_LOG_TSV)
        if not exists(log_fp):
            continue
        dataset_id, method_id, params_id = parse_assignment_results_dir(results_dir)
        database, novel_level, iteration = _dataset_parts_from_dataset_id(dataset_id)
        level_rows, _ = per_level_classification_ratios_from_log(log_fp)
        for level_row in level_rows:
            rows.append({
                "Dataset": database,
                "novel_level": novel_level,
                "iteration": iteration,
                "Method": method_id,
                "Parameters": params_id,
                **level_row,
            })

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    df["novel_level"] = df["novel_level"].astype("Int64")
    return df


def extract_per_level_classification_ratios(results_dirs):
    """Per-level classification ratios averaged across iterations.

    Groups by Dataset, novel_level, Method, Parameters and level, so novel-taxa
    simulations at different novel levels are never averaged together.
    """
    df = extract_per_level_classification_ratios_by_fold(results_dirs)
    if df.empty:
        return df
    return df.groupby(
        ["Dataset", "novel_level", "Method", "Parameters", "level"],
        as_index=False,
        dropna=False,
    )[CLASSIFICATION_RATIO_COLS].mean()


# Error ratios, where the best run has the lowest value.
LOWER_IS_BETTER_METRICS = frozenset({
    "misclassification_ratio",
    "overclassification_ratio",
    "underclassification_ratio",
})


def select_best_runs(df, metrics, group_cols=("Dataset",),
                     run_cols=("Method", "Parameters"), tolerance=1e-9):
    """Pick the best method + parameter run per group for each metric.

    *df* holds one row per fold at a single taxonomic level (or novel level).
    Each run's metric is averaged across its rows; metrics in
    ``LOWER_IS_BETTER_METRICS`` are minimised and all others maximised. Ties
    within *tolerance* go to the first run sorted by *run_cols*.

    Returns one row per group and metric with the group columns, ``metric``,
    ``direction`` (``highest`` / ``lowest``), the run columns, ``value``
    (mean score), ``n_folds`` and ``n_tied`` (runs sharing the best value).
    """
    group_cols, run_cols, metrics = list(group_cols), list(run_cols), list(metrics)
    grouped = df.groupby(group_cols + run_cols, dropna=False)
    means = grouped[metrics].mean()
    means["n_folds"] = grouped.size()
    means = means.reset_index()

    rows = []
    for group_key, group in means.groupby(group_cols, dropna=False, sort=True):
        if not isinstance(group_key, tuple):
            group_key = (group_key,)
        for metric in metrics:
            scored = group.dropna(subset=[metric])
            if scored.empty:
                continue
            lower = metric in LOWER_IS_BETTER_METRICS
            best_value = scored[metric].min() if lower else scored[metric].max()
            tied = scored[(scored[metric] - best_value).abs() <= tolerance]
            best = tied.sort_values(run_cols).iloc[0]
            rows.append({
                **dict(zip(group_cols, group_key)),
                "metric": metric,
                "direction": "lowest" if lower else "highest",
                **{col: best[col] for col in run_cols},
                "value": float(best_value),
                "n_folds": int(best["n_folds"]),
                "n_tied": len(tied),
            })
    return pd.DataFrame(rows)
