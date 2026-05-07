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
from os.path import join

import pandas as pd

from tax_credit.paths import (
    CLASSIFICATION_ACCURACY_LOG_TSV,
    QUERY_TAX_ASSIGNMENTS_TXT,
    QUERY_TAXA_TSV,
    parse_assignment_results_dir,
)
from tax_credit.simulation_names import parse_cv_dataset_id, parse_novel_dataset_id
from tax_credit.taxa_manipulator import export_list_to_file


def novel_taxa_classification_evaluation(results_dirs, expected_results_dir,
                                         summary_fp, test_type='novel-taxa'):
    '''Input glob of novel taxa results, receive a summary of accuracy results.
    results_dirs = list or glob of novel taxa observed results in format:
                    precomputed_results_dir/dataset_id/method_id/params_id/
    expected_results_dir = directory containing expected novel-taxa results in
                    format:
                    expected_results_dir/dataset_id/method_id/params_id/
    summary_fp = filepath to contain summary of results
    test_type = one of 'novel-taxa', 'cross-validated',
        'cross-validated-trad'

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
        else:
            raise ValueError(
                'test_type must be "novel-taxa", "cross-validated", or '
                '"cross-validated-trad"')

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
                    linecount = sum(col)
                    col_names.append("match_ratio")
                    cumulative_mismatches = sum(col[0:level+1])
                    if cumulative_mismatches < linecount:
                        score = (linecount - cumulative_mismatches) / linecount
                    else:
                        score = col[0]
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
