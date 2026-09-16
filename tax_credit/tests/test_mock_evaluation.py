#!/usr/bin/env python

import unittest

from tax_credit.eval_framework import compute_mock_results, evaluate_results
from tax_credit import mock_evaluation


class MockEvaluationModuleTests(unittest.TestCase):
    """Phase 2: mock evaluation logic lives in ``mock_evaluation``; eval keeps API."""

    def test_reexports_match_module(self):
        self.assertIs(evaluate_results, mock_evaluation.evaluate_results)
        self.assertIs(compute_mock_results, mock_evaluation.compute_mock_results)

    def test_pending_results_for_append(self):
        import pandas as pd

        mock_results = pd.DataFrame({
            'Dataset': ['a'],
            'Reference': ['r'],
            'Method': ['m'],
            'Parameters': ['p'],
        })
        results = [
            ('a', 'r', 'm', 'p', '/x'),
            ('b', 'r', 'm', 'p', '/y'),
        ]
        pending = mock_evaluation._pending_results_for_append(mock_results, results)
        self.assertEqual(pending, [results[1]])

    def test_biom_cache_reduces_mount_observations_for_shared_expected(self):
        """Phase 3: same dataset/reference reuses one collapsed expected table."""
        from unittest.mock import MagicMock, patch

        def biom_like_table(*_a, **_k):
            t = MagicMock()
            t.ids = MagicMock(
                side_effect=lambda axis=None: (
                    ['s1'] if axis == 'sample' else ['o1']))
            t.get_value_by_ids = MagicMock(return_value=1.0)
            return t

        result_tables = [
            ('ds', 'ref', 'm1', 'p1', '/path/a1.biom'),
            ('ds', 'ref', 'm2', 'p2', '/path/a2.biom'),
        ]
        expected_table_lookup = {'ds': {'ref': '/path/exp.biom'}}
        tlr = range(2, 4)

        def count_mount_calls(enable_cache):
            with patch(
                'tax_credit.eval_framework.mount_observations',
                side_effect=biom_like_table,
            ) as m_mount:
                with patch(
                    'tax_credit.mock_evaluation.load_table',
                    side_effect=biom_like_table,
                ):
                    with patch(
                        'tax_credit.eval_framework.compute_taxon_accuracy',
                        return_value=(0.5, 0.5),
                    ):
                        with patch(
                            'tax_credit.eval_framework.per_sequence_precision',
                            return_value=(-1.0, -1.0, -1.0),
                        ):
                            mock_evaluation.compute_mock_results(
                                result_tables,
                                expected_table_lookup,
                                'out.tsv',
                                '/mock',
                                taxonomy_level_range=tlr,
                                enable_biom_cache=enable_cache,
                            )
                return m_mount.call_count

        n_cached = count_mount_calls(True)
        n_uncached = count_mount_calls(False)
        self.assertLess(n_cached, n_uncached)
        self.assertEqual(n_uncached, 8)
        self.assertEqual(n_cached, 6)


if __name__ == '__main__':
    unittest.main()
