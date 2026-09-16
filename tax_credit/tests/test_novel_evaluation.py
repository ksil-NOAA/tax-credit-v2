#!/usr/bin/env python

import unittest
from os import makedirs
from os.path import join
from tempfile import TemporaryDirectory

import pandas as pd

from tax_credit import framework_functions
from tax_credit import novel_evaluation
from tax_credit.novel_evaluation import (
    _truncate_taxonomy_at_level,
    extract_per_level_classification_ratios,
    extract_per_level_classification_ratios_by_fold,
    per_level_classification_ratios_from_log,
    select_best_runs,
)

# reference lineage with an unknown (NA) order above assigned lower ranks
NA_TAXON = "Eukaryota;Chordata;Actinopteri;NA;Lutjanidae;Lutjanus;Lutjanus griseus"


class NovelEvaluationModuleTests(unittest.TestCase):
    """Phase 2: novel/CV assignment evaluation lives in ``novel_evaluation``."""

    def test_reexports_match_module(self):
        self.assertIs(
            framework_functions.novel_taxa_classification_evaluation,
            novel_evaluation.novel_taxa_classification_evaluation,
        )
        self.assertIs(
            framework_functions.extract_per_level_accuracy,
            novel_evaluation.extract_per_level_accuracy,
        )

    def test_truncate_taxonomy_at_level_keeps_internal_na(self):
        self.assertEqual(
            _truncate_taxonomy_at_level(NA_TAXON, 4),
            "Eukaryota;Chordata;Actinopteri;NA;Lutjanidae",
        )
        # truncating onto the NA rank leaves it trailing, so it is dropped
        self.assertEqual(
            _truncate_taxonomy_at_level(NA_TAXON, 3),
            "Eukaryota;Chordata;Actinopteri",
        )
        self.assertEqual(
            _truncate_taxonomy_at_level("A;B;C;D;E;F;NA", 6), "A;B;C;D;E;F"
        )

    def test_per_level_ratios_with_internal_na(self):
        # Wrong genus under an NA order: the error belongs at level 5 (genus).
        # Dropping the NA rank would shift it up to level 4 (family).
        observed = "Eukaryota;Chordata;Actinopteri;NA;Lutjanidae;Ocyurus;Ocyurus chrysurus"
        header = ["dataset", "level", "iteration", "method", "parameters",
                  "observed_taxonomy", "expected_taxonomy", "result",
                  "mismatch_level"]
        row = ["db1", "6", "0", "bt2-blca", "pi0.8", observed, NA_TAXON,
               "misclassification", "5"]
        with TemporaryDirectory() as tmp:
            log_fp = join(tmp, "classification_accuracy_log.tsv")
            with open(log_fp, "w") as f:
                f.write("\t".join(header) + "\n" + "\t".join(row) + "\n")
            rows, _ = per_level_classification_ratios_from_log(log_fp)

        by_level = {r["level"]: r for r in rows}
        self.assertEqual(by_level[3]["match_ratio"], 1.0)
        self.assertEqual(by_level[4]["match_ratio"], 1.0)
        self.assertEqual(by_level[5]["misclassification_ratio"], 1.0)
        self.assertEqual(by_level[6]["misclassification_ratio"], 1.0)

    def test_per_level_match_ratio_is_recall(self):
        # Two matches: one full 7-rank taxon (mismatch level 7) and one against a
        # genus-only reference (mismatch level 6). Summing mismatch_level_list
        # would give 0.5 at species; recall is 1.0.
        df = pd.DataFrame([{
            "Dataset": "db1", "level": 6, "iteration": "0",
            "Method": "m", "Parameters": "p",
            "mismatch_level_list": "[0, 0, 0, 0, 0, 0, 1, 1]",
            "Precision": "[0, 1, 1, 1, 1, 1, 1]",
            "Recall": "[0, 1, 1, 1, 1, 1, 1]",
            "F-measure": "[0, 1, 1, 1, 1, 1, 1]",
        }])
        pla = novel_evaluation.extract_per_level_accuracy(df)
        self.assertEqual(list(pla["match_ratio"]), [1.0] * 6)

    def test_classification_ratios_keep_novel_levels(self):
        header = ["dataset", "level", "iteration", "method", "parameters",
                  "observed_taxonomy", "expected_taxonomy", "result",
                  "mismatch_level"]
        # L5 read matches; L6 read is misclassified at class
        logs = {
            "db1-L5-iter0": ("A;B;C;D;E", "A;B;C;D;E", "match", "5"),
            "db1-L6-iter0": ("A;B;X;D;E;F", "A;B;C;D;E;F", "misclassification", "2"),
        }
        with TemporaryDirectory() as tmp:
            dirs = []
            for dataset_id, (obs, exp, result, mismatch) in logs.items():
                results_dir = join(tmp, dataset_id, "method1", "param1")
                makedirs(results_dir)
                row = ["db1", dataset_id[5], "0", "method1", "param1",
                       obs, exp, result, mismatch]
                with open(join(results_dir, "classification_accuracy_log.tsv"), "w") as f:
                    f.write("\t".join(header) + "\n" + "\t".join(row) + "\n")
                dirs.append(results_dir)

            by_fold = extract_per_level_classification_ratios_by_fold(dirs)
            averaged = extract_per_level_classification_ratios(dirs)

        self.assertEqual(set(by_fold["iteration"]), {"0"})
        self.assertEqual(sorted(averaged["novel_level"].unique()), [5, 6])
        class_level = averaged[averaged["level"] == 2].set_index("novel_level")
        self.assertEqual(class_level.loc[5, "match_ratio"], 1.0)
        self.assertEqual(class_level.loc[6, "misclassification_ratio"], 1.0)

    def test_select_best_runs(self):
        df = pd.DataFrame([
            {"Dataset": "db1", "Method": "a", "Parameters": "p1",
             "F-measure": 0.5, "misclassification_ratio": 0.1},
            {"Dataset": "db1", "Method": "a", "Parameters": "p1",
             "F-measure": 0.7, "misclassification_ratio": 0.3},
            {"Dataset": "db1", "Method": "b", "Parameters": "p1",
             "F-measure": 0.4, "misclassification_ratio": 0.2},
            {"Dataset": "db1", "Method": "b", "Parameters": "p2",
             "F-measure": 0.65, "misclassification_ratio": 0.2},
            {"Dataset": "db2", "Method": "b", "Parameters": "p1",
             "F-measure": 0.9, "misclassification_ratio": 0.0},
        ])
        best = select_best_runs(df, ["F-measure", "misclassification_ratio"])
        db1 = best[best["Dataset"] == "db1"].set_index("metric")

        # highest mean F-measure: b/p2 (0.65) beats a/p1 averaged over 2 folds (0.6)
        self.assertEqual(tuple(db1.loc["F-measure", ["Method", "Parameters"]]), ("b", "p2"))
        self.assertEqual(db1.loc["F-measure", "direction"], "highest")
        self.assertEqual(db1.loc["F-measure", "n_tied"], 1)
        # lowest misclassification: all three runs average 0.2; first sorted run wins
        self.assertEqual(tuple(db1.loc["misclassification_ratio", ["Method", "Parameters"]]), ("a", "p1"))
        self.assertEqual(db1.loc["misclassification_ratio", "direction"], "lowest")
        self.assertEqual(db1.loc["misclassification_ratio", "n_folds"], 2)
        self.assertEqual(db1.loc["misclassification_ratio", "n_tied"], 3)
        self.assertEqual(len(best[best["Dataset"] == "db2"]), 2)


if __name__ == '__main__':
    unittest.main()
