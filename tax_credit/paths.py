#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2014--, tax-credit development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""Filesystem layout and filename conventions for tax-credit analyses.

This module centralizes glob patterns and path parsing so directory depth and
basename contracts stay consistent across eval, simulation, and utilities.
"""

from __future__ import annotations

from pathlib import Path
from typing import NamedTuple

# --- Subdirectory segment under a dataset / reference pair -----------------
EXPECTED_SUBDIR = "expected"

# --- Default BIOM / table naming -------------------------------------------
DEFAULT_MOCK_RESULT_TABLE_PATTERN = "table*biom"
DEFAULT_EXPECTED_TABLE_PATTERN = "table.L{0}-taxa.biom"

FEATURE_TABLE_BIOM = "feature_table.biom"
MERGED_TABLE_BIOM = "merged_table.biom"

# --- Per-sequence and taxonomy sidecar files --------------------------------
TRUEISH_TAXONOMIES_TSV = "trueish-taxonomies.tsv"
REP_SEQS_TAX_ASSIGNMENTS_TXT = "rep_seqs_tax_assignments.txt"
TAXONOMY_TSV = "taxonomy.tsv"

# --- Novel / CV assignment evaluation --------------------------------------
QUERY_TAX_ASSIGNMENTS_TXT = "query_tax_assignments.txt"
QUERY_TAXA_TSV = "query_taxa.tsv"
CLASSIFICATION_ACCURACY_LOG_TSV = "classification_accuracy_log.tsv"

# --- CV / novel fold contents (FASTA + taxonomy TSV) ------------------------
REF_SEQS_FASTA = "ref_seqs.fasta"
QUERY_FASTA = "query.fasta"
REF_TAXA_TSV = "ref_taxa.tsv"


class MockResultTableParts(NamedTuple):
    """Segments for ``<root>/<dataset>/<reference>/<method>/<params>/<file>``."""

    dataset_id: str
    reference_id: str
    method_id: str
    parameter_id: str


class ExpectedTableParts(NamedTuple):
    """Segments for ``<root>/<dataset>/<reference>/expected/<file>``."""

    dataset_id: str
    reference_id: str


class AssignmentResultParts(NamedTuple):
    """Segments for ``<...>/<dataset_id>/<method_id>/<params_id>`` (novel/CV results)."""

    dataset_id: str
    method_id: str
    params_id: str


def mock_observed_tables_glob(
    start_dir: str, filename_pattern: str = DEFAULT_MOCK_RESULT_TABLE_PATTERN
) -> str:
    """Glob pattern for observed mock/CV-style BIOM tables under *start_dir*."""
    return str(Path(start_dir) / "*" / "*" / "*" / "*" / filename_pattern)


def expected_tables_glob(start_dir: str, filename: str) -> str:
    """Glob pattern for expected BIOM tables under *start_dir*."""
    return str(Path(start_dir) / "*" / "*" / EXPECTED_SUBDIR / filename)


def parse_mock_result_table_path(table_fp: str) -> MockResultTableParts:
    """Parse *table_fp* into dataset / reference / method / parameter ids."""
    parts = Path(table_fp).parts
    if len(parts) < 5:
        raise ValueError(
            "Observed result table path must end with "
            "<dataset>/<reference>/<method>/<parameters>/<file>; "
            "got {!r}".format(table_fp)
        )
    dataset_id, reference_id, method_id, parameter_id = (
        parts[-5],
        parts[-4],
        parts[-3],
        parts[-2],
    )
    return MockResultTableParts(
        dataset_id, reference_id, method_id, parameter_id
    )


def parse_expected_table_path(table_fp: str) -> ExpectedTableParts:
    """Parse *table_fp* into dataset and reference ids (…/expected/<file>)."""
    parts = Path(table_fp).parts
    if len(parts) < 4 or parts[-2] != EXPECTED_SUBDIR:
        raise ValueError(
            "Expected table path must end with "
            "<dataset>/<reference>/{}/<file>; got {!r}".format(
                EXPECTED_SUBDIR, table_fp
            )
        )
    return ExpectedTableParts(parts[-4], parts[-3])


def parse_assignment_results_dir(results_dir: str) -> AssignmentResultParts:
    """Parse a novel- or CV-assignment results directory (three trailing segments)."""
    parts = Path(results_dir).parts
    if len(parts) < 3:
        raise ValueError(
            "Results directory path must end with "
            "<dataset_id>/<method_id>/<params_id>; got {!r}".format(results_dir)
        )
    return AssignmentResultParts(parts[-3], parts[-2], parts[-1])


def parse_taxonomy_map_path_to_dataset_id(taxonomy_map_fp: str) -> str:
    """Infer *dataset_id* from ``…/<dataset>/<reference>/<method>/<param>/<map>``."""
    parts = Path(taxonomy_map_fp).parts
    if len(parts) < 5:
        raise ValueError(
            "Taxonomy map path must have at least five trailing segments "
            "(dataset/reference/method/parameter/file); got {!r}".format(
                taxonomy_map_fp
            )
        )
    return parts[-5]


def parse_result_leaf_dir_to_parts(method_dir: str) -> MockResultTableParts:
    """Parse ``…/<dataset>/<reference>/<method>/<param>`` (no filename segment)."""
    parts = Path(method_dir).parts
    if len(parts) < 4:
        raise ValueError(
            "Result leaf directory must end with "
            "<dataset>/<reference>/<method>/<parameters>; got {!r}".format(
                method_dir
            )
        )
    return MockResultTableParts(
        parts[-4], parts[-3], parts[-2], parts[-1]
    )
