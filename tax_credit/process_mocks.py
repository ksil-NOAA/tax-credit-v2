#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2016--, tax-credit development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""Mock community QIIME 2 pipeline helpers (Phase 5 facade).

Implementation is split into ``mockrobiota_extract``, ``mock_denoise``, and
``mock_transport``; this module re-exports the same public API as before.
"""

from tax_credit.mock_denoise import (
    batch_demux,
    demux_to_plot_quality,
    denoise_to_feature_table,
    denoise_to_phylogeny,
    load_demux_seqs,
    seqs_to_tree,
    visualize_qual,
)
from tax_credit.mock_transport import transport_to_repo
from tax_credit.mockrobiota_extract import (
    amend_biom_taxonomy_ids,
    clean_taxonomy_ids,
    extract_mockrobiota_data,
    extract_mockrobiota_dataset_metadata,
)

__all__ = [
    "amend_biom_taxonomy_ids",
    "batch_demux",
    "clean_taxonomy_ids",
    "demux_to_plot_quality",
    "denoise_to_feature_table",
    "denoise_to_phylogeny",
    "extract_mockrobiota_data",
    "extract_mockrobiota_dataset_metadata",
    "load_demux_seqs",
    "seqs_to_tree",
    "transport_to_repo",
    "visualize_qual",
]
