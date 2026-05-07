#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2014--, tax-credit development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


from unittest import TestCase, main
import json
from biom import Table
from tax_credit import mock_denoise
from tax_credit import mock_transport
from tax_credit import mockrobiota_extract
from tax_credit.process_mocks import amend_biom_taxonomy_ids
import tax_credit.process_mocks as process_mocks


class EvalFrameworkTests(TestCase):

    def test_amend_biom_taxonomy_ids(self):
        self.assertEqual(set(amend_biom_taxonomy_ids(self.table1,
                             clean_obs_ids=False).ids(axis='observation')),
                         {'k__Archaea;p__;c__;o__;f__;g__;s__',
                          'k__Bacteria;p__;c__;o__;f__;g__;s__',
                          'k__[Fungi];p__;c__;o__;f__;g__;s__'})
        # This also tests clean_taxonomy_ids()
        self.assertEqual(set(amend_biom_taxonomy_ids(self.table1,
                             clean_obs_ids=True).ids(axis='observation')),
                         {'k__Archaea;p__;c__;o__;f__;g__;s__',
                          'k__Bacteria;p__;c__;o__;f__;g__;s__',
                          'k__Fungi;p__;c__;o__;f__;g__;s__'})

    def setUp(self):
        _table1 = """{"id": "None",
                      "format": "Biological Observation Matrix 1.0.0",
                      "format_url": "http:\/\/biom-format.org",
                      "type": "OTU table",
                      "generated_by": "greg",
                      "date": "2013-08-22T13:10:23.907145",
                      "matrix_type": "sparse",
                      "matrix_element_type": "float",
                      "shape": [
                        3,
                        4
                      ],
                      "data": [
                        [
                          0,
                          0,
                          1
                        ],
                        [
                          0,
                          1,
                          2
                        ],
                        [
                          0,
                          2,
                          3
                        ],
                        [
                          0,
                          3,
                          4
                        ],
                        [
                          1,
                          0,
                          2
                        ],
                        [
                          1,
                          1,
                          0
                        ],
                        [
                          1,
                          2,
                          7
                        ],
                        [
                          1,
                          3,
                          8
                        ],
                        [
                          2,
                          0,
                          9
                        ],
                        [
                          2,
                          1,
                          10
                        ],
                        [
                          2,
                          2,
                          11
                        ],
                        [
                          2,
                          3,
                          12
                        ]
                      ],
                      "rows": [
                        {
                          "id": "k__Bacteria",
                          "metadata": {
                            "domain": "Bacteria"
                          }
                        },
                        {
                          "id": "k__Archaea",
                          "metadata": {
                            "domain": "Archaea"
                          }
                        },
                        {
                          "id": "k__[Fungi]",
                          "metadata": {
                            "domain": "[Fungi]"
                          }
                        }
                      ],
                      "columns": [
                        {
                          "id": "s1",
                          "metadata": {
                            "country": "Peru",
                            "pH": 4.2
                          }
                        },
                        {
                          "id": "s2",
                          "metadata": {
                            "country": "Peru",
                            "pH": 5.2
                          }
                        },
                        {
                          "id": "s3",
                          "metadata": {
                            "country": "Peru",
                            "pH": 5
                          }
                        },
                        {
                          "id": "s4",
                          "metadata": {
                            "country": "Peru",
                            "pH": 4.9
                          }
                        }
                      ]
                    }"""
        # table 1
        # OTU ID	   s1	s2	s3	s4
        # k__Archaea    1.0 2.0 3.0 4.0
        # k__Bacteria    2.0 0.0 7.0 8.0
        # k__[Fungi]    9.0 10.0    11.0    12.0

        self.table1 = Table.from_json(json.loads(_table1))


class ProcessMocksPhase5FacadeTests(TestCase):
    """Phase 5: ``process_mocks`` re-exports match submodule implementations."""

    def test_mockrobiota_extract_reexports(self):
        self.assertIs(
            process_mocks.extract_mockrobiota_dataset_metadata,
            mockrobiota_extract.extract_mockrobiota_dataset_metadata,
        )
        self.assertIs(
            process_mocks.extract_mockrobiota_data,
            mockrobiota_extract.extract_mockrobiota_data,
        )
        self.assertIs(
            process_mocks.amend_biom_taxonomy_ids,
            mockrobiota_extract.amend_biom_taxonomy_ids,
        )
        self.assertIs(
            process_mocks.clean_taxonomy_ids,
            mockrobiota_extract.clean_taxonomy_ids,
        )

    def test_mock_denoise_reexports(self):
        self.assertIs(process_mocks.batch_demux, mock_denoise.batch_demux)
        self.assertIs(
            process_mocks.demux_to_plot_quality,
            mock_denoise.demux_to_plot_quality,
        )
        self.assertIs(process_mocks.visualize_qual, mock_denoise.visualize_qual)
        self.assertIs(
            process_mocks.load_demux_seqs, mock_denoise.load_demux_seqs
        )
        self.assertIs(
            process_mocks.denoise_to_phylogeny,
            mock_denoise.denoise_to_phylogeny,
        )
        self.assertIs(
            process_mocks.denoise_to_feature_table,
            mock_denoise.denoise_to_feature_table,
        )
        self.assertIs(process_mocks.seqs_to_tree, mock_denoise.seqs_to_tree)

    def test_mock_transport_reexports(self):
        self.assertIs(
            process_mocks.transport_to_repo,
            mock_transport.transport_to_repo,
        )


if __name__ == "__main__":
    main()
