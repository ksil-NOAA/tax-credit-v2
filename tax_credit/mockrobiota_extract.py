#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2016--, tax-credit development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""Mockrobiota download, expected taxonomy copy, and BIOM taxonomy ID cleanup."""

from os import makedirs
from os.path import exists, join
from shutil import copyfile
from urllib.request import urlretrieve

from biom import load_table
from biom.cli.util import write_biom_table

from tax_credit.taxa_manipulator import import_taxonomy_to_dict


def extract_mockrobiota_dataset_metadata(mockrobiota_dir, communities):
    '''Extract mock community metadata from mockrobiota dataset-metadata.tsv
    files
    mockrobiota_dir: PATH to mockrobiota directory
    communities: LIST of mock communities to extract
    '''
    dataset_metadata_dict = dict()
    for community in communities:
        dataset_metadata = import_taxonomy_to_dict(
            join(mockrobiota_dir, "data", community, "dataset-metadata.tsv"))
        dataset_metadata_dict[community] = \
            (dataset_metadata['raw-data-url-forward-read'],
             dataset_metadata['raw-data-url-index-read'],
             dataset_metadata['target-gene'])
    return dataset_metadata_dict


def amend_biom_taxonomy_ids(biom_table,
                            empty_taxonomy=['k__', 'p__', 'c__', 'o__',
                                            'f__', 'g__', 's__'],
                            clean_obs_ids=True):
    '''Convert biom table taxonomy strings so that strings with incomplete
    taxonomies are filled out with ambiguous labels
    '''
    if clean_obs_ids is True:
        clean_taxonomy_ids(biom_table)

    new_ids = {}
    for taxa in biom_table.ids(axis='observation'):
        old_taxonomy = taxa.split(';')
        if len(old_taxonomy) < len(empty_taxonomy):
            new_taxonomy = empty_taxonomy
            for i in range(len(old_taxonomy)):
                new_taxonomy[i] = old_taxonomy[i]
            new_ids[taxa] = ';'.join(new_taxonomy)
        else:
            new_ids[taxa] = ';'.join(old_taxonomy)

    return biom_table.update_ids(new_ids, axis='observation')


def clean_taxonomy_ids(table, delete_chars='[]()'):
    new_ids = {obs_id: obs_id.translate(str.maketrans('', '', delete_chars))
               for obs_id in table.ids(axis="observation")}
    return table.update_ids(new_ids, axis='observation')


def extract_mockrobiota_data(communities, community_md, ref_dbs,
                             mockrobiota_dir, mock_data_dir,
                             expected_data_dir, biom_fn='table.L6-taxa.biom'):
    '''Extract sample metadata, raw data files, and expected taxonomy

    from mockrobiota, copy to new destination
    communities: LIST of mock communities to extract
    community_md: DICT of metadata for mock community.
        see extract_mockrobiota_dataset_metadata()
    ref_dbs = DICT mapping marker_gene to reference set names
    mockrobiota_dir = PATH to mockrobiota repo directory
    mock_data_dir = PATH to destination directory
    expected_data_dir = PATH to destination for expected taxonomy files
    '''
    for community in communities:
        # extract dataset metadata/params
        forward_read_url, index_read_url, marker_gene = community_md[community]
        ref_outdir, ref_indir, ref_version, otu_id = ref_dbs[marker_gene][0:4]

        # mockrobiota source directory
        mockrobiota_community_dir = join(mockrobiota_dir, "data", community)

        # new mock community directory
        community_dir = join(mock_data_dir, community)
        seqs_dir = join(community_dir, 'raw_seqs')
        if not exists(seqs_dir):
            makedirs(seqs_dir)
        # copy sample-metadata.tsv
        copyfile(join(mockrobiota_community_dir, 'sample-metadata.tsv'),
                 join(community_dir, 'sample-metadata.tsv'))
        # download raw data files
        for file_url_dest in [(forward_read_url, 'sequences.fastq.gz'),
                              (index_read_url, 'barcodes.fastq.gz')]:
            destination = join(seqs_dir, file_url_dest[1])
            if not exists(destination) and file_url_dest[0] != 'NA':
                try:
                    urlretrieve(file_url_dest[0], destination)
                except ValueError:
                    print('Error retrieving {0}'.format(file_url_dest[0]))

        # new directory containing expected taxonomy assignments at each level
        expected_taxa_dir = join(expected_data_dir, community,
                                 ref_outdir, "expected")
        if not exists(expected_taxa_dir):
            makedirs(expected_taxa_dir)
        # copy expected taxonomy.tsv and convert to biom
        exp_taxa_fp = join(expected_taxa_dir, 'expected-taxonomy.tsv')
        exp_biom_fp = join(expected_taxa_dir, biom_fn)
        copyfile(join(mockrobiota_community_dir, ref_indir,
                      ref_version, otu_id, 'expected-taxonomy.tsv'),
                 exp_taxa_fp)
        newbiom = amend_biom_taxonomy_ids(load_table(exp_taxa_fp))
        # add taxonomy ids (names) as observation metadata
        metadata = {sid: {'taxonomy': sid.split(';')}
                    for sid in newbiom.ids(axis='observation')}
        newbiom.add_metadata(metadata, 'observation')
        write_biom_table(newbiom, 'hdf5', exp_biom_fp)
