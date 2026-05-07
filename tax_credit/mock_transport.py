#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2016--, tax-credit development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""Copy processed mock community artifacts into the tax-credit repo layout."""

from os import makedirs
from os.path import basename, exists, join
from shutil import copyfile

import qiime2
from q2_types.feature_data import DNAIterator
from skbio import TreeNode, io


def transport_to_repo(communities,
                      mock_data_dir,
                      project_dir,
                      sample_type_dirname='mock-community',
                      rep_seqs_fn='rep_seqs.qza',
                      feature_table_fn='feature_table.qza',
                      tree_fn='phylogeny.qza',
                      sample_md_fn='sample-metadata.tsv',
                      biom_table_fn='feature_table.biom',
                      fasta_fn='rep_seqs.fna',
                      newick_fn='phylogeny.tre'):
    '''Copy essential mock community data to tax-credit repo

    communities: list
        list of dir names in mock_data_dir, a.k.a. names of mock communities
    mock_data_dir: path
        source directory containing mock communities dirs of results
    project_dir: path
        path to tax-credit repo directory
    sample_type_dirname: str
        name of destination directory to contain communities dirs. The analog
        of mock_data_dir in the repo, dirs for individual communities will be
        located in project_dir/data/sample_type_dirname/community
    rep_seqs_fn: str
        name of rep seqs FeatureData[Sequence] Artifact in community_dir
    feature_table_fn: str
        name of rep seqs FeatureTable[Frequency] Artifact in community_dir
    tree_fn: str
        name of Phylogeny[Rooted] Artifact in community_dir
    sample_md_fn: str
        name of metadata mapping file in community_dir
    biom_table_fn: str
        destination name of biom table in project_dir
    fasta_fn: str
        destination name of fasta file in project_dir
    newick_fn: str
        destination name of newick format tree in project_dir
    '''

    for community in communities:
        community_dir = join(mock_data_dir, community)

        # Define base dir destination for mock community directories
        repo_destination = join(project_dir, "data",
                                sample_type_dirname, community)
        if not exists(repo_destination):
            makedirs(repo_destination)

        # Files to move
        rep_seqs = join(community_dir, rep_seqs_fn)
        feature_table = join(community_dir, feature_table_fn)
        tree = join(community_dir, tree_fn)
        sample_md = join(community_dir, sample_md_fn)

        biom_table_fp = join(community_dir, biom_table_fn)
        rep_seqs_fp = join(community_dir, fasta_fn)
        tree_fp = join(community_dir, newick_fn)

        # Extract biom, tree, rep_seqs
        rep_seqs_fna = qiime2.Artifact.load(rep_seqs).view(DNAIterator)
        io.write(rep_seqs_fna.generator, format='fasta', into=rep_seqs_fp)

        if exists(tree):
            qiime2.Artifact.load(tree).view(TreeNode).write(tree_fp)

        # Move to repo:
        for f in [rep_seqs, feature_table, tree, sample_md,
                  biom_table_fp, rep_seqs_fp, tree_fp]:
            if exists(f):
                copyfile(f, join(repo_destination, basename(f)))
