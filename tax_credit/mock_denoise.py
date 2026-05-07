#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2016--, tax-credit development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""QIIME 2 demux, DADA2 denoise, feature table export, and phylogeny."""

from os import makedirs
from os.path import exists, join
from shutil import copyfile, rmtree

from biom import Table
from biom.cli.util import write_biom_table
import qiime2
from qiime2.plugins import feature_table, demux, dada2, alignment, phylogeny


def batch_demux(communities,
                mock_data_dir,
                demux_params,
                raw_seqs='raw_seqs',
                metadata_fn='sample-metadata.tsv',
                index_col='BarcodeSequence',
                demux_fn='demux-seqs.qza',
                summary_fn='demux_summary.qzv',
                qual_plot_fn='demux_plot_qual.qzv',
                n_qual_plots=1):
    '''raw fastq -> SampleData[SequencesWithQuality]

    Demultiplex raw fastq files, summarize demux, and plot qual scores on
    batch of files

        communities: list
            list of mock communities to extract
        mock_data_dir = filepath
            PATH to data dir containing communities
        demux_params = dict
            DICT of TUPLES listing demux parameters.
            {community : (rev_comp_barcodes, rev_comp_mapping_barcodes)}
        raw_seqs = str
            name of directory containing raw seqs, located in directory
            mock_data_dir/community. Must only contain the following files:
                sequences.fastq.gz
                barcodes.fastq.gz
        metadata_fn = filename
            fn of metadata map, located in mock_data_dir/community
        index_col = str
            column name of metadata column containing index/barcode sequences
        demux_fn = str
            filename to save demux artifact
        summary_fn = str
            filename to save demux summary visualization
        qual_plot_fn = str
            filename to save fastq quality plot visualization
        n_qual_plots = int
            number of fastq quality plots to print
    '''

    for community in communities:
        # extract dataset metadata/params
        community_dir = join(mock_data_dir, community)
        seqs_dir = join(community_dir, raw_seqs)
        sample_md = join(community_dir, metadata_fn)

        # demultiplex
        if demux_params[community][0] is True:
            demux_to_plot_quality(seqs_dir=seqs_dir,
                                  sample_md=sample_md,
                                  community_dir=community_dir,
                                  index_col=index_col,
                                  rc_barcodes=demux_params[community][1],
                                  rc_map_barcodes=demux_params[community][2],
                                  demux_fn=demux_fn,
                                  summary_fn=summary_fn,
                                  qual_plot_fn=qual_plot_fn,
                                  n_qual_plots=n_qual_plots)

        else:
            load_demux_seqs(community_dir, seqs_dir, demux_fn, sample_md)

        print("{0} complete".format(community))


def demux_to_plot_quality(seqs_dir,
                          sample_md,
                          community_dir,
                          index_col='BarcodeSequence',
                          rc_barcodes=False,
                          rc_map_barcodes=False,
                          demux_fn='demux-seqs.qza',
                          summary_fn='demux_summary.qzv',
                          qual_plot_fn='demux_plot_qual.qzv',
                          n_qual_plots=1):

    '''raw fastq -> SampleData[SequencesWithQuality]

    Demultiplex raw fastq files, summarize demux, and plot qual scores.

        seqs_dir = filepath
            directory containing fastq sequences
        sample_md = filepath
            filepath to sample metadata file
        community_dir: path
            destination directory to print results
        index_col = str
            column name of metadata column containing index/barcode sequences
        rc_barcodes = bool
            reverse complement barcodes prior to demultiplexing?
        rc_map_barcodes = bool
            reverse complement metadata map barcodes prior to demultiplexing?
        demux_fn = str
            filename to save demux artifact
        summary_fn = str
            filename to save demux summary visualization
        qual_plot_fn = str
            filename to save fastq quality plot visualization
        n_qual_plots = int
            number of fastq quality plots to print
    '''

    # import fastq to qiime2 artifact
    seq_artifact = qiime2.Artifact.import_data("RawSequences", seqs_dir)

    # demultiplex
    barcodes = qiime2.metadata.MetadataCategory.load(sample_md, index_col)
    demux_seqs = demux.methods.emp_single(
        seqs=seq_artifact, barcodes=barcodes, rev_comp_barcodes=rc_barcodes,
        rev_comp_mapping_barcodes=rc_map_barcodes)
    demux_seqs.per_sample_sequences.save(join(community_dir, demux_fn))

    visualize_qual(demux_seqs, community_dir, summary_fn,
                   qual_plot_fn, n_qual_plots)


def visualize_qual(demux_seqs, community_dir, summary_fn,
                   qual_plot_fn, n_qual_plots=1):

    '''visualize demux summary and fastq quality plots.'''

    # demultiplexing summary
    try:
        dsum = demux.visualizers.summarize(demux_seqs.per_sample_sequences)
        dsum.visualization.save(join(community_dir, summary_fn))
    except TypeError:
        # Fails if N=1 https://github.com/qiime2/q2-demux/issues/20
        print("Could not print demux summary: TypeError")

    # view fastq quality plots
    qualplot = dada2.visualizers.plot_qualities(
        n=n_qual_plots, demultiplexed_seqs=demux_seqs.per_sample_sequences)
    qualplot.visualization.save(join(community_dir, qual_plot_fn))


def load_demux_seqs(community_dir, seqs_dir, demux_fn, sample_md):
    '''Load demultiplexed sequences into artifact and write artifact.'''

    # extract sample name
    with open(sample_md, 'r') as md_in:
        md_in.readline()
        sample_id = md_in.readline().strip().split('\t')[0]

    # assemble artifact components
    tmpdir = join(seqs_dir, 'tmp')
    seq_fn = '{0}_1_L001_R1_001.fastq.gz'.format(sample_id)
    if not exists(tmpdir):
        makedirs(tmpdir)
    copyfile(join(seqs_dir, 'sequences.fastq.gz'),
             join(tmpdir, seq_fn))
    with open(join(tmpdir, 'MANIFEST'), 'w') as manifest:
        manifest.write('sample-id,filename,direction\n')
        manifest.write('{0},{1},forward'.format(sample_id, seq_fn))
    with open(join(tmpdir, 'metadata.yml'), 'w') as yml:
        yml.write('\n')

    seqs = qiime2.Artifact.import_data(
        "SampleData[SequencesWithQuality]", tmpdir)
    rmtree(tmpdir)

    seqs.save(join(community_dir, demux_fn))

    # Visualize currently fails on artifacts containing single sample
    # (seqs has no attribute 'per_sample_sequences')
    # visualize_qual(seqs, community_dir, summary_fn,
    #               qual_plot_fn, n_qual_plots)


def denoise_to_phylogeny(communities,
                         mock_data_dir,
                         trim_params,
                         demux_seqs_fn='demux-seqs.qza',
                         rep_seqs_fn='rep_seqs.qza',
                         feature_table_fn='feature_table.qza',
                         summary_fn='feature_table_summary.qzv'):

    '''SampleData[SequencesWithQuality] -> FeatureData[Sequence] +
                                           FeatureTable[Frequency]
        denoise fastqs with dada2, create feature table, rep_seqs,
        and view stats on batch of files.

        communities: LIST
            list of mock communities to extract
        mock_data_dir = filepath
            path to data dir containing communities
        trim_params = dict
            DICT of TUPLES listing dada2 trimming parameters.
            {community : (trim_left, trunc_len)}
        demux_seqs_fn = str
            filename of SampleData[SequencesWithQuality] Artifact to load
        rep_seqs_fn = str
            filename of representative sequences output Artifact
        feature_table_fn = str
            filename of feature table output Artifact
        summary_fn = str
            filename of feature table summary output visualization
    '''

    for community in communities:
        trim_left, trunc_len, buildtree = trim_params[community]
        community_dir = join(mock_data_dir, community)

        # denoise with dada2
        demux_seqs = qiime2.Artifact.load(join(community_dir, demux_seqs_fn))
        biom_table, rep_seqs = denoise_to_feature_table(
            demux_seqs, trim_left, trunc_len, community_dir)

        # Build phylogeny
        if buildtree is True:
            seqs_to_tree(rep_seqs, community_dir)

        print("{0} complete".format(community))


def denoise_to_feature_table(demux_seqs,
                             trim_left,
                             trunc_len,
                             community_dir,
                             rep_seqs_fn='rep_seqs',
                             feature_table_fn='feature_table.qza',
                             biom_table_fn='feature_table.biom',
                             summary_fn='feature_table_summary.qzv'):
    '''SampleData[SequencesWithQuality] -> FeatureData[Sequence] +
                                           FeatureTable[Frequency]
    denoise fastqs with dada2, create feature table, rep_seqs,
        and view stats.

        demux_seqs = SampleData[SequencesWithQuality]
            demultiplexed seqs output from qiime2.demux.methods.emp()
        trim_left = int
            trim X bases from 5' end
        trunc_len = int
            length to truncate all sequences
        community_dir: path
            destination directory to print results
        rep_seqs_fn = str
            filename of representative sequences output Artifact
        feature_table_fn = str
            filename of feature table output Artifact
        summary_fn = str
            filename of feature table summary output visualization
    '''
    biom_table, rep_seqs = dada2.methods.denoise_single(
        demux_seqs, trim_left=trim_left, trunc_len=trunc_len)
    # save Artifact
    rep_seqs.save(join(community_dir, rep_seqs_fn))

    # save biom Artifact
    biom_table.save(join(community_dir, feature_table_fn))
    biom_table_fp = join(community_dir, biom_table_fn)
    write_biom_table(biom_table.view(Table), 'hdf5', biom_table_fp)

    # summarize feature table
    feature_table_summary = feature_table.visualizers.summarize(biom_table)
    feature_table_summary.visualization.save(join(community_dir, summary_fn))

    return biom_table, rep_seqs


def seqs_to_tree(rep_seqs, community_dir, filename='phylogeny.qza'):
    '''FeatureData[Sequence] -> phylogeny

    rep_seqs: FeatureData[Sequence] Artifact
        representative sequences from dada2
    community_dir: path
        destination directory to print results
    '''
    aligned_seqs = alignment.methods.mafft(rep_seqs)
    m_aln = alignment.methods.mask(aligned_seqs.alignment)
    unrooted_tree = phylogeny.methods.fasttree(m_aln.masked_alignment)
    tree = phylogeny.methods.midpoint_root(unrooted_tree.tree)
    tree.rooted_tree.save(join(community_dir, filename))
