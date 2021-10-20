# vim: fdm=indent
'''
author:     Fabio Zanini
date:       04/02/21
content:    Try and demux samples based on gender/SNPs
'''
import os
import sys
import argparse
import numpy as np
import pandas as pd

sys.path.append('/home/fabio/university/postdoc/singlet')
import singlet

import matplotlib as mpl
mpl.rcParams['font.size'] = 7
import matplotlib.pyplot as plt
import seaborn as sns


def get_sex(samplenames):
    samples = ['iMSC', 'AdMSC']
    snv_df = []
    for sample in samples:
        fn_snv = f'../../data/sequencing/experiments/201016/snv_demux/{sample}.tsv'
        tmp = pd.read_csv(fn_snv, sep='\t', index_col=0, encoding='utf-8')
        tmp.index = tmp.index + '-' + sample
        snv_df.append(tmp)
    snv_df = pd.concat(snv_df)
    # Add missing samples (nobody is perfect)
    missing = set(samplenames) - set(snv_df.index)
    for cellID in missing:
        snv_df.loc[cellID] = 'miss'
    snv_df = snv_df.loc[samplenames]

    # We are lucky and can just set a replace
    gender = snv_df.replace({
        'SNG-0': 'F', 'SNG-1': 'M', 'SNG-2': 'M',
        'DBL-1': 'doublet', 'DBL-2': 'doublet', 'miss': 'unknown',

    })
    return gender


def add_sex(ds, inplace=False):
    sex = get_sex(ds.samplenames)
    if not inplace:
        return gender

    ds.obs['sex'] = sex
    ds.obs['sample_unique'] = ds.obs['sample'].astype(str) + '_' + ds.obs['sex']


if __name__ == '__main__':

    pa = argparse.ArgumentParser()
    pa.add_argument('--samples', default='iMSC,AdMSC')
    args = pa.parse_args()

    samples = args.samples.split(',')

    print('Load data from file')
    fn = '../../data/sequencing/experiments/201016_A00152_0317_AHMTN3DSXY/normalized.h5ad'
    ds0 = singlet.Dataset(
        dataset={'path': fn},
        )
    # Undo normalization
    ds0.counts *= 1e-6 * ds0.obs['coverage']
    ds0.obs['n_genes'] = (ds0.counts > 0).sum(axis=0)

    ds0.query_samples_by_metadata('sample in @samples', local_dict=locals(), inplace=True)

    print('Add chromosomes')
    genes_chroms = pd.read_csv(
        '../../data/genome/genes_chromosomes.tsv',
        sep='\t',
        encoding='utf-8',
        )
    genes_chroms.drop_duplicates('Gene name', inplace=True)
    genes_chroms.set_index('Gene name', inplace=True)
    genes_chroms = genes_chroms['Chromosome/scaffold name']

    missing = set(ds0.featurenames) - set(genes_chroms.index)
    for gene in missing:
        genes_chroms.loc[gene] = ''
    genes_chroms = genes_chroms.loc[ds0.featurenames]

    ds0.featuresheet['Chromosome'] = genes_chroms

    print('Check Y chromosome genes')
    genes_y = genes_chroms[genes_chroms == 'Y'].index
    dsy = ds0.query_features_by_name(genes_y)
    gender = pd.cut(
            dsy.counts.sum(axis=0),
            bins=[0, 1, 100, 1e6],
            include_lowest=True,
            labels=['F', 'unknown', 'M'],
            )
    ds0.obs['Y_chrom'] = dsy.counts.sum(axis=0)

    print('Check X chromosome genes')
    genes_x = genes_chroms[genes_chroms == 'X'].index
    dsx = ds0.query_features_by_name(genes_x)
    #gender = pd.cut(
    #        dsy.counts.sum(axis=0),
    #        bins=[0, 1, 100, 1e6],
    #        include_lowest=True,
    #        labels=['F', 'unknown', 'M'],
    #        )
    #ds0.obs['gender'] = gender
    ds0.obs['X_chrom'] = dsx.counts.sum(axis=0)

    print('Cross-check against SNV demux')
    snv_df = []
    for sample in samples:
        fn_snv = f'../../data/sequencing/experiments/201016/snv_demux/{sample}.tsv'
        tmp = pd.read_csv(fn_snv, sep='\t', index_col=0, encoding='utf-8')
        tmp.index = tmp.index + '-' + sample
        snv_df.append(tmp)
    snv_df = pd.concat(snv_df)
    # Add missing samples (nobody is perfect)
    missing = set(ds0.samplenames) - set(snv_df.index)
    for cellID in missing:
        snv_df.loc[cellID] = 'miss'
    snv_df = snv_df.loc[ds0.samplenames]
    ds0.obs['snv_cluster'] = snv_df

    # We are lucky and can just set a replace
    gender = snv_df.replace({
        'SNG-0': 'F', 'SNG-1': 'M', 'SNG-2': 'M',
        'DBL-1': 'doublet', 'DBL-2': 'doublet', 'miss': 'unknown',

    })
    ds0.obs['gender'] = gender

    print('Plot X, Y, cluster')
    dsp = ds0.split(['sample', 'snv_cluster'])
    clusters = ds0.obs['snv_cluster'].unique()
    orders = {
        'iMSC': ['SNG-1', 'SNG-0', 'DBL-2', 'miss'],
        'AdMSC': ['SNG-2', 'SNG-0', 'DBL-1', 'miss'],
    }
    fig, axs = plt.subplots(2, 2, figsize=(4, 4), sharex=True, sharey=True)
    cmap = {
        ('iMSC', 'SNG-1'): 'steelblue',
        ('iMSC', 'SNG-0'): 'deeppink',
        ('iMSC', 'miss'): 'grey',
        ('iMSC', 'DBL-2'): 'red',
        ('AdMSC', 'SNG-2'): 'steelblue',
        ('AdMSC', 'SNG-0'): 'deeppink',
        ('AdMSC', 'DBL-1'): 'red',
        ('AdMSC', 'miss'): 'grey',
    }
    labeld = {
        'SNG-1': 'boy', 'SNG-0': 'girl', 'DBL-2': 'doublet',
        'SNG-2': 'boy', 'DBL-1': 'doublet',
    }
    for sample, axs_row in zip(samples, axs):
        for ax, chrom in zip(axs_row, ['Y', 'X']):
            ax.set_title(chrom)
            for clu in orders[sample]:
                if (sample, clu) not in dsp:
                    continue
                x = np.sort(dsp[(sample, clu)].obs[f'{chrom}_chrom'])
                y = 1.0 - np.linspace(0, 1, len(x))
                ls = '-' if 'SNG' in clu else '--'
                ax.plot(
                        x, y,
                        lw=1.5 + 0.5 * clu.startswith('SNG'),
                        alpha=0.7, color=cmap[(sample, clu)],
                        label=labeld.get(clu, clu),
                        )
            if ax == axs[0, 0]:
                ax.legend()
            if ax == axs_row[0]:
                ax.set_ylabel('Fraction of cells > x')
            ax.grid(True)
            ax.set_xscale('log')
            ax.set_xlabel(f'# counts on chromosome {chrom}')
    fig.text(0.02, 0.75, samples[0])
    fig.text(0.02, 0.27, samples[1])
    fig.tight_layout(rect=(0.11, 0, 1, 1))


    if False:
        print('Plot Y vs X')
        x = dsx.counts.sum(axis=0)
        y = dsy.counts.sum(axis=0)
        idx = ds0.obs['n_genes'] > 2000
        fig, ax = plt.subplots()
        ax.scatter(x + 0.1, y + 0.1, alpha=0.3)
        ax.set_xlabel('X counts')
        ax.set_ylabel('Y counts')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.grid(True)
        fig.tight_layout()


