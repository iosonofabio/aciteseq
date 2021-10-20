# vim: fdm=indent
'''
author:     Fabio Zanini
date:       16/01/21
content:    Load and QC data
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

from demux_samples import add_sex
from read_antibody_counts import add_antibodies


def get_markers(dsa):
    cats = dsa.samplenames
    markers = {}
    for y in cats:
        oth = [x for x in cats if x != y]
        fc = np.log10(0.1 + dsa.counts[y]) - \
            np.log10(0.1 + dsa.counts[oth]).mean(axis=1)
        markers[y] = fc
    markers = pd.DataFrame(markers)
    return markers


def get_fc_pairs(dsa):
    cats = dsa.samplenames
    fcs = {}
    for i, y in enumerate(cats):
        for y2 in cats[:i]:
            fc = np.log10(0.1 + dsa.counts[y]) - \
                np.log10(0.1 + dsa.counts[y2])
            fcs[(y, y2)] = fc
            fcs[(y2, y)] = -fc
    fcs = pd.DataFrame(fcs)
    return fcs


if __name__ == '__main__':

    pa = argparse.ArgumentParser()
    pa.add_argument('--save', action='store_true')
    args = pa.parse_args()
    save = args.save

    print('Load data from file')
    fn = '../../data/sequencing/experiments/201016_A00152_0317_AHMTN3DSXY/normalized.h5ad'
    ds0 = singlet.Dataset(
        dataset={'path': fn},
        )
    ds0.obs['n_genes'] = (ds0.counts > 0).sum(axis=0)
    dsp = ds0.split('sample')

    print('Number of reads/genes')
    samples = ['iMSC', 'AdMSC']
    cmap = {'iMSC': 'purple', 'AdMSC': 'darkorange'}
    titlemap = {'coverage': 'Number of UMIs', 'n_genes': 'Number of genes'}
    xlmap = {'coverage': 3e2, 'n_genes': 3e1}
    fig, axs = plt.subplots(1, 2, figsize=(3, 1.5), sharey=True)
    for ax, col in zip(axs, ['coverage', 'n_genes']):
        for sample in samples:
            dsi = dsp[sample]
            x = np.sort(dsi.obs[col].values) + 0.1
            y = 1.0 - np.linspace(0, 1, len(x))
            ax.plot(x, y, lw=2, label=sample, color=cmap[sample])
        ax.grid(True)
        ax.set_xlim(left=xlmap[col])
        ax.set_xscale('log')
        ax.set_ylim(0.001, 1.01)
        ax.set_title(titlemap[col])
        if ax == axs[0]:
            ax.set_ylabel('Fraction of cells > x')
            ax.legend(loc='lower left')
    fig.tight_layout()
    if save:
        fig.savefig('../../figures/qc/n_umis_genes.png')

    ds = ds0.query_samples_by_metadata('coverage >= 1e4')

    fn_umap = '../../figures/qc/umap_coords.tsv'
    if not os.path.isfile(fn_umap):
        print('Feature selection')
        if False:
            features = ds.feature_selection.overdispersed_within_groups(
                'sample', inplace=False)

        features = ds.featurenames
        idx = ~features.str.startswith('MT-')
        idx &= ~features.str.startswith('RPS')
        idx &= ~features.str.startswith('RPL')

        if False:
            dsm = ds.counts.mean(axis=1)
            features = dsm.loc[idx].nlargest(500).index

        features = ds.featurenames
        features = features[idx]
        dsf = ds.query_features_by_name(features)

        #dsf.counts.log(inplace=True)
        features = dsf.feature_selection.overdispersed(n_features=1000)
        dsf = ds.query_features_by_name(features)
        #dsf.counts.unlog(inplace=True)

        print('PCA')
        dsc = dsf.dimensionality.pca(n_dims=25, robust=False, return_dataset='samples')

        print('Dimensionality reduction')
        vs = dsc.dimensionality.umap()
        ds.obs['vs1'] = vs.values.T[0]
        ds.obs['vs2'] = vs.values.T[1]
        if save:
            vs.to_csv(fn_umap, sep='\t', index=True)

        print('Get knn graph')
        #knn = ds.graph.knn(
        #        'samples', n_neighbors=10, return_kind='edges',
        #        threshold=0.5,
        #        )
        g = ds.graph.knn(
                'samples', n_neighbors=10, return_kind='igraph',
                threshold=0.5,
                )

        print('Clustering')
        ds.obs['cluster'] = membership = ds.cluster.leiden(
                'samples',
                g,
                resolution_parameter=0.001,
                ).astype(str)

        if save:
            vsc = vs.copy()
            vsc['cluster'] = membership
            vsc.to_csv(fn_umap, sep='\t', index=True)

    else:
        vsc = pd.read_csv(fn_umap, sep='\t', index_col=0)
        vs = vsc.iloc[:, :2]
        ds.obs['cluster'] = membership = vsc['cluster'].astype(str)

    clusters = membership.unique()
    ncl = len(clusters)

    print('Markers')
    dsa = ds.average('samples', 'cluster')
    idxe = (dsa.counts > 100).any(axis=1)
    fcd = get_markers(dsa)
    fcpd = get_fc_pairs(dsa)

    print('Get antibodies and sex')
    add_sex(ds, inplace=True)
    add_antibodies(ds, inplace=True)

    print('Plot dimensionality reduction')
    genes = ['sample_unique',
             'cluster',
             'coverage',
             'n_genes',
             #'CCR4',
             'NAMPT',
             #'CCR5',  # MISSING??
             'TAGLN',
             'ACTA2',
             'LOX',
             'TGFBI',
             'COL1A1',
             #'VIM', 'TIMP3', 'TPM2', 'DCN', 'XIST',
             #'CD36',
             'CCND1',
             #'KIF1C',
             #'IGLC3', 'JCHAIN',
             #'COL15A1',
             'TNFAIP6',
             'POSTN',
             'NEK10',
             'CDKN1A',
             'MKI67', #'CCNB1', 'CDC20',
             'HHIP', 'MCAM',
             #'DES',
             'BEX1',
             'MEST',
             'SFRP4',
             'KRT19',
             'COMP',
             'NEFM',
             'RPS4Y1',
             'ACTC1',
             'SOD3',
             #'NTF4',
             #'RBP1',
             'NECTIN3',
             #'SSR3',
             'TRIP10',
             #'COL6A6',
             #'IL13',
             #'BRI3',
             #'CAV1',
             #'DRAP1',
             #'EIF5A',
             'ID1',
             'SLC9A3R2',
             'SCUBE3',
             #'TUBA1B',
             'MGARP',
             'NTN4',
             'CXCL1',
             ]
    cmaps = {
        'sample': {'iMSC': 'purple', 'AdMSC': 'darkorange'},
        'sample_unique': {
            'iMSC_M': 'dodgerblue',
            'AdMSC_M': 'steelblue',
            'iMSC_F': 'purple',
            'AdMSC_F': 'deeppink',
            'iMSC_doublet': 'gold',
            'AdMSC_doublet': 'brown',
            },
        'cluster': dict(zip(clusters, sns.color_palette('husl', ncl))),
    }
    if True:
        fig, axs = plt.subplots(2, 2, figsize=(3, 3), sharex=True, sharey=True)
        axs = axs.ravel()
        for gene, ax in zip(genes, axs):
            ds.plot.scatter_reduced(
                    vs,
                    ax=ax,
                    cmap=cmaps.get(gene, 'viridis'),
                    color_log=(gene in ['coverage']) or (gene in ds.featurenames),
                    color_by=gene,
                    alpha=0.05,
                    s=10,
                    )
            ax.set_title(gene)
            if gene == 'cluster':
                for icst, cst in enumerate(clusters):
                    xsm, ysm = vs.loc[ds.samplesheet[gene] == cst].mean(axis=0)
                    ax.text(xsm, ysm, cst, ha='center', va='center')

        fig.tight_layout()

    if False:
        fig, axs = plt.subplots(4, 7, figsize=(8.5, 5), sharex=True, sharey=True)
        axs = axs.ravel()
        for gene, ax in zip(genes, axs):
            ds.plot.scatter_reduced(
                    vs,
                    ax=ax,
                    cmap=cmaps.get(gene, 'viridis'),
                    color_log=(gene in ['coverage']) or (gene in ds.featurenames),
                    color_by=gene,
                    alpha=0.05,
                    s=10,
                    )
            ax.set_title(gene)
            if gene == 'cluster':
                for icst, cst in enumerate(clusters):
                    xsm, ysm = vs.loc[ds.samplesheet[gene] == cst].mean(axis=0)
                    ax.text(xsm, ysm, cst, ha='center', va='center')
        fig.tight_layout()
        if save:
            fig.savefig('../../figures/qc/many_umaps.png')


    print('Dot plot')
    from scipy.spatial.distance import pdist
    from scipy.cluster.hierarchy import linkage, leaves_list
    mat = dsa.counts.loc[genes[4:]]
    pdis = pdist(mat.values)
    z = linkage(pdis, method='average', optimal_ordering=True)
    genes_order = np.array(genes[4:])[leaves_list(z)]
    pdis = pdist(mat.values.T)
    z = linkage(pdis, method='average', optimal_ordering=True)
    cluster_order = np.array(clusters)[leaves_list(z)]

    fig, axs = plt.subplots(
            1, 2, figsize=(6, 3), gridspec_kw={'width_ratios': [12, 1]},
            sharey=True,
            )
    ds.plot.dot_plot(
            group_by='cluster', group_order=cluster_order,
            plot_list=genes_order, ax=axs[0])
    axs[0].set_ylabel('Cluster #')

    dfab = ds.obs.groupby(['cluster', 'sample']).size().unstack().loc[cluster_order]
    dfab.plot.barh(color=cmaps['sample'], legend=None, ax=axs[1])
    axs[1].set_xlim(0.9, dfab.values.max() * 1.5)
    #axs[1].set_xscale('log')
    #axs[1].set_xticks([1e2], minor=True)
    axs[1].set_xlabel('# cells\nin sample')
    fig.tight_layout()
    if save:
        fig.savefig('../../figures/qc/dotplot_clusters.png')

    print('Print distributions of known pluripotency genes')
    from scipy.stats import gaussian_kde
    dspc = ds.split('cluster')
    genes_pp = [
        'POU5F1',  # OCT4
        'NANOG',
        'KLF4',
        'MYC',
    ]
    common_names = {'POU5F1': 'OCT4'}
    kind = 'cumulative'
    fig, axs = plt.subplots(2, 2, figsize=(3, 3), sharey=True, sharex=True)
    axs = axs.ravel()
    for gene, ax in zip(genes_pp, axs):
        for cl in range(ncl):
            x = np.log10(0.1 + dspc[str(cl)].counts.loc[gene].values)

            if kind == 'kde':
                xfit = np.linspace(-1, 4, 100)
                if x.max() > -1:
                    yfit = gaussian_kde(x, bw_method=0.3)(xfit)
                    yfit /= yfit.max()
                else:
                    yfit = np.zeros_like(xfit)
                    yfit[0] = 1
                yfit /= 1.1
                ax.fill_between(
                    xfit, cl, cl + yfit, color=cmaps['cluster'][str(cl)],
                    )
            else:
                x = np.sort(x)
                y = 1.0 - np.linspace(0, 1, len(x))
                ax.plot(x, y, color=cmaps['cluster'][str(cl)], lw=1,
                        label=str(cl))

        ax.grid(True)
        ax.set_yticks(np.arange(ncl))
        ax.set_yticklabels([str(x) for x in range(ncl)])
        ax.set_xticks([-1, 1, 3])
        ax.set_xticklabels(['$0$', '$10$', '$10^3$'])
        if gene in common_names:
            title = f'{gene}/'+common_names[gene]
        else:
            title = gene
        ax.set_title(title)

        if (ax == axs[0]) and (kind != 'kde'):
            ax.legend(ncol=2, fontsize=7)

    fig.tight_layout()

    #print('Avani wants to know some specific genes')
    #genesa = 

    plt.ion()
    plt.show()
