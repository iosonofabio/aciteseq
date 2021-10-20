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

import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append('/home/fabio/university/PI/projects/avani_cite-seq')
from aciteseq.ingest.load_iMS_data import load_data
from aciteseq.ingest.antibodies import get_antibody_controls


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
    pa.add_argument('--include-immune', action='store_true')
    pa.add_argument('--include-isotype-controls', action='store_true')
    pa.add_argument('--save', action='store_true')
    args = pa.parse_args()
    save = args.save

    ds = load_data()

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
        vsc = vsc.loc[ds.samplenames]
        vs = vsc.iloc[:, :2]
        ds.obs['umap1'] = vs.iloc[:, 0]
        ds.obs['umap2'] = vs.iloc[:, 1]
        ds.obs['cluster'] = membership = vsc['cluster'].astype(str)

    clusters = membership.unique()
    ncl = len(clusters)

    if False:
        print('Look at the distributions of a few immune genes (?)')
        genes = ['CD4', 'CD8A', 'CD3E', 'TRAC',
                 'IGHM', 'IGKC', 'TNF', 'OASL',
                 'CCL5', 'GAL']
        fig, axs = plt.subplots(2, 5, figsize=(10, 4), sharex=True, sharey=True)
        axs = axs.ravel()
        for gene, ax in zip(genes, axs):
            x = np.sort(ds.counts.loc[gene].values)
            y = (1.0 - np.linspace(0, 1, len(x))) * len(x)
            ax.plot(x + 0.1, y, lw=2, color='grey')
            ax.set_title(gene)
        ax.set_xscale('log')
        fig.tight_layout()

    if not args.include_immune:
        print('Exclude immune cells')
        # ~270 cells
        genes_immune = ['CCL5', 'TNF', 'TRAC', 'CD4', 'CD8A', 'IGHM', 'IGKC', 'OASL']
        immune_cells = ds.samplenames[ds.counts.loc[genes_immune].sum(axis=0) > 0]
        ds.exclude_samples_by_name(immune_cells, inplace=True)

        # ~325 cells
        ab_ctrls = get_antibody_controls()
        ab_immune = ab_ctrls['negative']
        immune_cells = ds.samplenames[(ds.obs[ab_immune] >= 10).any(axis=1)]
        ds.exclude_samples_by_name(immune_cells, inplace=True)

    if not args.include_isotype_controls:
        print('Exclude cells with >= 30 counts in isotype control Abs')
        # ~30 cells
        ab_iso = get_antibody_controls()['isotype']
        iso_cells = ds.samplenames[(ds.obs[ab_iso] >= 30).any(axis=1)]
        ds.exclude_samples_by_name(iso_cells, inplace=True)


    print('Add pseudocounts to antibodies')
    ds.obs[ds.antibodies] += 10

    # TODO: WHAT ABOUT GAL??

    print('Markers')
    dsa = ds.average('samples', 'cluster')
    idxe = (dsa.counts > 100).any(axis=1)
    fcd = get_markers(dsa)
    fcpd = get_fc_pairs(dsa)

    print('Plot dimensionality reduction')
    genes = ['sample_unique',
             'sample',
             'sex',
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
        'sex': {'M': 'dodgerblue', 'F': 'deeppink'},
        'sample_unique': {
            'iMSC_M': 'navy',
            'AdMSC_M': 'steelblue',
            'iMSC_F': 'purple',
            'AdMSC_F': 'pink',
            'iMSC_doublet': 'gold',
            'AdMSC_doublet': 'brown',
            },
        'cluster': dict(zip(clusters, sns.color_palette('husl', ncl))),
    }
    if False:
        fig, axs = plt.subplots(2, 3, figsize=(6, 4), sharex=True, sharey=True)
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
            if gene == 'sample':
                labels = ['iMSC', 'AdMSC']
                hs = [ax.scatter([], [], color=cmaps[gene][x]) for x in labels]
                ax.legend(hs, labels, loc='upper right')

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

    if False:
        print('Differential expression between iMSC and AdMSC by sex')
        genes_order = []
        for sex in ['M', 'F']:
            ds_sex = ds.query_samples_by_metadata('sex == @sex', local_dict=locals())
            dsp = ds_sex.split('sample')
            for key, dsi in dsp.items():
                dsi.subsample(n=200, inplace=True)
            comp = dsp['iMSC'].compare(dsp['AdMSC'], method='kolmogorov-smirnov-rich')
            comp.rename(
                    columns={'avg_self': 'avg_iMSC', 'avg_other': 'avg_AdMSC'},
                    inplace=True,
                    )
            genes_order += list(comp.loc[comp['avg_iMSC'] >= 100]
                           .nlargest(10, 'log2_fold_change')
                           .index)
            genes_order += list(comp.loc[comp['avg_AdMSC'] >= 100]
                           .nsmallest(10, 'log2_fold_change')
                           .index)

        cluster_order = ['iMSC_M', 'AdMSC_M', 'iMSC_F', 'AdMSC_F']
        fig, ax = plt.subplots(
                1, 1, figsize=(6, 3),
                )
        ds.plot.dot_plot(
                group_by='sample_unique', group_order=cluster_order,
                plot_list=genes_order, ax=ax)
        ax.set_ylabel('Sample')
        fig.tight_layout()


    if False:
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


    print('Focus on the female only for now')
    dsF = ds.query_samples_by_metadata('sex == "F"')
    dsF.antibodies = ds.antibodies
    dsF.obs['umap1'] = vs.loc[dsF.samplenames].iloc[:, 0]
    dsF.obs['umap2'] = vs.loc[dsF.samplenames].iloc[:, 1]

    if False:
        print('Plot a few embeddings')
        genes = ['sample', 'cluster', 'MKI67', 'ASPM', 'PLAC9', 'MT-ND4L']
        fig, axs = plt.subplots(2, 3, figsize=(6, 4), sharex=True, sharey=True)
        axs = axs.ravel()
        for gene, ax in zip(genes, axs):
            dsF.plot.scatter_reduced(
                    ('umap1', 'umap2'),
                    ax=ax,
                    cmap=cmaps.get(gene, 'viridis'),
                    color_log=(gene in ['coverage']) or (gene in dsF.featurenames),
                    color_by=gene,
                    alpha=0.05,
                    s=10,
                    )
            ax.set_title(gene)
            if gene == 'cluster':
                for icst, cst in enumerate(clusters):
                    xsm, ysm = dsF.obs[['umap1', 'umap2']].loc[dsF.samplesheet[gene] == cst].mean(axis=0)
                    ax.text(xsm, ysm, cst, ha='center', va='center')
            if gene == 'sample':
                labels = ['iMSC', 'AdMSC']
                hs = [ax.scatter([], [], color=cmaps[gene][x]) for x in labels]
                ax.legend(hs, labels, loc='upper right')

        fig.tight_layout()

    if False:
        print('Check antibody controls')
        # TODO: seurat uses Centered Log Transformed normalization, i.e.
        # x -> log(x)  - <log(x)>_{x over all cells}
        # (if I'm not mistaken:
        # https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html
        # https://github.com/satijalab/seurat/blob/9843b843ed0c3429d86203011fda00badeb29c2e/R/preprocessing.R#L2227
        # but that seems like a bad idea because it totally ignores the controls.
        # Asking Gammy
        ab_ctrls = get_antibody_controls()
        ab_allcontrols = sum(list(ab_ctrls.values()), [])
        ab_ctrls['non-control'] = [x for x in ds.antibodies if x not in ab_allcontrols]
        cats = ['isotype', 'negative', 'positive', 'non-control']
        fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(6, 6))
        axs = axs.ravel()
        for ax, abcat in zip(axs, cats):
            absi = ab_ctrls[abcat]
            # Sort by inverse avg abundance
            absi = dsF.obs[absi].mean(axis=0).sort_values(ascending=False).index
            print(absi)
            colors = sns.color_palette('husl', n_colors=len(absi))
            for ab, color in zip(absi, colors):
                #x = np.sort(dfa[ab])
                x = np.sort(dsF.obs[ab] + 0.1)
                y = (1.0 - np.linspace(0, 1, len(x))) * len(x)
                label = ab.split(';')[0].split(' ')[0]
                ax.plot(x, y, lw=2, label=label, color=color, alpha=0.7)
            ax.grid(True)
            if abcat not in ['isotype', 'non-control']:
                ax.legend(fontsize=6)
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_title(abcat)
        axs[0].set_ylabel('Number of cells\nwith more than x Ab UMIs')
        axs[2].set_ylabel('Number of cells\nwith more than x Ab UMIs')
        axs[2].set_xlabel('Ab counts [UMIs]')
        axs[3].set_xlabel('Ab counts [UMIs]')
        fig.tight_layout()

    if False:
        print('Correlate antibodies (excluding controls) with gene expression')
        ab_ctrls = get_antibody_controls()
        ab_allcontrols = sum(list(ab_ctrls.values()), [])
        ab_ctrls['non-control'] = [x for x in ds.antibodies if x not in ab_allcontrols] 

        corr = dsF.correlation.correlate_features_phenotypes(
                features='all', phenotypes=ab_ctrls['non-control'],
                ).fillna(0)

    if False:
        print('DEG within iMSC')
        dsFi = dsF.query_samples_by_metadata('sample == "iMSC"')
        dsFia = dsFi.average('samples', by='cluster', phenotypes=ds.antibodies)
        # Clusters 0, 2, 4 make sense
        comps = {}
        cls = ['0', '2', '4']
        for cl in cls:
            print(cl)
            other = [x for x in cls if x != cl]
            ds1 = dsFi.query_samples_by_metadata('cluster == @cl', local_dict=locals())
            ds2 = dsFi.query_samples_by_metadata('cluster in @other', local_dict=locals())
            ds1.subsample(100, inplace=True)
            ds2.subsample(100, inplace=True)
            comp = ds1.compare(ds2, method='kolmogorov-smirnov-rich')
            comp.rename(columns={'avg_self': f'avg_{cl}'}, inplace=True)
            comps[cl] = comp

        genes = []
        for cl in cls:
            tmp = list(comps[cl].nlargest(10, 'log2_fold_change').index)
            tmp += list(comps[cl].nsmallest(10, 'log2_fold_change').index)
            for gene in tmp:
                if gene not in genes:
                    genes.append(gene)

        fig, ax = plt.subplots(figsize=(3, 8))
        dsFi.plot.dot_plot(
                group_by='cluster', plot_list=genes, layout='vertical',
                group_order=cls,
                ax=ax,
                threshold=0.5,
                )
        fig.tight_layout()

    if True:
        genes = [
            'sample',
            #None, #'cluster',
            'MKI67',
            #'ACTC1', 'SCRG1',
            'SYNPO2',
            #'LMNA',
            #'PAX7',
            'MCAM']
        genes += [
            'CD146 (MCAM); P1H12-CCTTGGATAACATCA',
            #'CD195 (CCR5); J418F1-CCAAAGTAAGAGCCA',
            'CD54; HA58-CTGATAGACTTGAGT',
            #'ICAM1',
        ]
        fig, axs = plt.subplots(1, 6, figsize=(10, 1.8), sharex=True, sharey=True)
        axs = axs.T.ravel()
        for gene, ax in zip(genes, axs):
            dsF.plot.scatter_reduced(
                    ('umap1', 'umap2'),
                    ax=ax,
                    cmap=cmaps.get(gene, 'viridis'),
                    color_log=(gene in ['coverage']) or (gene in ds.featurenames) or (gene in ds.antibodies),
                    color_by=gene,
                    alpha=0.08,
                    s=15,
                    )
            ax.set_axis_off()
            if gene in ds.antibodies:
                ax.set_title('Ab: '+gene.split(';')[0])
            else:
                ax.set_title(gene)
            ax.set_ylim(-4, 8)
            if gene == 'cluster':
                for icst, cst in enumerate(clusters):
                    if cst not in cls + ['5', '6']:
                        continue
                    xsm, ysm = dsF.obs.loc[dsF.samplesheet[gene] == cst, ['umap1', 'umap2']].mean(axis=0)
                    ax.text(xsm, ysm, cst, ha='center', va='center')
            if gene == 'sample':
                for cst in ['iMSC', 'AdMSC']:
                    xsm, ysm = dsF.obs.loc[dsF.samplesheet[gene] == cst, ['umap1', 'umap2']].mean(axis=0)
                    ax.text(xsm, ysm, cst, ha='center', va='center',
                            bbox=dict(
                                facecolor=(1., 1., 1., 0.7),
                                edgecolor='none',
                                pad=4.0))
            fig.add_artist(plt.Rectangle((0.5, 0.03), 0.33, 0.95, facecolor='none', edgecolor='tomato', ls='--'))
            axs[1].arrow(
                    0.7, 0.9, -0.42, -0.7,
                    transform=ax.transAxes, length_includes_head=True,
                    facecolor='k', alpha=0.7, head_width=0.1, overhang=0.4)
            axs[2].arrow(
                    0.7 - 0.42, 0.2, 0.42, 0.7,
                    transform=ax.transAxes, length_includes_head=True,
                    facecolor='k', alpha=0.7, head_width=0.1, overhang=0.4)
            fig.text(0.20, 0.06, 'proliferation', fontsize=9, bbox=dict(
                                facecolor=(1., 1., 1., 0.2),
                                edgecolor='none',
                                pad=3.0))
            fig.text(0.37, 0.06, 'differentiation', fontsize=9, bbox=dict(
                                facecolor=(1., 1., 1., 0.2),
                                edgecolor='none',
                                pad=3.0))
        fig.tight_layout()
        fig.savefig('/home/fabio/university/PI/grants/AU_MRFF_2021/figures/Figure5_raw.png', dpi=300)



    if False:
        print('Differential expression and differential antibodies')
        dsFa = dsF.average('samples', 'cluster', phenotypes=dsF.antibodies)
        dsFa.antibodies = dsF.antibodies

        cls = dsFa.samplenames
        markers = {}
        for cl in cls:
            oth = [x for x in cls if x != cl]
            diff = (np.log2(0.1 + dsFa.counts['0']) -
                    np.log2(0.1 + dsFa.counts[oth].mean(axis=1)))
            markers[cl] = diff.nlargest(10)

        genes = [
            'AC078850.1', 'ACTG2', 'ANKRD1', 'BEX1', 'CD36', 'CEMIP',
            'FAM19A5', 'HIST1H2AC', 'IGFBP5', 'IL6', 'KCTD20', 'KRT14',
            'KRT34', 'KRTAP1-5', 'KRTAP2-3', 'LMCD1', 'PPP1R14A', 'SERPINB2',
            'SYNPO2', 'TINAGL1',
            ]
        dsF.plot.dot_plot(group_by='cluster', plot_list=genes, layout='vertical')

        dsFa = dsF.average('samples', 'sample', phenotypes=dsF.antibodies)
        dsFa.antibodies = dsF.antibodies
        fig, ax = plt.subplots()
        df = (dsFa.obs[ds.antibodies].T + 0.1)
        df.plot.scatter(x='iMSC', y='AdMSC', ax=ax)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.grid(True)
        fig.tight_layout()
        


    plt.ion()
    plt.show()
