# vim: fdm=indent
'''
author:     Fabio Zanini
date:       16/01/21
content:    Load and QC data
'''
import os
import sys
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import hstack

import anndata


if __name__ == '__main__':

    fdn_raw = '../../data/sequencing/raw/201016_A00152_0317_AHMTN3DSXY_summary/'
    samples = ['AdMSC', 'iMSC']

    data_raw = {}
    for sample in samples:
        print(sample)
        fdn = f'{fdn_raw}{sample}/summary/filtered_feature_bc_matrix/'
        fn_bcs = f'{fdn}barcodes.tsv.gz'
        fn_genes = f'{fdn}features.tsv.gz'
        fn_counts = f'{fdn}matrix.mtx.gz'

        print('Read barcodes')
        barcodes = pd.read_csv(fn_bcs, sep='\t', compression='gzip', header=None)
        barcodes.columns = ['barcode']
        barcodes['sample'] = sample
        print('Read features')
        features = pd.read_csv(fn_genes, sep='\t', compression='gzip', header=None)
        print('Read counts')
        counts = mmread(fn_counts)

        data_raw[sample] = {
                'barcodes': barcodes,
                'features': features,
                'counts': counts,
        }

    print('Set coverage')
    for sample, dic in data_raw.items():
        dic['barcodes']['coverage'] = np.asarray(dic['counts'].sum(axis=0))[0]

    print('Filter out zero genes')
    ngenes = dic['features'].shape[0]
    exp_cells = pd.Series(np.zeros(ngenes, np.int64))
    for sample, dic in data_raw.items():
        exp_cells += np.asarray((dic['counts'] > 0).sum(axis=1))[:, 0]
    idx = exp_cells.values.nonzero()[0]

    print('Concatenate')
    counts = hstack([data_raw[x]['counts'].tocsr()[idx] for x in samples])
    barcodes = pd.concat([data_raw[x]['barcodes'] for x in samples])
    barcodes['CellID'] = barcodes['barcode'] + '-' + barcodes['sample']
    barcodes.set_index('CellID', inplace=True)
    features = dic['features'].iloc[idx]

    print('Merge features with the same gene name')
    counts = counts.tolil()
    ngenes = len(features)
    idx = set(list(range(ngenes)))
    genes_dup = features[1].value_counts()
    genes_dup = genes_dup.index[genes_dup > 1]
    for gene in genes_dup:
        irows = (features[1] == gene).values.nonzero()[0]
        i1 = irows[0]
        for i2 in irows[1:]:
            counts[i1] += counts[i2]
            counts[i2] = 0
            idx.remove(i2)
    idx = sorted(idx)
    counts = counts.tocsr()[idx].tocsc().T
    genes = features[1].values[idx]

    print('Make output data structure')
    adata = anndata.AnnData(
        X=counts,
        obs=barcodes,
        )
    adata.var_names = pd.Index(genes, name='GeneName')

    print('Normalize to CPM (and convert to float)')
    import scanpy as sp
    sp.pp.normalize_total(adata, target_sum=1e6)

    print('Write to file')
    fn_out = '../../data/sequencing/experiments/201016_A00152_0317_AHMTN3DSXY/normalized.h5ad'
    adata.write_h5ad(fn_out)
