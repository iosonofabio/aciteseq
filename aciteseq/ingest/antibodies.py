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


def get_antibody_list():
    sample = 'iMSC'
    fdn_ab = f'../../data/sequencing/raw/201116_A00152_0332_AHNCV5DSXY/{sample}/umi_count'
    fn_abs = f'{fdn_ab}/features.tsv.gz'
    features = pd.read_csv(fn_abs, header=None, sep='\t')[0]
    return list(features)


def get_antibody_controls():
    features = get_antibody_list()
    ctrls = {
        'isotype': [x for x in features if 'sotype Ctrl' in x],
        }
    ctrls_substrings = {
        'positive': [
            'CD29',
            'CD73',
            'CD90',
            'CD105',
            'CD106',
            'CD140a',
            'CD140b',
            'CD146',
            'CD271',
            'Podoplanin',
            'CDH11',
            ],
        'negative': [
            'CD3',
            'CD4',
            'CD8a',
            'CD19',
            'IgG',
            'TCR',
            ],
        }
    for key, substrings in ctrls_substrings.items():
        tmp = []
        for sub in substrings:
            for abname in features:
                if (sub+';') in abname or (sub+' ') in abname:
                    break
            else:
                raise ValueError(f'Control not found: {sub}')
            tmp.append(abname)
        ctrls[key] = tmp

    return ctrls


def get_antibody_counts(sample, barcodes):
    from scipy.io import mmread

    fdn_ab = f'../../data/sequencing/raw/201116_A00152_0332_AHNCV5DSXY/{sample}/umi_count'
    fn_mat = f'{fdn_ab}/matrix.mtx.gz'
    fn_bcs = f'{fdn_ab}/barcodes.tsv.gz'
    fn_abs = f'{fdn_ab}/features.tsv.gz'

    bcs = pd.read_csv(fn_bcs, header=None, sep='\t')[0]
    features = pd.read_csv(fn_abs, header=None, sep='\t')[0]
    mat = np.asarray(mmread(fn_mat).todense())

    df = pd.DataFrame(
        mat.T,
        index=bcs,
        columns=features,
    )

    # Strip the dash
    barcodes = pd.Index(np.array(barcodes))
    if '-' in barcodes[0]:
        barcodes_new = barcodes.str.split('-').str.get(0)
    else:
        barcodes_new = barcodes
    barcodes_new = np.array(barcodes_new)

    # Zero-pad missing cells and sort as input
    missing = set(barcodes_new) - set(bcs)
    for bc in missing:
        df.loc[bc] = 0
    df = df.loc[barcodes_new]

    # Reinstate the dash
    df.index = pd.Index(barcodes, name='CellID')
    df.columns.name = 'Antibody'

    return df


def add_antibodies(ds, inplace=False, normalize=False):
    samples = ds.obs['sample'].unique()
    res = []
    for sample in samples:
        bcs = ds.obs.loc[ds.obs['sample'] == sample, 'barcode']
        cids = bcs.index
        resi = get_antibody_counts(sample, bcs)
        resi.index = cids
        res.append(resi)
    res = pd.concat(res)
    res = res.loc[ds.samplenames]

    if not inplace:
        return res

    ds.obs['coverage_antibodies'] = cov = res.sum(axis=1)
    for col in res.columns:
        if normalize:
            ds.obs[col] = 1.0 * res[col] / cov.values
        else:
            ds.obs[col] = res[col]

    ds.antibodies = list(res.columns)



