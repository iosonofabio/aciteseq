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

    barcoded = {
        key: dsi.obs['barcode'].str.split('-').str.get(0) for key, dsi in dsp.items()
    }

    for key, bcs in barcoded.items():
        bcs.to_csv(
            f'../../data/sequencing/experiments/201016_A00152_0317_AHMTN3DSXY/barcodes_{key}.csv',
            index=False, header=None,
            )

