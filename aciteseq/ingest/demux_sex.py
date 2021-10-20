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
        return sex

    ds.obs['sex'] = sex
    ds.obs['sample_unique'] = ds.obs['sample'].astype(str) + '_' + ds.obs['sex']
