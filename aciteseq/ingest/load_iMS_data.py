# vim: fdm=indent
'''
author:     Fabio Zanini
date:       18/02/21
content:    Tools to load the data properly
'''
import os
import sys

sys.path.append('/home/fabio/university/postdoc/singlet')
import singlet

from .demux_sex import add_sex
from .antibodies import add_antibodies


def load_data(exclude_low_quality=True, sex=('M', 'F')):

    print('Load data from file')
    fn = '../../data/sequencing/experiments/201016_A00152_0317_AHMTN3DSXY/normalized.h5ad'
    ds = singlet.Dataset(
        dataset={'path': fn},
        )
    ds.counts *= 0.01  # cptt
    ds.counts._normalized = 'counts per ten thousand'
    ds.obs['n_genes'] = (ds.counts > 0).sum(axis=0)

    if exclude_low_quality:
        ds.query_samples_by_metadata('coverage >= 1e4', inplace=True)

    print('Add sex')
    add_sex(ds, inplace=True)

    if set(sex) != set(['M', 'F', 'doublet', 'miss']):
        ds.query_samples_by_metadata(
                'sex in @sex', local_dict=locals(), inplace=True)

    print('Add antibodies')
    add_antibodies(ds, inplace=True)

    return ds


if __name__ == '__main__':

    ds = load_data()
