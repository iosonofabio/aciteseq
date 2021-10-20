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
from aciteseq.ingest.antibodies import get_antibody_list, get_antibody_controls


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

    features = get_antibody_list()
    ctrls = get_antibody_controls()
