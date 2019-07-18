#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Tom van Steijn, Royal HaskoningDHV


import pandas as pd


def read_topsys(csvfile):
    column_names = ['ilay', 'row', 'col', 'head', 'rbot', 'cond', 'inffact', 'topsys']
    topsys_data = pd.read_csv(csvfile,
        index_col=None,
        header=None,
        names=column_names,
        delim_whitespace=True,
        na_values=[-999.],
        )
    return topsys_data

def read_wells(csvfile):
    wel_data = pd.read_csv(csvfile)
    return wel_data