#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Tom van Steijn, Royal HaskoningDHV


import pandas as pd

from enum import Enum
import logging
import os


log = logging.getLogger(os.path.basename(__file__))


class TableSchema(Enum):
    topsys = 1
    wells = 2


def read_table(csvfile, schema):
    if schema is TableSchema.topsys:
        return read_topsys(csvfile)
    elif schema is TableSchema.wells:
        return read_wells(csvfile)


def read_topsys(csvfile):
    log.debug('reading {f.name:}'.format(f=csvfile))
    column_names = ['layer', 'row', 'col', 'head', 'rbot', 'cond', 'inffact', 'topsys']
    topsys_data = pd.read_csv(csvfile,
        index_col=None,
        header=None,
        names=column_names,
        delim_whitespace=True,
        na_values=[-999.],
        )
    return topsys_data


def read_wells(csvfile):
    log.debug('reading {f.name:}'.format(f=csvfile))
    wel_data = pd.read_csv(csvfile)
    return wel_data