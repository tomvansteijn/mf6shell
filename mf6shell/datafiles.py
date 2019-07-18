#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Tom van Steijn, Royal HaskoningDHV

import numpy as np

import logging
import os

log = logging.getLogger(os.path.basename(__file__))


def save_int(datfile, array):
    save_array(datfile, array, fmt='%i', delimiter=' ')


def save_float(datfile, array):
    save_array(datfile, array, fmt='%15.6E')


def save_array(datfile, array, fmt, delimiter=''):
    np.savetxt(datfile, array, fmt=fmt, delimiter=delimiter)
