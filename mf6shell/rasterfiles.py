#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Tom van Steijn, Royal HaskoningDHV


import rasterio

import logging
import os

log = logging.getLogger(os.path.basename(__file__))


def read_profile(rasterfile):
    with rasterio.open(rasterfile) as src:
        return src.profile


def read_array(rasterfile, masked=True,  band=1):
    with rasterio.open(rasterfile) as src:
        return src.read(band, masked=masked)


def write_array(rasterfile, values, profile):
    with rasterio.open(rasterfile, 'w', **profile) as dst:
        return dst.write(values, 1)