#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Tom van Steijn, Royal HaskoningDHV


import rasterio

import logging
import os

log = logging.getLogger(os.path.basename(__file__))


def read_raster(rasterfile, masked=True,  band=1):
    log.debug('reading {f.name:}'.format(f=rasterfile))
    with rasterio.open(rasterfile) as src:
        return src.read(band, masked=masked)


def write_raster(rasterfile, values, profile):
    log.debug('writing {f.name:}'.format(f=rasterfile))
    with rasterio.open(rasterfile, 'w', **profile) as dst:
        return dst.write(values, 1)