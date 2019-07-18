#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Tom van Steijn, Royal HaskoningDHV

from mf6shell.rasterfiles import write_array

from affine import Affine
import numpy as np
import flopy

import logging
import os


log = logging.getLogger(os.path.basename(__file__))

def export_heads(headsfile, rasterfolder, transform,
    noflow=1e30,
    driver='GTiff',
    epsg=28992,
    ):
    # create output folder
    rasterfolder.mkdir(exist_ok=True)

    hds = flopy.utils.binaryfile.HeadFile(headsfile)
    heads = np.ma.masked_equal(hds.get_data(), noflow)

    for ilay, iheads in enumerate(heads[::2]):
        log.debug('exporting heads layer {ilay:d}'.format(ilay=ilay + 1))
        rasterfile = rasterfolder / 'heads_l{ilay:02d}.tif'.format(
            ilay=ilay + 1,
            )

        width, height = iheads.shape
        profile = {
            'driver': driver,
            'width': width,
            'height': height,
            'count': 1,
            'transform': Affine.from_gdal(*transform),
            'dtype': iheads.dtype,
            }

        write_array(rasterfile, iheads, profile)



