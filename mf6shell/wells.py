#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Tom van Steijn, Royal HaskoningDHV

from affine import Affine
import rasterio

def get_wells_within_grid(wel_data, grid, xfield, yfield):
    # copy
    wel_data = wel_data.copy()

    # transform
    fwd = Affine.from_gdal(*grid.transform)

    # transform xy to row,col
    wel_data.loc[:, 'row'], wel_data.loc[:, 'col'] = (
        rasterio.transform.rowcol(
            fwd,
            wel_data.loc[:, xfield],
            wel_data.loc[:, yfield],
            ))

    # select wells in grid
    in_model = (
        (wel_data.loc[:, 'row'] >= 0) & (wel_data.loc[:, 'row'] < grid.nrow) &
        (wel_data.loc[:, 'col'] >= 0) & (wel_data.loc[:, 'col'] < grid.ncol)
        )
    return wel_data.loc[in_model, :]

