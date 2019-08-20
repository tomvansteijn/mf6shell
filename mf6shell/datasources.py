#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Tom van Steijn, Royal HaskoningDHV

from mf6shell.parameters import ConstantParameter, RasterParameter, TableParameter

from enum import Enum


class DataSourceFormat(Enum):
    constant = 1
    raster = 2
    topsys_table = 3
    wells_table = 4


class DataSource(object):
    def __init__(self,
        name,
        fmt,
        value=None,
        filepath=None,
        layer=None,
        layered=False,
        ):
        self.name = name
        self.format = DataSourceFormat[fmt]
        self.value = value
        self.filepath = filepath
        self.layer = layer
        self.layered = layered

    def to_parameter(self, layer=None):
        if layer is None:
            layer = self.layer
        if self.layered and layer is not None:
            filepath_for_parameter = self.filepath.format(
                    layer=layer + 1,
                    )
        else:
            filepath_for_parameter = self.filepath
            
        if self.format is DataSourceFormat.constant:
            return ConstantParameter(
                name=self.name,
                value=self.value,
                layer=layer,
                )
        elif self.format is DataSourceFormat.raster:
            return RasterParameter(
                name=self.name,
                filepath=filepath_for_parameter,
                layer=layer,
                )
        elif self.format is DataSourceFormat.topsys_table:
            return TableParameter(
                name=self.name,
                filepath=filepath_for_parameter,
                schema='topsys',
                layer=layer,
                )
        elif self.format is DataSourceFormat.wells_table:
            return TableParameter(
                name=self.name,
                filepath=filepath_for_parameter,
                schema='wells',
                layer=layer,
                )
        else:
            raise ValueError(
                'unknown datasource format {:d}'.format(self.format.value)
                )


