#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Tom van Steijn, Royal HaskoningDHV

from mf6shell.mixins import AsDictMixin, CopyMixin, ReprMixin
from mf6shell.adofiles import read_ado
from mf6shell.csvfiles import TableSchema, read_table
from mf6shell.rasterfiles import read_raster

from typing import Iterable
from pathlib import Path
import logging
import os

log = logging.getLogger(os.path.basename(__file__))


class Parameter(AsDictMixin, CopyMixin, ReprMixin):
    def __init__(self,
        name,
        layer=None,
        ):
        self.name = name
        self.layer = layer

    def __str__(self):
        return ('{s.__class__.__name__:}('
            'name={s.name:}, '
            'layer={s.layer:}'
            ')').format(s=self)

    def is_constant(self):
        return False

    def get_value(self):
        raise NotImplementedError('Not implemented in this class')


class ConstantParameter(Parameter):
    def __init__(self, name, value, layer=None):
        super().__init__(name, layer)
        self.value = value

    def is_constant(self):
        return True

    def get_value(self):
        return self.value


class FileParameter(Parameter):
    def __init__(self, name, filepath, layer=None):
        super().__init__(name, layer)
        self.filepath = Path(filepath)


class RasterParameter(FileParameter):
    def get_value(self):
        return read_raster(self.filepath)


class AdoParameter(FileParameter):
    def __init__(self, name, filepath, block_name, layer=None):
        super().__init__(name, filepath, layer)
        self.block_name = block_name

    def get_value(self):
        return read_ado(self.filepath, self.block_name)


class TableParameter(FileParameter):
    def __init__(self, name, filepath, schema, layer=None):
        super().__init__(name, filepath, layer)
        self.schema = TableSchema[schema]

    def get_value(self):
        return read_table(self.filepath, self.schema)
