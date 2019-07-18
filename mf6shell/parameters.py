#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Tom van Steijn, Royal HaskoningDHV

from mf6shell.mixins import AsDictMixin, CopyMixin, ReprMixin
from mf6shell.csvfiles import read_topsys, read_wells
from mf6shell.rasterfiles import read_array

from typing import Iterable
from pathlib import Path
from enum import Enum
import logging
import os

log = logging.getLogger(os.path.basename(__file__))


class ParameterFormat(Enum):
    constant = 1
    raster = 2
    table = 3


class Parameter(AsDictMixin, CopyMixin, ReprMixin):
    def __init__(self,
        name,
        format_,
        file=None,
        value=None,
        ilay=None,
        layered=False,
        ) -> None:
        self.name = name
        self.format = ParameterFormat[format_]
        self.file = Path(file) if file is not None else None
        self.value = value
        self.ilay = ilay        
        self.layered = layered

    def __str__(self):
        return ('{s.__class__.__name__:}('
            'name={s.name:}, '
            'format={s.format:}, '
            'ilay={s.ilay:}'
            ')').format(s=self)

    def is_constant(self):
        return self.format is ParameterFormat.constant

    def expand(self, nlay) -> Iterable['Parameter']:
        for ilay in range(nlay):

            # get a fresh copy
            iparameter = self.copy()

            # format filename using layer number
            iparameter.file = (
                iparameter.file.parent / iparameter.file.name.format(
                    ilay=ilay + 1,
                    )
                )

            # set ilay
            iparameter.ilay = ilay + 1

            # set layered to false
            iparameter.layered = False

            # yield fresh parameter
            yield iparameter

    def get_data(self):
        if self.format is ParameterFormat.raster:
            return read_array(self.file)
        elif self.format is ParameterFormat.table:
            if self.name == 'topsys':
                return read_topsys(self.file)
            elif self.name == 'wells':
                return read_wells(self.file)
            else:
                pass
        else:
            pass

            