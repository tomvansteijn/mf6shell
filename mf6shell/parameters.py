#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Tom van Steijn, Royal HaskoningDHV

from pathlib import Path
import logging
import os

log = logging.getLogger(os.path.basename(__file__))


def parameters_from_datasources(datasources, nlay) -> Iterable[Parameter]:
    for datasource in datasources:
        if datasource['layered']:
             for ilay in range(nlay):
                pass



class Parameter(object):
    def __init__(self, name, file=None, value=None, layered=False):
        self.name = name
        self.file = Path(file) if file is not None else None
        self.value = value
        self.layered = layered