#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Tom van Steijn, Royal HaskoningDHV


class Model(object):
    def __init__(self, name, grid):
        self.name = name
        self.grid = grid


class Quasi3DModel(Model):
    def __init__(self, name, grid, nlay, parameters):
        super().__init__(name, grid)
        self.nlay = nlay
        self.parameters = parameters

    @property
    def nlay3d(self):
        return self.nlay*2 - 1
    