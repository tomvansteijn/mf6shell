#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Tom van Steijn, Royal HaskoningDHV


import pickle


class Grid(object):
    def __init__(self):
        pass


class RegularGrid(Grid):
    def __init__(self,
        nrow = 1,
        ncol = 1,
        delr = 1.0,
        delc = 1.0,
        xmin = 0.0,
        ymin = 0.0,
        ):
        self.nrow = nrow
        self.ncol = ncol
        self.delr = delr
        self.delc = delc
        self.xmin = xmin
        self.ymin = ymin

    @property
    def size(self):
        return self.nrow * self.ncol

    @property
    def shape(self):
        return self.nrow, self.ncol

    @shape.setter
    def shape(self, value):
        self.nrow, self.ncol = value

    @property
    def width(self):
        return self.ncol * self.delc

    @property
    def height(self):
        return self.nrow * self.delr

    @property
    def xmax(self):
        return self.xmin + self.width

    @xmax.setter
    def xmax(self, value):
        self.ncol = int((value - self.xmin) / self.delc)

    @property
    def ymax(self):
        return self.ymin + self.height

    @ymax.setter
    def ymax(self, value):
        self.nrow = int((value - self.ymin) / self.delr)

    @property
    def extent(self):
        return self.xmin, self.xmax, self.ymin, self.ymax

    @extent.setter
    def extent(self, value):
        self.xmin, self.xmax, self.ymin, self.ymax = value

    @property
    def transform(self):
        """Affine transformation consistent with GDAL"""
        return self.xmin, self.delc, 0., self.ymax, 0., -self.delr


class PolygonNode(object):
    def __init__(self, number, x, y, verticenumbers):
        self.number = number
        self.x = x
        self.y = y
        self.verticenumbers = verticenumbers

    def to_cell2d(self):
        return [
            self.number,
            self.x,
            self.y,
            len(self.verticenumbers),
            ] + self.verticenumbers
            


class PolygonVertex(object):
    def __init__(self, number, x, y):
        self.number = number
        self.x = x
        self.y = y

    def to_vertex(self):
        return (
            self.number,
            self.x,
            self.y,
            )


class PolygonGrid(Grid):
    def __init__(self,
            nodes,
            vertices,
            nia,
            boundary_nodes,
            river_nodes,
            source_nodes,
            ):
        self.nodes = [PolygonNode(*n) for n in nodes]
        self.vertices = [PolygonVertex(*v) for v in vertices]
        self.nia = nia
        self.boundary_nodes = boundary_nodes
        self.river_nodes = river_nodes
        self.source_nodes = source_nodes

    @staticmethod
    def from_pickle(picklefile):
        with open(picklefile, 'rb') as f:
            return pickle.load(f)

    def to_pickle(self, picklefile):
        with open(picklefile, 'wb') as f:
            pickle.dump(self, f)
