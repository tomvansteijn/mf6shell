#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Tom van Steijn, Royal HaskoningDHV

from mf6shell.grid import RegularGrid


class TestGrid(object):
    def test_size(self):
        g = RegularGrid(nrow=6, ncol=7)
        assert g.size == 42
        
    def test_shape(self):
        g = RegularGrid(nrow=6, ncol=7)
        assert g.shape == (6, 7)

    def test_width(self):
        g = RegularGrid(nrow=6, ncol=7, delc=20.)
        assert g.width == 140.

    def test_height(self):
        g = RegularGrid(nrow=6, ncol=7, delr=50.)
        assert g.height == 300.

    def test_xmax(self):
        g = RegularGrid(nrow=6, delc=5., ncol=7, xmin=10.)
        assert g.xmax == 45.

    def test_set_xmax(self):
        g = RegularGrid(nrow=6, delc=5., ncol=7, xmin=10.)
        g.xmax = 30
        assert g.ncol == 4

    def test_ymax(self):
        g = RegularGrid(nrow=6, delr=2., ncol=7, ymin=67.)
        assert g.ymax == 79.

    def test_set_ymax(self):
        g = RegularGrid(nrow=6, delr=2., ncol=7, ymin=67.)
        g.ymax = 107.
        assert g.nrow == 20

    def test_extent(self):
        g = RegularGrid(nrow=6, ncol=7, delr=2., delc=5., xmin=10., ymin=67.)
        assert g.extent == (10., 45., 67., 79.)

    def test_set_extent(self):
        g = RegularGrid(nrow=6, ncol=7, delr=2., delc=5., xmin=10., ymin=67.)
        g.extent = 10., 50., 67., 157.
        assert g.shape == (45, 8)

    def test_transform(self):
        g = RegularGrid(nrow=6, delr=2., delc=5., xmin=10., ymin=67.)
        assert g.transform == (10., 5., 0., 79., 0., -2.,)