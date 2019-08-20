#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Tom van Steijn, Royal HaskoningDHV

from mf6shell.adofiles import polygongrid_from_teo

import numpy as np
import pytest

import shutil
import os

@pytest.fixture
def teofile(tmpdir):
    datadir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
    teofilename = r'grid.teo'
    teofile = os.path.join(datadir, teofilename)
    testfile = tmpdir.join(teofilename)
    shutil.copyfile(teofile, testfile)
    return testfile


def test_polygongrid_from_teo(self, teofile):
    grid = polygongrid_from_teo(teofile)


def main():
    grid = polygongrid_from_teo(r'data\grid.teo', close_cell_polygons=True)

    import flopy
    from pathlib import Path

    # name
    name = 'test_disv'

    # workspace
    workspace = Path(r'test_disv')

    # exe name
    mf_exe = Path(r'..\bin\mf6.0.4\bin\mf6.exe')

    # create workspace directory
    workspace.mkdir(exist_ok=True)

    # Create the Flopy simulation object
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        exe_name=str(mf_exe), 
        version='mf6',
        sim_ws=str(workspace))

    # Create the Flopy temporal discretization object
    tdis = flopy.mf6.modflow.mftdis.ModflowTdis(sim,
        pname='tdis',
        time_units='DAYS',
        nper=1, 
        perioddata=[(1.0, 1, 1.0)],
        )

    # Create the Flopy groundwater flow (gwf) model object
    model_nam_file = '{}.nam'.format(name)
    gwf = flopy.mf6.ModflowGwf(sim,
        modelname=name, 
        model_nam_file=model_nam_file,
        save_flows=True,
        )

    cell2d = [n.to_cell2d() for n in grid.nodes]
    vertices = [v.to_vertex() for v in grid.vertices]    

    disv = flopy.mf6.modflow.mfgwfdisv.ModflowGwfdisv(gwf,
        pname='disv', nlay=1,
        ncpl=len(cell2d), nvert=len(vertices),
        vertices=vertices, cell2d=cell2d,
        top=10., botm=0.,
        idomain=1,
        length_units='METERS',
        )

    disv.write()

    disv.export(str(workspace / 'disv.shp'))


if __name__ == '__main__':
    import logging
    logging.basicConfig(level=logging.DEBUG)
    main()

