#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Tom van Steijn, Royal HaskoningDHV

import mf6shell as mfs

from pathlib import Path
import argparse
import logging
import yaml
import os

log = logging.getLogger(os.path.basename(__file__))


def get_parser():
    '''get argumentparser and add arguments'''
    parser = argparse.ArgumentParser(
        'run MF6 model specified in yaml file',
        )

    # Command line arguments
    parser.add_argument('inputfile', type=str,
        help=('YAML input file containing keyword arguments'))
    return parser


def run(**kwargs) -> None:
    # unpack input from kwargs
    name = kwargs['name']
    workspace = Path(kwargs['workspace'])
    exe_name = Path(kwargs['exe_name'])
    grid = kwargs['grid']
    nlay = kwargs['nlay']
    datasources = kwargs['datasources']

    # create output directory if it does not exist
    workspace.mkdir(exist_ok=True)

    # create grid object
    grid = mfs.grid.RegularGrid(**grid)

    # read parameters from datasources
    parameters = mfs.parameters.parameters_from_datasources(datasources, nlay)

    # create basemodel object
    basemodel = mfs.model.Quadi3DModel(name, grid, nlay, parameters)

    # define modflow 6 model
    mf6model = mfs.modflow.MF6Model(basemodel, workspace, exe_name)

    # create external data files
    mf6model.create_external_data_files()

    # create modflow 6 simulation packages
    mf6model.create_simulation_packages()

    # create modflow 6
    mf6model.create_model_packages()

    # write packages
    mf6model.write_packages()

    # run model
    mf6model.run()    


def main():
    # arguments from input file
    args = get_parser().parse_args()
    with open(args.inputfile) as y:
        kwargs = yaml.load(y)
        kwargs['inputfile'] = args.inputfile
    run(**kwargs)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()