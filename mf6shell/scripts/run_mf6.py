#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Tom van Steijn, Royal HaskoningDHV

from mf6shell.grid import RegularGrid
from mf6shell.parameters import Parameter
from mf6shell.model import Quasi3DModel
from mf6shell.modflow import Modflow6Model
from mf6shell.export import export_heads

from collections import ChainMap
from pathlib import Path
import argparse
import logging
import yaml
import os

log = logging.getLogger(os.path.basename(__file__))


def get_parser() -> argparse.ArgumentParser:
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
    run_options = kwargs['run_options']
    config = kwargs['config']

    # create output directory if it does not exist
    workspace.mkdir(exist_ok=True)

    # create grid object
    grid = RegularGrid(**grid)

    # read parameters, expand by layer and add to parameters dict
    parameters = {}
    for datasource in datasources:
        parameter = Parameter(**datasource)
        if parameter.layered:
            parameters[parameter.name] = {}
            for iparameter in parameter.expand(nlay):
                parameters[parameter.name][iparameter.ilay] = iparameter
        else:
            parameters[parameter.name] = parameter

    # create basemodel object
    basemodel = Quasi3DModel(name, grid, nlay, parameters)

    # define modflow 6 model
    mf6model = Modflow6Model(basemodel, workspace, exe_name,
        options=config.get('package_options'),
        )

    if run_options.get('write', False):   
        # create external data files
        mf6model.write_data()

        # create modflow 6 simulation packages
        mf6model.create_simulation_packages()

        # create modflow 6
        mf6model.create_model_packages()

        # write packages
        mf6model.write_simulation()
    else:
        # load from name file??
        pass

    # run model
    if run_options.get('run', False):
        mf6model.run_simulation()

    # export
    if run_options.get('export', False):
        exportfolder = mf6model.workspace / 'heads'
        export_heads(mf6model.workspace / mf6model.headsfile,
            exportfolder,
            grid.transform,
            )


def main():
    # arguments from input file
    args = get_parser().parse_args()
    with open(args.inputfile) as y:
        kwargs = yaml.load(y)
        kwargs['inputfile'] = args.inputfile

    # read default config
    scripts_folder = os.path.dirname(os.path.realpath(__file__))
    defaultconfigfile = os.path.join(os.path.dirname(scripts_folder),
        'defaultconfig.yaml')
    with open(defaultconfigfile) as y:
        defaultconfig = yaml.load(y)

    # get user config from input file
    userconfig = kwargs.get('config') or {}

    # chain config
    kwargs['config'] = ChainMap(userconfig, defaultconfig)

    run(**kwargs)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()