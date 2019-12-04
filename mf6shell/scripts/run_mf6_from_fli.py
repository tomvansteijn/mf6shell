#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Tom van Steijn, Royal HaskoningDHV

from mf6shell.adofiles import (
    get_modflow_name, get_layer_from_name, voronoi_from_teo
    )
from mf6shell.flifiles import FliFile
from mf6shell.grid import PolygonGrid
from mf6shell.parameters import Parameter
from mf6shell.model import Quasi3DModel
from mf6shell.modflow import Modflow6DisvModelWriter
from mf6shell.export import heads_to_ado

from collections import ChainMap
from pathlib import Path
import argparse
import logging
import yaml
import os

log = logging.getLogger(os.path.basename(__file__))


def get_parser():
    '''get argumentparser and add arguments'''
    parser = argparse.ArgumentParser(
        'run MF6 model specified in triwaco fli and teo file',
        )

    # Command line arguments
    parser.add_argument('gridfile', type=str,
        help=('grid definition in gridfile'))
    parser.add_argument('flifile', type=str,
        help=('parameter listing in fli file'))
    parser.add_argument('workspace', type=str,
        help=('model workspace folder path'))
    parser.add_argument('exe_name', type=str, 
        help=('filepath Modflow 6 executable'))
    parser.add_argument('-w', '--write', action='store_true',
        help=('write data files and packages'))
    parser.add_argument('-r', '--run', action='store_true',
        help=('run model'))
    parser.add_argument('-e', '--export', action='store_true',
        help=('export heads'))
    return parser


def run(**kwargs):
    # unpack input from kwargs
    flifile = Path(kwargs['flifile'])
    gridfile = Path(kwargs['gridfile'])
    workspace = Path(kwargs['workspace'])
    exe_name = Path(kwargs['exe_name'])
    config = kwargs['config']
    run_options = {
        'write': kwargs['write'],
        'run': kwargs['run'],
        'export': kwargs['export'],
        }

    # create output directory if it does not exist
    workspace.mkdir(exist_ok=True)

    # read grid from pickle or teo file
    if gridfile.suffix == '.pickle':
        grid = PolygonGrid.from_pickle(gridfile)        
    elif gridfile.suffix == '.teo':
        grid = voronoi_from_teo(gridfile)
        grid.to_pickle(gridfile.with_suffix('.pickle'))
    else:
        raise ValueError(
            'unknown gridfile extension \'{f.suffix}\''.format(
                f=gridfile,
                )
            )

    # read fli file
    fli = FliFile.from_file(flifile)

    # read parameters and add to parameters dict
    parameters = {}
    for fli_parameter in fli.parameters:
        # take a copy to adjust name, layer and filepath
        parameter = fli_parameter.copy()

        # update name and layer number        
        parameter.name, parameter.layer = get_layer_from_name(parameter.name)
        parameter.name = get_modflow_name(parameter.name)

        # update parameter filepath
        parameter.filepath = (flifile.parent / parameter.filepath).resolve()

        # append to list if layered parameter
        if parameter.layer is not None:
            if parameter.name not in parameters:
                parameters[parameter.name] = []
            parameters[parameter.name].append(parameter)
        else:
            parameters[parameter.name] = parameter


    # sort parameters by layer number
    for parameter_name in parameters:
        try: 
            parameters[parameter_name].sort(key=lambda p: p.layer)
        except AttributeError:
            pass

    # create basemodel object
    basemodel = Quasi3DModel(
        fli.header.name,
        grid,
        fli.settings.npak,
        parameters,
        )

    # update solver options based on fli
    options = config.get('package_options')
    options['ims'].update(
        {
            'inner_maximum': fli.calc_options.itmi,
            'outer_maximum': fli.calc_options.itmo,
            },
        )

    # define modflow 6 model
    mf6model = Modflow6DisvModelWriter(basemodel, workspace, exe_name,
        options=options,
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
        # raise NotImplementedError()
        pass

    # run model
    if run_options.get('run', False):
        mf6model.run_simulation()

    # export
    if run_options.get('export', False):
        adofile = mf6model.workspace / 'heads' / 'mf6.flo'
        heads_to_ado(mf6model.workspace / mf6model.headsfile,
            adofile,
            )


def main():
    # arguments from input file
    args = get_parser().parse_args()
    kwargs = vars(args)

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
    logging.basicConfig(level=logging.DEBUG)
    main()