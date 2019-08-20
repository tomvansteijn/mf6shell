#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Tom van Steijn, Royal HaskoningDHV

import re

HEADERPATTERN = (
    r'^\[(?P<timeunit>[a-z]+).(?P<lengthunit>[a-z]+)\] (?P<name>\w+)'
    )


class FliHeader(object):
    def __init__(self, timeunit, lengthunit, name):
        self.timeunit = timeunit
        self.lengthunit = lengthunit
        self.name = name


class FliSettings(object):
    def __init__(self, npak, iffr, ifss, ifsf, iftr, nsar, rlax):
        self.npak = npak
        self.iffr = iffr
        self.ifss = ifss
        self.ifsf = ifsf
        self.iftr = iftr
        self.nsar = nsar
        self.rlax = rlax


class FliParameter(object):
    def __init__(self, name, filepath, block_name):
        self.name = name
        self.filepath = filepath
        self.block_name = block_name


class FliCalculationOptions(object):
    def __init__(self, itmi, itmo, feps, fepssc, fepshb, precon):
        self.itmi = itmi
        self.itmo = itmo
        self.feps = feps
        self.fepssc = fepssc
        self.fepshb = fepshb
        self.precon = precon


class FliPrintOutput(object):
    def __init__(self, iprt, irst):
        self.iprt = iprt
        self.irst = irst


class FliFile(object):
    def __init__(self,
        header,
        settings,
        parameters,
        calc_options,
        print_output,
        ):
        self.header = FliHeader(**header)
        self.settings = FliSettings(**settings)
        self.parameters = [FliParameter(**p) for p in parameters]
        self.calc_options = FliCalculationOptions(**calc_options)
        self.print_output = FliPrintOutput(**print_output)

    @classmethod
    def from_file(cls, flifile):
        with open(flifile) as f:
            # lines generator
            lines = (l.rstrip('\n') for l in f)

            # read units and name
            header_line = next(lines)
            m = re.search(HEADERPATTERN, header_line)
            header = {
                'timeunit': m.group('timeunit'),
                'lengthunit': m.group('lengthunit'),
                'name': m.group('name'),              
                }

            # read settings
            settings_line = next(lines)
            settings = {
                'npak': int(settings_line[:5]),
                'iffr': int(settings_line[5:10]),
                'ifss': int(settings_line[10:15]),
                'ifsf': int(settings_line[15:20]),
                'iftr': int(settings_line[20:25]),
                'nsar': int(settings_line[25:30]),
                'rlax': float(settings_line[30:]),
                }

            # skip 4 lines
            # end of sum input
            # end of sources input
            # end of boundary input
            # end of river input
            for iline in range(4):
                next(lines)

            # read parameters
            parameters = []
            while True:
                line = next(lines)
                if line == 'end':
                    break
                parameter = {}          
                parameter['name'] = next(lines)
                parameter['filepath'] = next(lines)
                parameter['block_name'] = next(lines)

                parameters.append(
                    parameter
                    )

            # read calculation options
            calc_options_line = next(lines)
            calc_options = {
                'itmi': int(calc_options_line[:5]),
                'itmo': int(calc_options_line[5:10]),
                'feps': float(calc_options_line[10:20]),
                'fepssc': float(calc_options_line[20:30]),
                'fepshb': float(calc_options_line[30:40]),
                'precon': int(calc_options_line[40:]),
                }

            # read print output flags
            print_output_line = next(lines)
            print_output = {
                'iprt': int(print_output_line[:5]),
                'irst': int(print_output_line[5:]),
                }

            return cls(
                header,
                settings,
                parameters,
                calc_options,
                print_output,
                )






