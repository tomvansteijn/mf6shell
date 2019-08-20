#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Tom van Steijn, Royal HaskoningDHV

from mf6shell.boundary import get_square_boundary, get_idomain_boundary
from mf6shell.datafiles import save_float, save_int, save_array
from mf6shell.wells import get_wells_within_grid

import numpy as np
import flopy

from typing import Tuple, Dict
from pathlib import Path
import logging
import os


log = logging.getLogger(os.path.basename(__file__))


class Modflow6ModelWriter(object):
    def __init__(self,
        model,
        workspace,
        exe_name,
        datafolder=None,
        packages=None,
        data=None,
        maxbound=None,
        options=None,
        namefile=None,
        headsfile=None,
        budgetfile=None,
        ) -> None:

        self.model = model
        self.workspace = Path(workspace)
        self.exe_name = Path(exe_name)

        if datafolder is not None:
            self.datafolder = Path(datafolder)
        else:
            self.datafolder = self.workspace / 'data'

        self.packages = packages or {}
        self.data = data or {}
        self.maxbound = maxbound or {}
        self.options = options or {}

        # input and output filenames
        self.namefile = (
            namefile or '{name:}.nam'.format(name=self.model.name)
            )
        self.headsfile = (
            headsfile or '{name:}.hds'.format(name=self.model.name)
            )
        self.budgetfile = (
            budgetfile or '{name:}.cbc'.format(name=self.model.name)
            )

    def _get_file_ext(self, datfile):
        return {'filename': str(datfile.relative_to(self.workspace))}

    def write_top(self, layer=0, fill_value=0.) -> None:
        # read from parameter
        top = self.model.parameters['top'][layer].get_value()

        # fill masked with fill value
        top = top.filled(fill_value)

        # write to data file
        topdatfile = self.datafolder / 'top.dat'
        save_float(topdatfile, top)

        # add file reference to self.data
        self.data['top'] = self._get_file_ext(topdatfile)

    def write_botm(self, fill_value=1e-6) -> None:
        self.data['botm'] = []
        for layer in range(self.model.nlay):
            botm = (
                self.model.parameters['bot'][layer].get_value()
                )

            # fill masked with fill value
            botm = botm.filled((2*layer + 1) * -fill_value)
            botmdatfile = self.datafolder / 'botm_l{layer:02d}.dat'.format(
                layer=layer*2 + 1,
                )
            save_float(botmdatfile, botm)
            self.data['botm'].append(self._get_file_ext(botmdatfile))
            if (layer + 1) < self.model.nlay:
                botm = (
                    self.model.parameters['top'][layer + 1].get_value()
                    )

                # fill masked with fill value
                botm = botm.filled((2*layer + 2) * -fill_value)

                botmdatfile = self.datafolder / 'botm_l{layer:02d}.dat'.format(
                    layer=layer*2 + 2,
                    )
                save_float(botmdatfile, botm)
                self.data['botm'].append(self._get_file_ext(botmdatfile))

    def write_idomain(self) -> None:
        if self.model.parameters['idomain'].is_constant():
            self.data['idomain'] = self.model.parameters['idomain'].value
        else:
            # read from parameter
            idomain = self.model.parameters['idomain'].get_value()

            # fill masked with fill value
            idomain = idomain.filled(fill_value)

            # write to data file
            idomaindatfile = self.datafolder / 'idomain.dat'
            save_int(idomaindatfile, idomain)

            # add file reference to self.data
            self.data['idomain'] = self._get_file_ext(idomaindatfile)

    def write_kh(self, fill_value=1e-6) -> None:
        self.data['kh'] = []
        for layer in range(self.model.nlay):
            kd = self.model.parameters['kd'][layer].get_value()

            # convert to kh
            top = self.model.parameters['top'][layer].get_value()
            bot = self.model.parameters['bot'][layer].get_value()
            kh = kd / (top - bot)

            # fill with low value
            kh = kh.filled(fill_value)

            khdatfile = self.datafolder / 'kh_l{layer:02d}.dat'.format(
                layer=layer*2 + 1,
                )
            save_float(khdatfile, kh)
            self.data['kh'].append(self._get_file_ext(khdatfile))

            if (layer + 1) < self.model.nlay:
                # dummy values
                self.data['kh'].append(fill_value)

    def write_kv(self, fill_value=1e6) -> None:
        self.data['kv'] = []
        for layer in range(self.model.nlay):
            # dummy values
            self.data['kv'].append(fill_value)

            if (layer + 1) < self.model.nlay:
                c = self.model.parameters['c'][layer].get_value()

                # convert to kv            
                bot = self.model.parameters['bot'][layer].get_value()
                top = self.model.parameters['top'][layer + 1].get_value()
                kv = (bot - top) / c

                # fill with high value
                kv = kv.filled(fill_value)

                kvdatfile = self.datafolder / 'kv_l{layer:02d}.dat'.format(
                    layer=layer*2 + 2,
                    )
                save_float(kvdatfile, kv)
                self.data['kv'].append(self._get_file_ext(kvdatfile))

    def write_start(self, fill_value=0.) -> None:
        self.data['start'] = []
        for layer in range(self.model.nlay):
            start = self.model.parameters['start'][layer].get_value()

            # fill with fill_value
            start = start.filled(fill_value)

            startdatfile = self.datafolder / 'start_l{layer:02d}.dat'.format(
                layer=layer*2 + 1,
                )
            save_float(startdatfile, start)
            self.data['start'].append(self._get_file_ext(startdatfile))
            if (layer + 1) < self.model.nlay:
                self.data['start'].append(self._get_file_ext(startdatfile))

    def write_recharge(self, fill_value=0.) -> None:
        # read from parameter
        recharge = self.model.parameters['recharge'].get_value()

        # fill masked with fill value
        recharge = recharge.filled(fill_value)

        # write to data file
        rechargedatfile = self.datafolder / 'recharge.dat'
        save_float(rechargedatfile, recharge)

        # add file reference to self.data
        self.data['recharge'] = self._get_file_ext(rechargedatfile)

    def write_chd(self, fill_value=0.) -> None:
        if self.model.parameters['idomain'].is_constant():
            chd_row, chd_col = get_square_boundary(*self.model.grid.shape)
        else:
            chd_row, chd_col = get_idomain_boundary(idomain)

        chd_data = []
        for layer in range(self.model.nlay3d):
            start = self.model.parameters['start'][layer//2].get_value()

            # fill with fill_value
            start = start.filled(fill_value)

            # get boundary values
            chd_values = start[chd_row, chd_col]

            # expand layer to vector
            layer_c = np.ones_like(chd_row, dtype=np.int) * layer

            # combine in record array
            chd_array = np.rec.fromarrays(
                [layer_c + 1, chd_row + 1, chd_col + 1, chd_values],
                names=['layer', 'row', 'col', 'value'],
                )

            # append
            chd_data.append(chd_array)

        # concatenate
        chd_data = np.concatenate(chd_data)

        # write to data file
        chddatfile = self.datafolder / 'chd.dat'
        save_array(chddatfile, chd_data, fmt='  %i %i %i %16.8f')

        # add file reference to self.data
        self.data['chd'] = self._get_file_ext(chddatfile)

        # add nrows to maxbound
        self.maxbound['chd'] = len(chd_data)

    def write_topsys(self, fill_value=0.):
        topsys_data = self.model.parameters['topsys'].get_value()

        # drain package
        drn_columns = ['layer', 'row', 'col', 'head', 'cond']

        # select drain data 
        select_drn = topsys_data.loc[:, 'inffact'] < 1.

        # return if no selection
        if select_drn.sum() == 0:
            return

        # select and copy rows
        drn_data = topsys_data.loc[select_drn, :].copy()

        # update conductivity with infiltration factor
        drn_data.loc[:, 'cond'] = (
            drn_data.loc[:, 'cond'] * (1. - drn_data.loc[:, 'inffact'])
            )

        # clip conductivity at zero
        drn_data.loc[:, 'cond'] = drn_data.loc[:, 'cond'].clip(lower=0.)

        # convert layer number to full3d (1 -> 1, 2 -> 3, 3 -> 5, etc.)
        drn_data.loc[:, 'layer'] = drn_data.loc[:, 'layer'] * 2 - 1

        # write to data file
        drndatfile = self.datafolder / 'drn.dat'
        save_array(drndatfile, drn_data.loc[:, drn_columns],
            fmt='  %i %i %i %16.8f %16.8f',
            )

        # add file reference to self.data
        self.data['drn'] = self._get_file_ext(drndatfile)

        # add nrows to maxbound
        self.maxbound['drn'] = len(drn_data)

        # ghb package
        ghb_columns = ['layer', 'row', 'col', 'head', 'cond']
        
        # select ghb data 
        select_ghb = (
            topsys_data.loc[:, 'rbot'].isna() &
            (topsys_data.loc[:, 'inffact'] > 0.)
            )

        # return if no selection
        if select_ghb.sum() == 0:
            return

        # select and copy rows
        ghb_data = topsys_data.loc[select_ghb, :].copy()
        
        # update conductivity with infiltration factor
        ghb_data.loc[:, 'cond'] = (
            ghb_data.loc[:, 'cond'] * ghb_data.loc[:, 'inffact']
            )

        # clip conductivity at zero
        ghb_data.loc[:, 'cond'] = ghb_data.loc[:, 'cond'].clip(lower=0.)

        # convert layer number to full3d (1 -> 1, 2 -> 3, 3 -> 5, etc.)
        ghb_data.loc[:, 'layer'] = ghb_data.loc[:, 'layer'] * 2 - 1

        # write to data file
        ghbdatfile = self.datafolder / 'ghb.dat'
        save_array(ghbdatfile, ghb_data.loc[:, ghb_columns],
            fmt='  %i %i %i %16.8f %16.8f',
            )

        # add file reference to self.data
        self.data['ghb'] = self._get_file_ext(ghbdatfile)

        # add nrows to maxbound
        self.maxbound['ghb'] = len(ghb_data)

        # riv package
        riv_columns = ['layer', 'row', 'col', 'head', 'cond', 'rbot']
        
        # select riv data 
        select_riv = (
            topsys_data.loc[:, 'rbot'].notna() &
            (topsys_data.loc[:, 'inffact'] > 0.)
            )

        # return if no selection
        if select_riv.sum() == 0:
            return

        # select and copy rows
        riv_data = topsys_data.loc[select_riv, :].copy()
        
        # update conductivity with infiltration factor
        riv_data.loc[:, 'cond'] = (
            riv_data.loc[:, 'cond'] * riv_data.loc[:, 'inffact']
            )

        # clip conductivity at zero
        riv_data.loc[:, 'cond'] = riv_data.loc[:, 'cond'].clip(lower=0.)

        # convert layer number to full3d (1 -> 1, 2 -> 3, 3 -> 5, etc.)
        riv_data.loc[:, 'layer'] = riv_data.loc[:, 'layer'] * 2 - 1

        # write to data file
        rivdatfile = self.datafolder / 'riv.dat'
        save_array(rivdatfile, riv_data.loc[:, riv_columns],
            fmt='  %i %i %i %16.8f %16.8f %16.8f',
            )

        # add file reference to self.data
        self.data['riv'] = self._get_file_ext(rivdatfile)

        # add nrows to maxbound
        self.maxbound['riv'] = len(riv_data)


    def write_wel(self) -> None:
        wel_data = self.model.parameters['wells'].get_value()

        # get row, column in grid
        wel_data = get_wells_within_grid(wel_data, self.model.grid, 'x', 'y')

        # select layer numbers, row, column & pumping rates
        wel_data = wel_data.loc[:, ['layer', 'row', 'col', 'q_assigned']]

        # convert layer number to full3d (1 -> 1, 2 -> 3, 3 -> 5, etc.)
        wel_data.loc[:, 'layer'] = wel_data.loc[:, 'layer'] * 2 - 1
        wel_data.loc[:, ['row', 'col']] += 1

        # write to data file
        weldatfile = self.datafolder / 'wel.dat'
        save_array(weldatfile, wel_data.values, fmt='  %i %i %i %16.8f')

        # add file reference to self.data
        self.data['wel'] = self._get_file_ext(weldatfile)

        # add nrows to maxbound
        self.maxbound['wel'] = len(wel_data)


    def write_data(self) -> None:
        # create data folder
        self.datafolder.mkdir(exist_ok=True)

        # write data per package
        self.write_top()
        self.write_botm()
        self.write_idomain()
        self.write_kh()
        self.write_kv()
        self.write_start()
        self.write_recharge()
        self.write_chd()
        self.write_topsys()
        self.write_wel()

    def create_simulation_packages(self) -> None:
        # create the Flopy simulation object
        self.packages['sim'] = flopy.mf6.MFSimulation(
            sim_name=self.model.name,
            exe_name=str(self.exe_name), 
            sim_ws=str(self.workspace),
            **self.options.get('sim', {}),
            )

        # create the Flopy temporal discretization object
        self.packages['tdis'] = flopy.mf6.modflow.mftdis.ModflowTdis(
            self.packages['sim'],
            **self.options.get('tdis', {}),
            )

        # create the Flopy groundwater flow (gwf) model object
        self.packages['gwf'] = flopy.mf6.ModflowGwf(
            self.packages['sim'],
            modelname=self.model.name, 
            model_nam_file=self.namefile,
            **self.options.get('gwf', {}),
            )

        # create the Flopy iterative model solver (ims) Package object
        self.packages['ims'] = flopy.mf6.modflow.mfims.ModflowIms(
            self.packages['sim'],
            **self.options.get('ims', {}),
            )

    def create_model_packages(self) -> None:
        # create the DIS package
        self.packages['dis'] = flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
            self.packages['gwf'],
            nlay=self.model.nlay3d,
            nrow=self.model.grid.nrow, ncol=self.model.grid.ncol,
            delr=self.model.grid.delr, delc=self.model.grid.delc,
            xorigin=self.model.grid.xmin, yorigin=self.model.grid.ymin,
            top=self.data['top'], botm=self.data['botm'],
            idomain=self.data['idomain'],
            **self.options.get('dis', {}),
            )

        # create the NPF package
        self.packages['npf'] = flopy.mf6.modflow.mfgwfnpf.ModflowGwfnpf(
            model=self.packages['gwf'],
            k=self.data['kh'],
            k22=self.data['kh'],
            k33=self.data['kv'],
            **self.options.get('npf', {}),
            )

        # create the initial conditions package
        self.packages['ic'] = flopy.mf6.modflow.mfgwfic.ModflowGwfic(
            model=self.packages['gwf'], 
            strt=self.data['start'],
            **self.options.get('ic', {}),
            )

        # create the RCH package
        self.packages['rch'] = flopy.mf6.ModflowGwfrcha(
            model=self.packages['gwf'], 
            recharge=[self.data['recharge'],],
            **self.options.get('rch', {}),
            )

        # create the CHD package
        self.packages['chd'] = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(
            model=self.packages['gwf'], 
            stress_period_data={0: self.data['chd']},
            maxbound=self.maxbound['chd'],
            **self.options.get('chd', {}),
            )

        # create the DRN package
        if 'drn' in self.data:
            self.packages['drn'] = flopy.mf6.modflow.mfgwfdrn.ModflowGwfdrn(
                self.packages['gwf'],
                stress_period_data={0: self.data['drn']},
                maxbound=self.maxbound['drn'],
                **self.options.get('drn', {}),
                )

        # create the GHB package
        if 'ghb' in self.data:
            self.packages['ghb'] = flopy.mf6.modflow.mfgwfghb.ModflowGwfghb(
                self.packages['gwf'],
                stress_period_data={0: self.data['ghb']},
                maxbound=self.maxbound['ghb'],
                **self.options.get('ghb', {}),
                )

        # create the RIV package
        if 'riv' in self.data:
            self.packages['riv'] = flopy.mf6.modflow.mfgwfriv.ModflowGwfriv(
                self.packages['gwf'],
                stress_period_data={0: self.data['riv']},
                maxbound=self.maxbound['riv'],
                **self.options.get('riv', {}),
                )

        # create WEL package
        if 'wel' in self.data:
            self.packages['wel'] = flopy.mf6.modflow.mfgwfwel.ModflowGwfwel(
                self.packages['gwf'],
                stress_period_data={0: self.data['wel']},
                maxbound=self.maxbound['wel'],               
                **self.options.get('wel', {}),
                )

        # create the output control package
        self.packages['oc'] = flopy.mf6.modflow.mfgwfoc.ModflowGwfoc(
            self.packages['gwf'],            
            head_filerecord=[self.headsfile, ],
            budget_filerecord=[self.budgetfile, ],
            **self.options.get('oc', {}),
            )

    def write_simulation(self) -> None:
        self.packages['sim'].write_simulation()

    def run_simulation(self) -> Tuple[bool, list]:
        return self.packages['sim'].run_simulation()



