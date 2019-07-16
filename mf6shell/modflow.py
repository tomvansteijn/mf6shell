#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Tom van Steijn, Royal HaskoningDHV

import flopy

from pathlib import Path
import logging
import os


log = logging.getLogger(os.path.basename(__file__))


class Modflow6Model(object):
    def __init__(self,
        basemodel,
        workspace,
        exe_name,
        packages=None,
        datafiles=None,
        options=None,
        ) -> None:
        self.basemodel = basemodel
        self.workspace = Path(workspace)
        self.exe_name = Path(exe_name)
        self.packages = packages or {}
        self.datafiles = datafiles or {}
        self.options = options or {}

    def write_external_datafiles(self) -> None:
        pass

    def create_simulation_packages(self) -> None:
        # Create the Flopy simulation object
        self.packages['sim'] = flopy.mf6.MFSimulation(
            sim_name=self.basemodel.name,
            exe_name=str(self.exe_name), 
            sim_ws=str(self.workspace),
            **self.options.get('sim', {}),
            )

        # Create the Flopy temporal discretization object
        self.packages['tdis'] = flopy.mf6.modflow.mftdis.ModflowTdis(
            self.packages['sim'],
            **self.options.get('tdis', {}),
            )

        # Create the Flopy groundwater flow (gwf) model object
        self.packages['gwf'] = flopy.mf6.ModflowGwf(
            self.packages['sim'],
            modelname=self.basemodel.name, 
            model_nam_file='{name:}.nam'.format(name=self.basemodel.name),
            **self.options.get('gwf', {}),
            )

        # Create the Flopy iterative model solver (ims) Package object
        self.packages['ims'] = flopy.mf6.modflow.mfims.ModflowIms(
            self.packages['sim'],
            **self.options.get('ims', {}),
            )

    def create_model_packages(self) -> None:
        # Create the DIS package
        self.packages['dis'] = flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
            self.packages['gwf'],
            nlay=self.basemodel.nlay3d,
            nrow=self.basemodel.grid.nrow, ncol=self.basemodel.grid.ncol,
            delr=self.basemodel.grid.delr, delc=self.basemodel.grid.delc,
            top=top_ext, botm=botm_ext,
            idomain=idomain_ext,
            **self.options.get('dis', {}),
            )

        # Create the NPF package
        self.packages['npf'] = flopy.mf6.modflow.mfgwfnpf.ModflowGwfnpf(
            model=self.packages['gwf'],
            k=kh_ext,
            k22=kh_ext,
            k33=kv_ext,
            **self.options.get('npf', {}),
            )

        # Create the initial conditions package
        self.packages['ic'] = flopy.mf6.modflow.mfgwfic.ModflowGwfic(
            model=self.packages['gwf'], 
            strt=start_ext,
            **self.options.get('ic', {}),
            )

        # Create the CHD package
        self.packages['chd'] = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(
            model=self.packages['gwf'], 
            maxbound=len(chd_data),
            stress_period_data=chd_ext,
            **self.options.get('dis', {}),
            )

        # initialize the RCH package
        self.packages['rch'] = flopy.mf6.ModflowGwfrcha(
            model=self.packages['gwf'], 
            recharge=recharge_ext,
            **self.options.get('rch', {}),
            )

        # initialize the DRN package
        self.packages['drn'] = flopy.mf6.modflow.mfgwfdrn.ModflowGwfdrn(
            self.packages['gwf'],
            maxbound=drn_data.shape[0],
            stress_period_data=drn_ext,
            **self.options.get('drn', {}),
            )

        # initialize the GHB package
        self.packages['ghb'] = flopy.mf6.modflow.mfgwfghb.ModflowGwfghb(
            self.packages['gwf'],
            maxbound=ghb_data.shape[0],
            stress_period_data=ghb_ext,
            **self.options.get('ghb', {}),
            )

        # initialize the RIV package
        self.packages['riv'] = flopy.mf6.modflow.mfgwfriv.ModflowGwfriv(
            self.packages['gwf'],
            maxbound=riv_data.shape[0],
            stress_period_data=riv_ext,
            **self.options.get('riv', {}),
            )

        # initialize WEL package
        self.packages['wel'] = flopy.mf6.modflow.mfgwfwel.ModflowGwfwel(
            self.packages['gwf'],
            maxbound=len(wel_data),
            stress_period_data=wel_ext,
            **self.options.get('wel', {}),
            )

        # Create the output control package
        headfile = '{}.hds'.format(self.basemodel.name)
        head_filerecord = [headfile]
        budgetfile = '{}.cbb'.format(self.basemodel.name)
        budget_filerecord = [budgetfile]
        saverecord = [('HEAD', 'ALL'), 
                      ('BUDGET', 'ALL')]
        printrecord = [('HEAD', 'LAST')]
        self.packages['oc'] = flopy.mf6.modflow.mfgwfoc.ModflowGwfoc(
            self.packages['gwf'],            
            head_filerecord=head_filerecord,
            budget_filerecord=budget_filerecord,
            saverecord=saverecord, 
            printrecord=printrecord,
            **self.options.get('oc', {}),
            )

    def write_simulation(self) -> None:
        self.packages['sim'].write_simulation()

    def run_model(self) -> bool, list:
        return self.packages['sim'].run_model()



