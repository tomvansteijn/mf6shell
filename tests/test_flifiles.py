#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Tom van Steijn, Royal HaskoningDHV

from mf6shell.flifiles import FliFile

import numpy as np
import pytest

import shutil
import os

@pytest.fixture
def flifile(tmpdir):
    datadir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
    flifilename = r'flairs.fli'
    flifile = os.path.join(datadir, flifilename)
    testfile = tmpdir.join(flifilename)
    shutil.copyfile(flifile, testfile)
    return testfile


class TestFliFile(object):
    def test_fli_from_file(self, flifile):
        fli = FliFile.from_file(flifile)
        fli