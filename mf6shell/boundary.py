#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Tom van Steijn, Royal HaskoningDHV

from scipy.ndimage.morphology import binary_erosion
import numpy as np


def get_square_boundary(nrow, ncol):
    is_boundary = np.zeros((nrow, ncol), dtype=np.bool)
    is_boundary[0, :] = True
    is_boundary[-1, :] = True
    is_boundary[:, 0] = True
    is_boundary[:, -1] = True
    
    return np.where(is_boundary)


def get_idomain_boundary(idomain):
    structure = np.ones((3, 3), dtype=np.bool)
    is_boundary = (idomain == 1) & ~binary_erosion(idomain, structure)
    return np.where(is_boundary)