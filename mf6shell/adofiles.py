#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Tom van Steijn, Royal HaskoningDHV

from mf6shell.grid import PolygonGrid

import numpy as np
import adopy

import logging
import os

log = logging.getLogger(os.path.basename(__file__))


def read_ado(adofile, block_name, masked=True, nodata=-999.):
    log.debug('reading {f.name:}'.format(f=adofile))
    with adopy.open(adofile) as src:
        for block in src.read():
            if block.name == block_name:
                return np.ma.masked_equal(block.values, nodata)


def sort_around_point(xy, xyp, clockwise=False):
    '''sort points around center point'''
    x, y = xy.transpose()
    xp, yp = xyp
    theta = np.arctan2(y - yp, x - xp)
    if clockwise:
        theta *= -1
    return np.argsort(theta)


def add_first_to_end(sequence):
    '''add the first item to the end of a sequence'''
    iter_seq = iter(sequence)
    first = next(iter_seq)
    yield first
    for item in iter_seq:
        yield item
    yield first


def polygongrid_from_teo(teofile, close_cell_polygons=False):
    log.info('creating polygon grid from teo file')
    with adopy.open_grid(teofile) as src:
        teo = src.read()

    # get vertices
    center_coords = teo.get_center_coords()

    # get nodes and vertices per node
    nodes = []
    vertices = []
    iv = 0 
    for nodenumber, node_coords in enumerate(teo.get_node_coords()):
        log.debug('creating node number {nodenumber:d} from TEO'.format(
            nodenumber=nodenumber,
            ))
        node_elems = teo.get_elements_for_node(nodenumber)
        node_vertex_coords = center_coords[node_elems]

        boundary_vertex_coords = []

        # add this node if node is on boundary
        if teo.is_boundary_node(nodenumber):
            boundary_vertex_coords.append(node_coords)

        # get nodes of surrounding elements
        for node_elem in node_elems:
            for neighbornumber in teo.get_nodes_for_element(node_elem):

                # skip node itself
                if neighbornumber == nodenumber:
                    continue

                # add midpoint between node and neighbor
                neighbor_coords = teo.get_node_coords(neighbornumber)
                midpoint_coords = np.mean(
                    (node_coords, neighbor_coords),
                    axis=0)
                boundary_vertex_coords.append(midpoint_coords)

        # stack and concatenate
        boundary_vertex_coords = np.stack(boundary_vertex_coords)
        node_vertex_coords = np.concatenate(
            (node_vertex_coords, boundary_vertex_coords),
            axis=0,
            )
   
        # sort vertices clockwise
        if teo.is_boundary_node(nodenumber):
            midpoint_coords = node_vertex_coords.mean(axis=0)
        else:
            midpoint_coords = node_coords        
        node_vertex_coords = node_vertex_coords[
            sort_around_point(
                node_vertex_coords,
                midpoint_coords,
                clockwise=True,
                ),
            ]

        if close_cell_polygons:
            node_vertex_coords = add_first_to_end(node_vertex_coords)

        # vertices
        node_vertices = [
            (i + iv, xv, yv) for i, (xv, yv) in
            enumerate(node_vertex_coords)
            ]
        iv += len(node_vertices) 

        # append nodes
        node_x, node_y = node_coords
        node_vertex_numbers = [iv for iv, *v in node_vertices]
        nodes.append(
            (
                nodenumber,
                node_x,
                node_y,
                node_vertex_numbers,
                )
            )

        # extend vertices
        vertices.extend(node_vertices)
 
    return PolygonGrid(nodes, vertices, teo.boundary_nodes)