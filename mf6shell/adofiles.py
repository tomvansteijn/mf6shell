#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Tom van Steijn, Royal HaskoningDHV

from mf6shell.grid import PolygonGrid

import numpy as np
import adopy

import logging
import re
import os

log = logging.getLogger(os.path.basename(__file__))


def get_layer_from_name(name):
    if name.startswith('RP'):
        return name, None
    m = re.search(r'(?P<name>[A-Z]+)(?P<layer>\d+)', name)
    if m is None:
        return name, None
    else:
        name = m.group('name')
        layer = int(m.group('layer'))
        return name, layer


def get_modflow_name(name):
    return {
        'RL': 'top',
        'TH': 'bot',
        'TX': 'kd',
        'CL': 'c',
        'HH': 'start',
        'BH': 'chd',
        'RP1': 'recharge',
        'SQ': 'wel',
        }.get(name, name)


def read_ado(adofile, block_name, masked=True, nodata=-999.):
    log.debug('reading {f.name:}'.format(f=adofile))
    with adopy.open(adofile) as src:
        for block in src.read():
            if block.name == block_name:
                return np.ma.masked_equal(block.values, nodata)


def sort_around_point(xy, xyp=(0., 0.), clockwise=False):
    '''sort points around center point'''
    x, y = xy.transpose()
    xp, yp = xyp
    theta = np.arctan2(y - yp, x - xp)
    if clockwise:
        theta *= -1
    return np.argsort(theta)


def pseudoangle(dy, dx):
    '''arctan2 is faster'''
    return np.sign(dy) * (1. - dx / (np.fabs(dx) + np.fabs(dy)))


def pseudo_sort_around_point(xy, xyp=(0., 0.), clockwise=False):
    x, y = xy.transpose()
    xp, yp = xyp
    theta = pseudoangle(y - yp, x - xp)
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


def voronoi_from_teo(teofile, close_cell_polygons=False):
    log.info('creating voronoi grid from teo file')
    with adopy.open_grid(teofile) as src:
        teo = src.read()

    # midpoints per element
    x1 = teo.x_nodes[teo.elem1]
    x2 = teo.x_nodes[teo.elem2]
    x3 = teo.x_nodes[teo.elem3]
    y1 = teo.y_nodes[teo.elem1]
    y2 = teo.y_nodes[teo.elem2]
    y3 = teo.y_nodes[teo.elem3]
    (mx1, my1), (mx2, my2), (mx3, my3) = (
            (np.mean([x1, x2], axis=0), np.mean([y1, y2], axis=0)),
            (np.mean([x2, x3], axis=0), np.mean([y2, y3], axis=0)),
            (np.mean([x3, x1], axis=0), np.mean([y3, y1], axis=0)),
            )

    with np.errstate(divide='ignore', invalid='ignore'):
        # opposite  reciprocal slopes per element 1, 2
        s1, s2 = (x1 - x2) / (y2 - y1), (x2 - x3) / (y3 - y2)

        # intercepts per element 1, 2
        b1, b2 = my1 - s1*mx1, my2 - s2*mx2

        # intersection points per element
        dn = s1 - s2
        xc, yc = (b2 - b1) / dn, (s1*b2 - s2*b1) / dn

        # update where line is flat
        xc = np.where(np.isinf(s1),
            mx1,
            xc,
            )
        yc = np.where(np.isinf(s1),
            s2*mx1 + b2,
            yc,
            )

        xc = np.where(np.isinf(s2),
            mx2,
            xc,
            )
        yc = np.where(np.isinf(s2),
            s1*mx2 + b1,
            yc,
            )

    center_coords = np.stack((xc, yc), axis=1)

    # get nodes and vertices per node
    nodes = []
    vertices = []
    iv = 0 
    for nodenumber, node_coords in enumerate(teo.get_node_coords()):

        node_elems = teo.get_elements_for_node(nodenumber)
        node_vertex_coords = center_coords[node_elems]
       
        # get nodes of surrounding elements
        if teo.is_boundary_node(nodenumber):
            boundary_vertex_coords = []
            boundary_vertex_coords.append(node_coords)

            for node_elem in node_elems:
                elem_nodes = teo.get_nodes_for_element(node_elem)
                n1 = next(elem_nodes)
                n2 = next(elem_nodes)
                n3 = next(elem_nodes)
                if teo.is_boundary_node(n1) and teo.is_boundary_node(n2):
                    boundary_vertex_coords.append(
                        (mx1[node_elem], my1[node_elem])
                        )
                if teo.is_boundary_node(n2) and teo.is_boundary_node(n3):
                    boundary_vertex_coords.append(
                        (mx2[node_elem], my2[node_elem])
                        )
                if teo.is_boundary_node(n3) and teo.is_boundary_node(n1):
                    boundary_vertex_coords.append(
                        (mx3[node_elem], my3[node_elem])
                        )

            # stack and concatenate
            boundary_vertex_coords = np.stack(boundary_vertex_coords)
            node_vertex_coords = np.concatenate(
                (node_vertex_coords, boundary_vertex_coords),
                axis=0,
                )
   
        # sort vertices clockwise
        if teo.is_boundary_node(nodenumber):
            midpoint = np.mean(node_vertex_coords, axis=0)  
        else:
            midpoint = node_coords   
        node_vertex_coords = node_vertex_coords[
            sort_around_point(
                node_vertex_coords,
                midpoint,
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

    return PolygonGrid(nodes,vertices,
        teo.nia,
        teo.boundary_nodes,
        teo.river_nodes,
        teo.source_nodes,
        )


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
 
    return PolygonGrid(nodes,vertices,
        teo.nia,
        teo.boundary_nodes,
        teo.river_nodes,
        teo.source_nodes,
        )