from skel3d.io import npy 
from skel3d.util import MAHelper
from skel3d.geom3d import Point, Line

from time import time
import sys
import math
import numpy as np
from povi import App

import click

def cluster_spokes(ma,ma_idx):
    cross = np.cross(ma.D['ma_f1'][ma_idx], ma.D['ma_f2'][ma_idx])
    l = Line([Point(v) for v in cross])
    # l.t is a unit vector in one direction of the line
    x = np.empty(len(cross))
    for i, v in enumerate(cross):
        x[i] = np.dot(l.t,v)

    return x > 0, cross


@click.command(help='Classify and visualise spokes for a MAT sheet')
@click.argument('input_npy', type=click.Path(exists=True))
@click.argument('cluster_id', type=int)
@click.argument('sheet_id', type=int)
@click.option('-n', '--near_clip', default=0.1, type=float, help='Near clip value.')
@click.option('-f', '--far_clip', default=100.0, type=float, help='Far clip value.')
def mat(input_npy, cluster_id, sheet_id, near_clip, far_clip):
    c = App(near_clip=near_clip, far_clip=far_clip)
    
    t0=time()
    datadict = npy.read(input_npy)
    ma = MAHelper(datadict)
    min_count = 5
    # print(("{} points loaded from file in {} s".format(ma.m, time()-t0)))

    # import ipdb; ipdb.set_trace()
    master_g = ma.D['ma_segment_graph']
    master_g = master_g.subgraph_edges(master_g.es.select(adj_count_gt=min_count))
    # contract_edges(g, contract_thres)

    ma.D['ma_segment'] = np.zeros(ma.m*2,dtype=np.int64)
    for v in master_g.vs:
        ma.D['ma_segment'][v['ma_idx']] = v.index
        v['s_id'] = v.index
    
    graphs = master_g.clusters().subgraphs()
    sheet = graphs[cluster_id].vs[sheet_id]
    ma_idx = np.array(sheet['ma_idx'])

    clustermap, cross = cluster_spokes(ma, ma_idx)
    print(clustermap)

    side_A = ma.D['ma_f1'][ma_idx]
    side_A[~clustermap] =  ma.D['ma_f2'][ma_idx][~clustermap]

    side_B = ma.D['ma_f1'][ma_idx]
    side_B[clustermap] =  ma.D['ma_f2'][ma_idx][clustermap]


    ma_coords = ma.D['ma_coords'][ma_idx]

    # import ipdb; ipdb.set_trace()

    c.add_data_source_line(
        name = 'cluster A',
        coords_start = ma_coords,
        coords_end = side_A+ma_coords,
        color=(.9,0,0)
    )

    c.add_data_source_line(
        name = 'cluster B',
        coords_start = ma_coords,
        coords_end = side_B+ma_coords,
        color=(0,0,0.9)
    )

    c.add_data_source_line(
        name = 'cross',
        coords_start = ma_coords[clustermap],
        coords_end = cross[clustermap]+ma_coords[clustermap],
        color=(0.9,0.9,0.9),
        options = ['alternate_vcolor']
    )


    c.run()

if __name__ == '__main__':
    mat()