import math, sys
from time import time
import numpy as np
from pointio import io_npy
from ma_util import MAHelper
from povi import App
from graph import *
from region_growing import *
from geometry import *
from povi import App

def view(ma):
    # ref_count = timeit(count_refs)
    min_link_adj = 5
    max_r=190.
    # ma.g = ma.D['ma_segment_graph']
    
    

    min_count = 5
    contract_thres = 15
    g = ma.D['ma_segment_graph'].copy()
    g = g.subgraph(g.vs.select(ma_theta_mean_lt=math.radians(100), up_angle_gt=math.radians(40)))
    g = g.subgraph_edges(g.es.select(adj_count_gt=min_count))
    contract_edges(g, contract_thres)
    
    graphlib = get_graph_library()
    this_mapping = g.get_subisomorphisms_vf2(graphlib['flatcube_top'])[2]

    ma_idx, normals, labels = gf_flatcube_top(g, this_mapping, ma)
    normals = np.array(normals, dtype=np.float32)

    normals_0 = normals[labels==0]
    normals_1 = normals[labels==1]
    # print normals
    c = App()
    # c.add_data_source_line(
    #   name = 'normals 0',
    #   coords_start = np.zeros(normals_0.shape, dtype=np.float32),
    #   coords_end = normals_0,
    #   color = (0,1,0)
    # )

    # c.add_data_source_line(
    #   name = 'normals 1',
    #   coords_start = np.zeros(normals_1.shape, dtype=np.float32),
    #   coords_end = normals_1,
    #   color = (1,0,0)
    # )
    # c.viewerWindow.data_center = (0,0,0)
    # print c.viewerWindow.data_center, c.viewerWindow.data_width

    c.add_data_source(
        name = 'Surface points',
        opts=['splat_disk', 'with_normals'],
        points=ma.D['coords'][ma_idx], normals=ma.D['normals'][ma_idx],
    )

    c.add_data_source(
        name = 'MAT points',
        opts=['splat_point', 'with_intensity'],
        points=ma.D['ma_coords'][ma_idx], 
        category=ma.D['ma_segment'][ma_idx].astype(np.float32),
        colormap='random'
    )
        
    
    c.add_data_source_line(
      name = 'normals a',
      coords_start = np.concatenate([ma.D['coords'],ma.D['coords']])[ma_idx] + np.concatenate([ma.D['normals'],ma.D['normals']])[ma_idx],
      coords_end = np.concatenate([ma.D['coords'],ma.D['coords']])[ma_idx]
    )
    # c.add_data_source_line(
    #   name = 'Secondary spokes',
    #   coords_start = ma.D['ma_coords'],
    #   coords_end = np.concatenate([ma.D['coords'][ma.D['ma_qidx_in']],ma.D['coords'][ma.D['ma_qidx_out']]])
    # )

    # c.add_data_source_line(
    #   name = 'Primary spokes',
    #   coords_start = ma.D['ma_coords'],
    #   coords_end = np.concatenate([ma.D['coords'],ma.D['coords']])
    # )
    # c.add_data_source_line(
    #   name = 'Secondary spokes',
    #   coords_start = ma.D['ma_coords'],
    #   coords_end = np.concatenate([ma.D['coords'][ma.D['ma_qidx_in']],ma.D['coords'][ma.D['ma_qidx_out']]])
    # )

    c.run()

if __name__ == '__main__':
    if len(sys.argv)>1:
        INFILE = sys.argv[-1]
    # import ipdb;ipdb.set_trace()
    datadict = io_npy.read_npy(INFILE)
    ma = MAHelper(datadict, origin=True)

    g = ma.D['ma_segment_graph']
    for v in g.vs:
        v['up_angle'] = np.sign(v['ma_bisec_mean'])[2] * np.arccos(np.dot(v['ma_bisec_mean'], [0,0,1] ))

    view(ma)