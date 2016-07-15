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
from itertools import chain

from pyqtgraph import PlotWidget

class TestApp(App):

    def __init__(self, ma, args=[]):
        super(TestApp, self).__init__(args)
        self.ma = ma
        self.plotWindow = PlotWidget()
        self.plotWindow.addLegend()
        # self.plotWindow.plot(np.random.normal(size=100), name="Data 1")
    
    def run(self): 
        self.plotWindow.show()
        super(TestApp, self).run()
        

def view(ma, vids):
    # ref_count = timeit(count_refs)
    max_r=190.
    # ma.g = ma.D['ma_segment_graph']
    
    c = TestApp(ma)

    min_count = 15
    contract_thres = 15
    g = ma.D['ma_segment_graph'].copy()
    # g = g.subgraph(g.vs.select(ma_theta_mean_lt=math.radians(100), up_angle_gt=math.radians(40)))
    g = g.subgraph_edges(g.es.select(adj_count_gt=min_count))
    # contract_edges(g, contract_thres)

    
    graphlib = get_graph_library()
    vertex_clustering = g.clusters()
    # import ipdb;ipdb.set_trace()
    for that_id in vids:
        this_g = vertex_clustering.subgraph(that_id)

        for i, ma_idx in enumerate(this_g.vs['ma_idx']):
            vals=ma.D['ma_radii'][ma_idx]
            ## compute standard histogram
            y,x = np.histogram(vals, bins=50)

            ## Using stepMode=True causes the plot to draw two lines for each sample.
            ## notice that len(x) == len(y)+1
            color = tuple(np.random.uniform(0.3,1.0,3)*255) + (255,)
            c.plotWindow.plot(x, y, stepMode=True, fillLevel=0, pen={'color': color, 'width': 2}, name='radius '+str(i))

        this_m = build_map(this_g, ma)

        adj_rel_start = []
        adj_rel_end = []
        # this_g = g.subgraph(this_mapping)
        for hn in this_m.ns:
            hn['coords_mean'] = np.mean(ma.D['coords'][hn['s_idx']], axis=0)
        for e in this_m.es:
            if e.kind == 'match':
                source, target = e.nodes
                adj_rel_start.append(source['coords_mean'])
                adj_rel_end.append(target['coords_mean'])
        color = np.random.uniform(0.3,1.0,3)
        p = c.add_data_source_line(
            name = 'map edges {}'.format(that_id),
            coords_start = np.array(adj_rel_start),
            coords_end = np.array(adj_rel_end),
            color = (0,1,0),
            is_visible=False
        )

        adj_rel_start = []
        adj_rel_end = []
        for i in range(len(this_m.ns)/2):
            hn = this_m.ns[i*2]
            source, target = hn, hn.twin
            adj_rel_start.append(source['coords_mean'])
            adj_rel_end.append(target['coords_mean'])
        color = np.random.uniform(0.3,1.0,3)
        p = c.add_data_source_line(
            name = 'map twin links {}'.format(that_id),
            coords_start = np.array(adj_rel_start),
            coords_end = np.array(adj_rel_end),
            color = (1,1,1),
            is_visible=False
        )

        try:
            planes = polyhedral_reconstruct(this_m, ma)
            for i, (coords, normals) in enumerate(planes):
                c.add_data_source_triangle(
                    name = 'roof surface _'+str(that_id)+' '+str(i),
                    coords = coords[0:6],
                    normals = normals[0:6],
                    color = (1.0,0.1,0.1),
                    is_visible = True
                )
        except Exception:
            print('polyhedral_reconstruct failed')
        # for i,pts in enumerate(pointsets):
        #     c.add_data_source(
        #         name = 'Surface points _'+' vid '+ ' - ' +str(i),
        #         opts=['splat_disk', 'with_normals'],
        #         points=ma.D['coords'][list(pts)], 
        #         normals=ma.D['normals'][list(pts)],
        #     )
        
    ma_ids = []
    for ma_idx in this_g.vs['ma_idx']:
        ma_ids += ma_idx
    # ma_idx = ma_idx[0] + ma_idx[1] +ma_idx[2] +ma_idx[3]
    s_idx = np.mod(ma_ids, ma.m)
    
    # normals = np.array(normals, dtype=np.float32)

    # normals_0 = normals[labels==0]
    # normals_1 = normals[labels==1]
    # print normals
    
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

    adj_rel_start = []
    adj_rel_end = []
    # this_g = g.subgraph(this_mapping)
    for e in this_g.es:
        adj_rel_start.append(this_g.vs[e.source]['ma_coords_mean'])
        adj_rel_end.append(this_g.vs[e.target]['ma_coords_mean'])
    # import ipdb; ipdb.set_trace()
    # color = np.random.rand(3)
    # color[np.random.random_integers(0,2)] = np.random.uniform(0.5,1.0,1)
    color = np.random.uniform(0.3,1.0,3)
    p = c.add_data_source_line(
        name = 'matched graph',
        coords_start = np.array(adj_rel_start),
        coords_end = np.array(adj_rel_end),
        color = (0.,0.9,0.9),
        is_visible=True
    )

    c.add_data_source(
        name = 'Surface points all',
        opts=['splat_disk', 'with_normals'],
        points=ma.D['coords'], 
        normals=ma.D['normals'],
    )

    c.add_data_source(
        name = 'Surface points',
        opts=['splat_disk', 'with_normals'],
        points=np.concatenate([ma.D['coords'][s_idx], ma.D['coords'][ma.D['ma_qidx'][ma_ids]]]), 
        normals=np.concatenate([ma.D['normals'][s_idx], ma.D['normals'][ma.D['ma_qidx'][ma_ids]]]),
    )
    c.add_data_source_line(
      name = 'Surface normals',
      coords_start = np.concatenate([ma.D['coords'][s_idx], ma.D['coords'][ma.D['ma_qidx'][ma_ids]]]) + np.concatenate([ma.D['normals'][s_idx], ma.D['normals'][ma.D['ma_qidx'][ma_ids]]]),
      coords_end = np.concatenate([ma.D['coords'][s_idx], ma.D['coords'][ma.D['ma_qidx'][ma_ids]]]),
      color = (1,1,0)
    )

    for v in g.vs():
		ma.D['ma_segment'][ v['ma_idx'] ] = v.index	
    c.add_data_source(
        name = 'MAT points',
        opts=['splat_point', 'with_intensity'],
        points=ma.D['ma_coords'][ma_ids], 
        category=ma.D['ma_segment'][ma_ids].astype(np.float32),
        colormap='random'
    )

    f =ma.D['ma_segment'] != 0
    c.add_data_source(
        name = 'MAT points all',
        opts=['splat_point', 'with_intensity'],
        points=ma.D['ma_coords'][f], 
        category=ma.D['ma_segment'][f].astype(np.float32),
        colormap='random'
    )
        
    
    # c.add_data_source_line(
    #   name = 'normals a',
    #   coords_start = np.concatenate([ma.D['coords'],ma.D['coords']])[ma_ids] + np.concatenate([ma.D['normals'],ma.D['normals']])[ma_ids],
    #   coords_end = np.concatenate([ma.D['coords'],ma.D['coords']])[ma_ids],
    #   color=(1,1,0)
    # )
    # c.add_data_source_line(
    #   name = 'Secondary spokes',
    #   coords_start = ma.D['ma_coords'],
    #   coords_end = np.concatenate([ma.D['coords'][ma.D['ma_qidx_in']],ma.D['coords'][ma.D['ma_qidx_out']]])
    # )

    c.add_data_source_line(
      name = 'Primary spokes',
      coords_start = ma.D['ma_coords'][ma_ids],
      coords_end = np.concatenate([ma.D['coords'],ma.D['coords']])[ma_ids]
    )
    c.add_data_source_line(
      name = 'Secondary spokes',
      coords_start = ma.D['ma_coords'][ma_ids],
      coords_end = np.concatenate([ma.D['coords'],ma.D['coords']])[ma.D['ma_qidx']][ma_ids]
    )

    c.add_data_source_line(
        name = 'Bisectors',
        coords_start = ma.D['ma_coords'][ma_ids],
        coords_end = ma.D['ma_bisec'][ma_ids]+ma.D['ma_coords'][ma_ids],
        color=(0,1,0)
    )

    c.run()

if __name__ == '__main__':
    vids = [4]
    #box: 4,6,8 (3,10 broken in part, 7 very broken)
    #simple gable: 9
    #exterior: 5
    if len(sys.argv)>1:
        vids = [int(sys.argv[-1])]
        # INFILE = sys.argv[-1]
    # import ipdb;ipdb.set_trace()
    INFILE = "/Users/ravi/git/mat_util/Random3Dcity/NPY"
    datadict = io_npy.read_npy(INFILE)
    ma = MAHelper(datadict, origin=True)

    g = ma.D['ma_segment_graph']
    for v in g.vs:
        v['up_angle'] = np.sign(v['ma_bisec_mean'])[2] * np.arccos(np.dot(v['ma_bisec_mean'], [0,0,1] ))

    view(ma, vids)