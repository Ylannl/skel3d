import math, sys
from time import time
import numpy as np

from itertools import chain
from PyQt5.QtCore import Qt
from pyqtgraph import PlotWidget

from povi import App, Layer, LinkedLayer
from skel3d.io import npy
from skel3d.util import MAHelper
from skel3d.graph import *
from skel3d.segmentation import *
from skel3d.polyhedralise import *

class TestApp(App):

    def __init__(self, ma, args=[]):
        super(TestApp, self).__init__(args)
        self.ma = ma
        
        # self.plotWindow.plot(np.random.normal(size=100), name="Data 1")
    
    def addGraphWindow(self, g): 
        self.plotWindow = GraphWindow(g, master_app=self)
        # self.plotWindow.addLegend()
        self.plotWindow.show()
        # super(TestApp, self).run()
        
class GraphWindow(PlotWidget):
    def __init__(self, g, master_app, parent=None):
        self.layer_manager = master_app.layer_manager
        self.master_app = master_app
        self.avtive_sheet = 0
        self.g=g
        self.iterator = self.sheet_loop()
        super(GraphWindow, self).__init__(parent)
    
    def sheet_loop(self):
        while True:
            for v in self.g.vs:
                yield v

    def keyPressEvent(self, event):
        key = event.key()
        repeat = event.isAutoRepeat()
        # print('keypressevent!')
        if key == Qt.Key_R:
            self.layer_manager['Surface'].mask(None)
            self.layer_manager.layers['MAT'].mask(None)
            self.master_app.viewerWindow.render()
        if key == Qt.Key_N:
            self.clear()
            v = self.iterator.next()
            
            ma_idx = v['ma_idx']
            f = np.zeros(ma.m*2, dtype=bool)
            f[ma_idx] = True
            self.layer_manager['MAT'].mask(f)
            
            f = np.zeros(ma.m, dtype=bool)
            f[np.mod(ma_idx, ma.m)] = True
            f[self.master_app.ma.D['ma_qidx'][ma_idx]] = True
            self.layer_manager.layers['Surface'].mask(f)

            self.master_app.viewerWindow.center_view(np.mean(self.master_app.ma.D['ma_coords'][ma_idx], axis=0))
            self.master_app.viewerWindow.render()
            
            # import ipdb/;ipdb.set_trace()
            i = v.index
            vals=ma.D['ma_radii'][ma_idx]
            b=v['ma_bisec_mean']
            v['ma_theta_mean']
            ## compute standard histogram
            y,x = np.histogram(vals, bins=50)

            ## Using stepMode=True causes the plot to draw two lines for each sample.
            ## notice that len(x) == len(y)+1
            color = tuple(np.random.uniform(0.3,1.0,3)*255) + (255,)
            self.plot(x, y, stepMode=True, fillLevel=0, pen={'color': color, 'width': 2}, name='radius '+str(i))



def view(ma, vids):
    # ref_count = timeit(count_refs)
    max_r=190.
    # ma.g = ma.D['ma_segment_graph']
    
    c = TestApp(ma)

    min_count = 30
    contract_thres = 20
    g = ma.D['ma_segment_graph'].copy()
    # g = g.subgraph(g.vs.select(ma_theta_mean_lt=math.radians(100), up_angle_gt=math.radians(40)))
    g = g.subgraph_edges(g.es.select(adj_count_gt=min_count))
    # contract_edges(g, contract_thres)

    
    graphlib = get_graph_library()
    vertex_clustering = g.clusters()
    # import ipdb;ipdb.set_trace()
    for that_id in [8]:
        this_g = vertex_clustering.subgraph(that_id)

        # c.addGraphWindow(this_g)
        
        try:
            
            this_m = build_map(this_g, ma)
            layer_thisg = c.add_layer(Layer(name='cluster '+str(that_id)))

            adj_rel_start = []
            adj_rel_end = []
            # this_g = g.subgraph(this_mapping)
            for hn in this_m.ns:
                hn['coords_mean'] = np.mean(ma.D['coords'][hn['s_idx']], axis=0)
            # import ipdb;ipdb.set_trace()
            for e in this_m.es:
                if e.kind == 'match':
                    source, target = e.nodes
                    adj_rel_start.append(source['coords_mean'])
                    adj_rel_end.append(target['coords_mean'])
            color = np.random.uniform(0.3,1.0,3)
            p = layer_thisg.add_data_source_line(
                name = 'map edges {}'.format(that_id),
                coords_start = np.array(adj_rel_start),
                coords_end = np.array(adj_rel_end),
                color = (0,1,0),
                is_visible=True,
                options = ['alternate_vcolor']
            )

            adj_rel_start = []
            adj_rel_end = []
            for i in range(len(this_m.ns)//2):
                hn = this_m.ns[i*2]
                source, target = hn, hn.twin
                adj_rel_start.append(source['coords_mean'])
                adj_rel_end.append(target['coords_mean'])
            color = np.random.uniform(0.3,1.0,3)
            p = layer_thisg.add_data_source_line(
                name = 'map twin links {}'.format(that_id),
                coords_start = np.array(adj_rel_start),
                coords_end = np.array(adj_rel_end),
                color = (1,1,1),
                is_visible=False
            )
        # except Exception as e:
        #     continue

        # try:
            planes = polyhedral_reconstruct(this_m, ma)
            for i, (coords, normals) in enumerate(planes):
                layer_thisg.add_data_source_triangle(
                    name = 'plane '+str(that_id)+' '+str(i),
                    coords = coords,
                    normals = normals,
                    color = (0.88,1.0,1.0),
                    is_visible = False,
                    # draw_type='line_loop'
                    draw_type='triangles'
                )
        except Exception as e:
            print('polyhedral_reconstruct failed')
            raise
        # for i,pts in enumerate(pointsets):
        #     layer_thisg.add_data_source(
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
        p = layer_thisg.add_data_source_line(
            name = 'matched graph',
            coords_start = np.array(adj_rel_start),
            coords_end = np.array(adj_rel_end),
            color = (0.,0.9,0.9),
            is_visible=False
        )

    layer_ma = c.add_layer(LinkedLayer(name='MAT'))
    layer_s = c.add_layer(LinkedLayer(name='Surface'))

    layer_s.add_data_source(
        name = 'Surface points',
        opts=['splat_disk', 'with_normals', 'fixed_color'],
        points=ma.D['coords'], 
        normals=ma.D['normals'],
        color = (.4,.4,1.)
    )
    layer_s.add_data_source_line(
      name = 'Surface normals',
      coords_start = ma.D['coords'] + ma.D['normals'],
      coords_end = ma.D['coords'],
      color = (1,1,0)
    )
    f_s = np.zeros(ma.m, dtype=bool)
    f_s[s_idx] = True
    f_s[ma.D['ma_qidx'][ma_ids]] = True
    layer_s.mask(f_s)
    

    for v in g.vs():
        ma.D['ma_segment'][ v['ma_idx'] ] = v.index
    # f =ma.D['ma_segment'] != 0
    
    # compute mat 'normals'
    # cross product of spokes is perpendicular to bisector and tangent to sheet
    vec_coplanar = np.cross(ma.D['ma_f1'],ma.D['ma_f2'])
    # now compute this cross product to find a vector in the normal direction of the plane that we want to reconstruct
    ma_n = np.cross(vec_coplanar, ma.D['ma_bisec'])
    ma_n = ma_n / np.linalg.norm(ma_n, axis=1)[:,None]

    layer_ma.add_data_source(
        name = 'MAT points',
        opts=['splat_disk', 'with_normals', 'with_intensity'],
        points=ma.D['ma_coords'], 
        normals=ma_n,
        category=ma.D['ma_segment'].astype(np.float),
        colormap='random',
        default_mask=ma.D['ma_segment'] != 0
    )             
    layer_ma.add_data_source_line(
      name = 'Primary spokes',
      coords_start = ma.D['ma_coords'],
      coords_end = np.concatenate([ma.D['coords'],ma.D['coords']])
    )
    layer_ma.add_data_source_line(
      name = 'Secondary spokes',
      coords_start = ma.D['ma_coords'],
      coords_end = np.concatenate([ma.D['coords'],ma.D['coords']])[ma.D['ma_qidx']]
    )
    layer_ma.add_data_source_line(
        name = 'Bisectors',
        coords_start = ma.D['ma_coords'],
        coords_end = ma.D['ma_bisec']+ma.D['ma_coords'],
        color=(.2,.2,1)
    )
    f = np.zeros(2*ma.m, dtype=bool)
    f[ma_ids] = True
    layer_ma.mask(f)

    c.viewerWindow.center_view(center=np.mean(ma.D['coords'][f_s], axis=0))
    c.run()

if __name__ == '__main__':
    vids = [4]
    #box: 4,6,8 (3,10 broken in part)
    #sloped box: 7
    #simple gable: 9
    #exterior: 5
    if len(sys.argv)>1:
        vids = [int(sys.argv[-1])]
        # INFILE = sys.argv[-1]
    # import ipdb;ipdb.set_trace()
    INFILE = "/Users/ravi/git/mat_util/Random3Dcity/1/NPY"
    # INFILE = "/Users/ravi/git/mat_util/test_cases/sloped_gable/NPY"
    datadict = npy.read(INFILE)
    ma = MAHelper(datadict, origin=True)

    g = ma.D['ma_segment_graph']
    for v in g.vs:
        v['up_angle'] = np.sign(v['ma_bisec_mean'])[2] * np.arccos(np.dot(v['ma_bisec_mean'], [0,0,1] ))

    view(ma, vids)