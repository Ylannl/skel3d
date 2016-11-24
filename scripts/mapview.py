from skel3d.io import npy 
from skel3d.util import MAHelper
from skel3d.polyhedralise import *

from time import time
import sys
import math
import numpy as np
from povi import App, Layer, LinkedLayer

# import click

from PyQt5.QtWidgets import QApplication

app = QApplication.instance() # retrieves the ipython qt application if any
if app is None:
    app = QApplication([]) # create one if standalone execution
povi_app=None

# class Interactor:
#     def __init__(self, app):
#         self.app = app

#     def set_cluster(self):


#     def draw_updateMAP(self):

#     def draw_polyhedron(self):


def get_MAPData(map, ma):
    adj_rel_start = []
    adj_rel_end = []

    for hn in this_m.ns:
        hn['coords_mean'] = np.mean(ma.D['coords'][hn['s_idx']], axis=0)
    # import ipdb;ipdb.set_trace()
    for e in this_m.es:
        if e.kind == 'match':
            source, target = e.nodes
            adj_rel_start.append(source['coords_mean'])
            adj_rel_end.append(target['coords_mean'])

    adj_rel_start_node = []
    adj_rel_end_node = []
    for i in range(len(this_m.ns)//2):
        hn = this_m.ns[i*2]
        source, target = hn, hn.twin
        adj_rel_start_node.append(source['coords_mean'])
        adj_rel_end_node.append(target['coords_mean'])
    return {'edges':(adj_rel_start, adj_rel_end), 'nodes':(adj_rel_start_node, adj_rel_end_node)}

def draw_MAP(povi_app, map, ma):

    this_m = build_map(that_g, ma)
    map_data = get_MAPData(this_m, ma)

    # self.app.layer_manager['Clusters'][name]
    if not 'Cluster' in povi_app.layer_manager.layers:
        layer_thisg = povi_app.add_layer(Layer(name='Cluser'))
        color = np.random.uniform(0.3,1.0,3)
        p = layer_thisg.add_data_source_line(
            name = 'map edges',
            coords_start = np.array(map_data['edges'][0]),
            coords_end = np.array(map_data['edges'][1]),
            color = (0,1,0),
            is_visible=True,
            options = ['alternate_vcolor']
        )

        p = layer_thisg.add_data_source_line(
            name = 'map twin links',
            coords_start = np.array(map_data['nodes'][0]),
            coords_end = np.array(map_data['nodes'][1]),
            color = (1,1,1),
            is_visible=False
        )
    else:
        layer_thisg = self.app.layer_manager['Clusters']
        p=layer_thisg['map edges']
        

        vertices = np.empty((m*2,n), dtype=np.float32)
        vertices[0::2] = map_data['edges'][0]
        vertices[1::2] = map_data['edges'][1]]

        p.updateAttribute


def mapview(input_npy='../Random3Dcity/1/NPY', near_clip=0.1, far_clip=100.):

    povi_app = App(qapplication_instance=app, near_clip=near_clip, far_clip=far_clip)
    
    datadict = npy.read(input_npy)
    ma = MAHelper(datadict)
    that_id = 12
    that_g = ma.D['ma_clusters'][that_id]

    layer_s = povi_app.add_layer(LinkedLayer(name='Surface'))
    layer_s.add_data_source(
        name = 'Surface points',
        opts=['splat_disk', 'with_normals'],
        points=datadict['coords'], normals=datadict['normals']
    )
    layer_s.add_data_source_line(
      name = 'point normals',
      coords_start = datadict['coords'] + datadict['normals'],
      coords_end = datadict['coords'],
      color = (1,1,0)
    )


    # import ipdb; ipdb.set_trace()



    povi_app.run()
    # import ipdb; ipdb.set_trace()
    return povi_app


# if __name__ == '__main__':
povi_app=mapview()