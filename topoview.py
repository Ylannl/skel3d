import math, sys
from time import time
import numpy as np
from pointio import io_npy
from ma_util import MAHelper
from povi import App
from graph import *
from region_growing import *
from geometry import *
from povi import App, Layer, LinkedLayer

from PyQt5 import uic
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QWidget, QToolBox, QTreeWidgetItem

# todo:
# - clipping planes

class MatApp(App):

    def __init__(self, ma, args=[]):
        super(MatApp, self).__init__(args)
        self.ma = ma
        self.segment_filter = ma.D['ma_segment']>0
        self.radius_value = 150.

    def run(self):
        self.layers.add_layer(Layer(name='Clusters'))
        self.dialog = ToolsDialog(self)
        self.draw_clusters()
        self.update_radius(150.)
        super(MatApp, self).run()

    def filter_linkcount(self, value):
        f = ma.D['seg_link_adj'][:,2] >= value
        self.layers['Other']['Adjacency relations'].updateAttributes(filter=np.repeat(f,2))
        # f = ma.D['seg_link_flip'][:,2] >= value
        # self.viewerWindow.data_programs['Flip relations'].updateAttributes(filter=np.repeat(f,2))
        self.viewerWindow.render()

    def filter_component_all(self, toggle):
        for gp in self.layers['Clusters']:
            gp.is_visible = True
        self.layers['Surface'].mask()
        # self.update_radius(self.radius_value)
        if toggle==True:
            self.filter_component(index=self.dialog.ui.comboBox_component.currentIndex())

    def filter_component(self, index):
        for program in self.layers['Clusters']:
            program.is_visible=False
        name = self.dialog.ui.comboBox_component.itemText(index)
        self.layers['Clusters'][name].is_visible = True
        
        g = self.layers['Clusters'][name].graph
        self.viewerWindow.center_view(np.mean(g.vs['ma_coords_mean'], axis=0))
        self.viewerWindow.render()

        ma_idx = np.concatenate(g.vs['ma_idx'])
        if ma_idx.sum() <1:return
        # update mat points
        f = np.zeros(self.ma.m*2, dtype=bool)
        f[ma_idx] = True
        self.layers['MAT'].mask(f)
        # update coords
        # find indices of all surface points related to these mat points

        f = np.zeros(self.ma.m, dtype=bool)
        f[np.mod(ma_idx, self.ma.m)] = True
        f[self.ma.D['ma_qidx'][ma_idx]] = True
        self.layers['Surface'].mask(f)
        
        self.viewerWindow.render()

    def update_radius(self, value):
        self.radius_value = value
        self.radius_filter = self.ma.D['ma_radii'] <= self.radius_value 
        f=np.logical_and(self.segment_filter, self.radius_filter)
        self.layers['MAT'].mask(f)
        return

    def draw_clusters(self):
        self.layers['Clusters'].clear()
        self.dialog.ui.comboBox_component.clear()
        self.dialog.ui.groupBox_component.setChecked(False)
        self.filter_component_all(False)

        min_count = self.dialog.ui.spinBox_linkcount.value()
        contract_thres = self.dialog.ui.doubleSpinBox_contractthres.value()
        g = self.ma.D['ma_segment_graph'].copy()
        # g = g.subgraph(g.vs.select(ma_theta_mean_lt=math.radians(100), up_angle_gt=math.radians(40)))
        g = g.subgraph_edges(g.es.select(adj_count_gt=min_count))
        contract_edges(g, contract_thres)
        
        # self.graphs = []
        # graphlib = get_graph_library()
        # for mapping in g.get_subisomorphisms_vf2(graphlib['flatcube_top']):
        #     self.graphs.append(g.subgraph(mapping))
        
        self.graphs = g.clusters().subgraphs()

        i=0
        for g in self.graphs:
            adj_rel_start = []
            adj_rel_end = []

            if 0<g.ecount():#<1000:
                for e in g.es:
                    adj_rel_start.append(g.vs[e.source]['ma_coords_mean'])
                    adj_rel_end.append(g.vs[e.target]['ma_coords_mean'])
                # import ipdb; ipdb.set_trace()
                # color = np.random.rand(3)
                # color[np.random.random_integers(0,2)] = np.random.uniform(0.5,1.0,1)
                color = np.random.uniform(0.3,1.0,3)
                p = self.layers['Clusters'].add_data_source_line(
                    name = 'cluster {}'.format(i),
                    coords_start = np.array(adj_rel_start),
                    coords_end = np.array(adj_rel_end),
                    color = tuple(color),
                    is_visible=True
                )
                i+=1
                p.graph = g
        self.viewerWindow.render()

        self.layers['Clusters'].is_visible=True

        # populate comboBox_component
        self.dialog.ui.comboBox_component.insertItems(0, [name for name in self.layers['Clusters'].programs.keys()])


class ToolsDialog(QWidget):
    def __init__(self, app, parent=None):
        super(ToolsDialog, self).__init__(parent)
        self.ui = uic.loadUi('tools.ui', self)
        self.app = app

        for layer in self.app.layers:
            item = QTreeWidgetItem([layer.name], 0)
            for program in layer:
                item.addChild(QTreeWidgetItem([program.name], 0))
            self.ui.treeWidget_layers.addTopLevelItem(item)
            self.ui.treeWidget_layers.expandItem(item)
            item.setSelected(True)

        # self.ui.doubleSpinBox_filterRadius.valueChanged.connect(self.app.update_radius)
        self.ui.spinBox_linkcount.valueChanged.connect(self.app.filter_linkcount)
        # self.ui.doubleSpinBox_contractthres.valueChanged.connect(self.app.doubleSpinBox_contractthres)
        self.ui.pushButton_regraph.clicked.connect(self.app.draw_clusters)
        self.ui.groupBox_component.clicked.connect(self.app.filter_component_all)
        self.ui.comboBox_component.activated.connect(self.app.filter_component)
        # self.ui.listWidget_layers.itemSelectionChanged.connect(self.app.set_layer_selection)
        self.ui.treeWidget_layers.itemSelectionChanged.connect(self.app.set_layer_selection)
        # import ipdb; ipdb.set_trace()

    def slot_tcount(self, value):
        print('tcount', value)

# INFILE = 'data/scan_npy'
INFILE = "/Users/ravi/git/masbcpp/rdam_blokken_npy"
# INFILE = "/Volumes/Data/Data/pointcloud/AHN2_matahn_samples/ringdijk_opmeer_npy"
# INFILE = "/Volumes/Data/Data/pointcloud/AHN2_matahn_samples/denhaag_a12_npy"

def timeit(func):
    t0 = time()
    r = func()
    print(("Executed function %s took %f s" % (func.__name__, time()-t0)))
    return r

def assign_seg_point():
    """find for each coord the set of segments it is linked to"""
    pdict = {}

    for i, segment in enumerate(ma.D['ma_segment']):
        coord_id = i# % ma.m
        fp = coord_id
        fq = ma.D['ma_qidx'][coord_id]

        for idx in [fp,fq]:
            if idx in pdict:
                pdict[idx].add(segment)
            else:
                pdict[idx] = set([segment])
    del pdict[-1]

    datadict['segment_count'] = np.array([ len(s) for s in list(pdict.values()) ], dtype=np.int32)
    io_npy.write_npy(INFILE, datadict, ['segment_count'])

def count_refs():
    """count the number of times each coord is used as feature point for a medial ball"""

    pdict = {}
    for i in np.arange(ma.m*2):
        coord_id = i# % ma.m
        fp = coord_id %ma.m
        fq = ma.D['ma_qidx'][coord_id]

        for idx in [fp,fq]:
            if idx in pdict:
                pdict[idx] += 1
            else:
                pdict[idx] = 1
    del pdict[-1]

    return np.array(list(pdict.values()), dtype=np.int32)

def view(ma):
    # ref_count = timeit(count_refs)
    # min_link_adj = 5
    # max_r=190.
    # ma.segment_centers_dict = segment_centers_dict

    # seg_centers = np.array([v[1] for v in list(segment_centers_dict.values())], dtype=np.float32)
    # seg_cnts = np.array([v[0] for v in list(segment_centers_dict.values())], dtype=np.float32)
    # seg_ids = np.array([k for k in list(segment_centers_dict.keys())], dtype=np.float32)

    flip_rel_start = np.zeros((len(ma.D['seg_link_flip']),3), dtype=np.float32)
    flip_rel_end = np.zeros((len(ma.D['seg_link_flip']),3), dtype=np.float32)
    i=0
    for s,e in ma.D['seg_link_flip'][:,:2]:
        flip_rel_start[i] = ma.g.vs[s]['ma_coords_mean']
        flip_rel_end[i] = ma.g.vs[e]['ma_coords_mean']
        i+=1

    adj_rel_start = np.zeros((len(ma.D['seg_link_adj']),3), dtype=np.float32)
    adj_rel_end = np.zeros((len(ma.D['seg_link_adj']),3), dtype=np.float32)
    i=0
    # f = ma.D['seg_link_adj'][:,2] > min_link_adj
    for s,e in ma.D['seg_link_adj'][:,:2]:
        adj_rel_start[i] = ma.g.vs[s]['ma_coords_mean']
        adj_rel_end[i] = ma.g.vs[e]['ma_coords_mean']
        i+=1
    
    c = MatApp(ma)
    layer_s = c.add_layer(LinkedLayer(name='Surface'))
    # layer_ma = c.add_layer(LinkedLayer(name='MAT'))
    layer_ma = c.add_layer(LinkedLayer(name='MAT'))
    layer_misc = c.add_layer(Layer(name='Other'))

    layer_s.add_data_source(
        name = 'Surface points',
        opts=['splat_disk', 'with_normals'],
        points=ma.D['coords'], normals=ma.D['normals'],
    )
    layer_s.add_data_source_line(
      name = 'Surface normals',
      coords_start = ma.D['coords'] + ma.D['normals'],
      coords_end = ma.D['coords'],
      color = (1,1,0)
    )

    if 'ma_segment' in ma.D:
        # f = np.logical_and(ma.D['ma_radii'] < max_r, ma.D['ma_segment']>0)
        layer_ma.add_data_source(
            name = 'MAT points segmented',
            opts=['splat_point', 'with_intensity'],
            points=ma.D['ma_coords'], 
            category=ma.D['ma_segment'].astype(np.float32),
            colormap='random'
        )
    else:
        # f = ma.D['ma_radii_in'] < max_r
        layer_ma.add_data_source(
            name = 'interior MAT',
            opts = ['splat_point', 'blend'],
            points=ma.D['ma_coords_in']
        )
        # f = ma.D['ma_radii_out'] < max_r
        layer_ma.add_data_source(
            name = 'exterior MAT',
            opts = ['splat_point', 'blend'],
            points=ma.D['ma_coords_out']
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
        
    # v_up = np.array([0,0,1],dtype=np.float)
    # biup_angle = np.arccos(np.sum(ma.D['ma_bisec']*v_up, axis=1))
    # f_exterior = np.logical_and(biup_angle < math.pi/2, ma.D['ma_theta'] < math.radians(175)) 
    # f_interior = np.logical_and(biup_angle > math.pi/2, ma.D['ma_theta'] < math.radians(175))

    # c.add_data_source(
    #     name='int',
    #     opts = ['splat_point', 'blend'],
    #     points=ma.D['ma_coords'][f_interior]
    # )
    # c.add_data_source(
    #     name='ext',
    #     opts = ['splat_point', 'blend'],
    #     points=ma.D['ma_coords'][f_exterior]
    # ) 
    # import ipdb; ipdb.set_trace()

    # f = seg_cnts!=1
    # c.add_data_source(
    #     name = 'Segment centers',
    #     opts = ['splat_point'],
    #     points = seg_centers[f]
    # )
    # max_r=150
    # f = np.logical_and(ma.D['ma_radii'] < max_r, ma.D['ma_segment']==0)
    f = ma.D['ma_segment']==0
    layer_misc.add_data_source(
        name = 'MAT points unsegmented',
        opts = ['splat_point', 'blend'],
        points=ma.D['ma_coords'][f]
    )

    if len(flip_rel_start)>0:
        layer_misc.add_data_source_line(
            name = 'Flip relations',
            coords_start = flip_rel_start,
            coords_end = flip_rel_end
        )

    if len(adj_rel_start)>0:
        # f = seg_cnts!=1
        layer_misc.add_data_source_line(
            name = 'Adjacency relations',
            coords_start = adj_rel_start,
            coords_end = adj_rel_end,
            color = (0,1,0)
        )

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
