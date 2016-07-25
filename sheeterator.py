import math, sys
from time import time
import numpy as np
from pointio import io_npy
from ma_util import MAHelper
from povi import App, Layer, LinkedLayer, ToolsDialog
from graph import *
from region_growing import *
from geometry import *
from povi import App
from itertools import chain
from PyQt5.QtCore import Qt
from pyqtgraph import PlotWidget

class TestApp(App):

    def __init__(self, ma, args=[]):
        super(TestApp, self).__init__(args)
        self.ma = ma

    def run(self):
        self.addGraphWindow()
        self.layer_manager.add_layer(Layer(name='Clusters', is_aggregate=True))
        # here all layers need to have been added
        self.dialog = ToolsWindow(self)
        self.draw_clusters()
        super(TestApp, self).run()
        
        # self.plotWindow.plot(np.random.normal(size=100), name="Data 1")

    def draw_clusters(self):

        min_count = 15 #self.dialog.ui.spinBox_linkcount.value()
        contract_thres = 15 #self.dialog.ui.doubleSpinBox_contractthres.value()
        # g = g.subgraph(g.vs.select(ma_theta_mean_lt=math.radians(100), up_angle_gt=math.radians(40)))
        g = self.ma.D['ma_segment_graph']
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
                color = np.random.uniform(0.3,1.0,3)
                p = self.layer_manager['Clusters'].add_data_source_line(
                    name = 'cluster {}'.format(i),
                    coords_start = np.array(adj_rel_start),
                    coords_end = np.array(adj_rel_end),
                    color = tuple(color),
                    is_visible=True
                )
                i+=1
                p.graph = g
        self.viewerWindow.render()

        self.layer_manager['Clusters'].is_visible=True
        self.active_graph = None

        # populate comboBox_component
        self.dialog.ui.comboBox_clusters.insertItems(0, ['All']+[name for name in self.layer_manager['Clusters'].programs.keys()])

    def toggle_selection(self, toggle):
        if toggle==True:
            self.filter_cluster(index=self.dialog.ui.comboBox_clusters.currentIndex())
        else:
            self.filter_idx()
        # self.viewerWindow.render()

    def filter_cluster(self, index):
        if index == 0:
            self.filter_idx(None)
            self.active_graph = None
            self.dialog.ui.comboBox_sheets.clear()
            return

        for program in self.layer_manager['Clusters']:
            program.is_visible=False
        name = self.dialog.ui.comboBox_clusters.itemText(index)
        self.layer_manager['Clusters'][name].is_visible = True
        
        g = self.layer_manager['Clusters'][name].graph
        self.active_graph = g
        if self.dialog.ui.checkBox_clusterCenterView.isChecked():
            self.viewerWindow.center_view(np.mean(g.vs['ma_coords_mean'], axis=0))
            self.viewerWindow.render()

        ma_idx = np.concatenate(g.vs['ma_idx'])
        # if ma_idx.sum() <1:return
        self.plotWindow.show_sheet(ma_idx)
        self.filter_idx(ma_idx)

        # populate sheet list
        self.dialog.ui.comboBox_sheets.clear()
        self.dialog.ui.comboBox_sheets.insertItems(0, ['All']+['Sheet '+str(v.index) for v in g.vs])

    def filter_sheet(self, index):
        if index == 0:
            self.filter_cluster(self.dialog.ui.comboBox_clusters.currentIndex())
        else:
            ma_idx = self.active_graph.vs[index-1]['ma_idx']
            # if ma_idx.sum() <1:return
            self.plotWindow.show_sheet(ma_idx)
            self.filter_idx(ma_idx)
        

    def filter_idx(self, ma_idx=None):
        # update mat points
        if ma_idx is None:
            self.layer_manager['MAT'].mask()
            self.layer_manager['Surface'].mask()
        else:    
            f = np.zeros(self.ma.m*2, dtype=bool)
            f[ma_idx] = True
            self.layer_manager['MAT'].mask(f)

            # find indices of all surface points related to these mat points
            f = np.zeros(self.ma.m, dtype=bool)
            f[np.mod(ma_idx, self.ma.m)] = True
            f[self.ma.D['ma_qidx'][ma_idx]] = True
            self.layer_manager['Surface'].mask(f)
        
        self.viewerWindow.render()
    
    def addGraphWindow(self): 
        self.plotWindow = GraphWindow(master_app=self)
        # self.plotWindow.addLegend()
        self.plotWindow.show()
        # super(TestApp, self).run()

class ToolsWindow(ToolsDialog):
    def __init__(self, app):
        super(ToolsWindow, self).__init__(app, ui_path='sheeterator.ui', parent=None)

    def connectUI(self):
        super(ToolsWindow, self).connectUI()

        self.ui.comboBox_clusters.activated.connect(self.app.filter_cluster)
        self.ui.comboBox_sheets.activated.connect(self.app.filter_sheet)
        self.ui.groupBox_cluster.clicked.connect(self.app.toggle_selection)

        
class GraphWindow(PlotWidget):
    def __init__(self, master_app, parent=None):
        self.layer_manager = master_app.layer_manager
        self.master_app = master_app
        super(GraphWindow, self).__init__(parent)

    def show_sheet(self, ma_idx):
        # self.clear()

        # import ipdb/;ipdb.set_trace()
        # i = v.index
        vals=ma.D['ma_radii'][ma_idx]
        ## compute standard histogram
        y,x = np.histogram(vals, bins=50)

        ## Using stepMode=True causes the plot to draw two lines for each sample.
        ## notice that len(x) == len(y)+1
        # color = tuple(np.random.uniform(0.3,1.0,3)*255) + (255,)
        color = (0,0,200,255)
        self.plot(x, y, stepMode=True, fillLevel=0, pen={'color': color, 'width': 2}, name='radius '+str(42), clear=True)
        
        vals=ma.D['ma_theta'][ma_idx]
        ## compute standard histogram
        y,x = np.histogram(vals, bins=50)

        ## Using stepMode=True causes the plot to draw two lines for each sample.
        ## notice that len(x) == len(y)+1
        # color = tuple(np.random.uniform(0.3,1.0,3)*255) + (255,)
        color = (0,200,0,255)
        self.plot(x, y, stepMode=True, fillLevel=0, pen={'color': color, 'width': 2}, name='theta '+str(42))
    


def view(ma, vid):
    # ref_count = timeit(count_refs)
    max_r=190.
    # ma.g = ma.D['ma_segment_graph']
    
    c = TestApp(ma)   

    layer_ma = c.add_layer(LinkedLayer(name='MAT'))
    layer_s = c.add_layer(LinkedLayer(name='Surface'))
    layer_misc = c.add_layer(LinkedLayer(name='Other'))

    layer_s.add_data_source(
        name = 'Surface points',
        opts=['splat_disk', 'with_normals'],
        points=ma.D['coords'], 
        normals=ma.D['normals'],
    )
    layer_s.add_data_source_line(
      name = 'Surface normals',
      coords_start = ma.D['coords'] + ma.D['normals'],
      coords_end = ma.D['coords'],
      color = (1,1,0)
    )

    for v in g.vs():
		ma.D['ma_segment'][ v['ma_idx'] ] = v.index
    # f =ma.D['ma_segment'] != 0
    layer_ma.add_data_source(
        name = 'MAT points',
        opts=['splat_point', 'with_intensity'],
        points=ma.D['ma_coords'], 
        category=ma.D['ma_segment'].astype(np.float32),
        colormap='random'
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

    
    layer_misc.add_data_source(
        name = 'Unsegmented',
        opts=['splat_point','fixed_color', 'blend'],
        points=ma.D['coords'],
        color=(.6,.6,.6)
    )
    layer_misc.mask(ma.D['ma_segment'][:ma.m] == 0)     

    # c.viewerWindow.center_view(center=np.mean(ma.D['coords'][f_s], axis=0))
    c.run()

if __name__ == '__main__':
    vids = 4
    #box: 4,6,8 (3,10 broken in part)
    #sloped box: 7
    #simple gable: 9
    #exterior: 5
    if len(sys.argv)>1:
        # vids = [int(sys.argv[-1])]
        INFILE = sys.argv[-1]
    else:
    # INFILE = "/Users/ravi/git/mat_util/Random3Dcity/NPY"
        INFILE = "/Users/ravi/git/mat_util/test_cases/sloped_gable/NPY"
    datadict = io_npy.read_npy(INFILE)
    ma = MAHelper(datadict, origin=True)

    g = ma.D['ma_segment_graph']
    for v in g.vs:
        v['up_angle'] = np.sign(v['ma_bisec_mean'])[2] * np.arccos(np.dot(v['ma_bisec_mean'], [0,0,1] ))

    view(ma, vids)