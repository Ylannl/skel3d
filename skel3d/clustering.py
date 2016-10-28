from __future__ import division
import math, sys
from time import time
import numpy as np

from scipy.stats import linregress

# from skel3d.io import npy
from .util import angle
# from skel3d.graph import *
# from skel3d.segmentation import *
# from skel3d.polyhedralise import *

from .geom3d import *

def get_clusters(ma, min_count = 50):
    ## Remove insignificant edges in graph
    contract_thres = 15 #self.dialog.ui.doubleSpinBox_contractthres.value()
    # g = g.subgraph(g.vs.select(ma_theta_mean_lt=math.radians(100), up_angle_gt=math.radians(40)))
    master_g = ma.D['ma_segment_graph']
    master_g = master_g.subgraph_edges(master_g.es.select(adj_count_gt=min_count))
    # contract_edges(g, contract_thres)

    ## Update segment indices based on graph ids
    ma.D['ma_segment'] = np.zeros(ma.m*2,dtype=np.int64)
    for v in master_g.vs:
        ma.D['ma_segment'][v['ma_idx']] = v.index
        v['s_id'] = v.index
    # flipdic = find_flip_relations(ma)
    
    # self.graphs = []
    # graphlib = get_graph_library()
    # for mapping in g.get_subisomorphisms_vf2(graphlib['flatcube_top']):
    #     self.graphs.append(g.subgraph(mapping))
    
    ## Find clusters, ie. connected component analysis
    ma.D['ma_clusters'] = master_g.clusters().subgraphs()


def classify_clusters(ma):
    i=0
    for g in ma.D['ma_clusters']:
        adj_rel_start = []
        adj_rel_end = []

        if 0<g.ecount():#<1000:

            ## Classify clusters as inside or outside
            # attempt to distinghuish between interior and exterior sheets based on vertical component of bisectors
            color = (1.,0.,0.)
            name = 'cluster {}'.format(i)
            zsum = 0
            ma_idx = np.concatenate(g.vs['ma_idx'])
            
            # select points with good sepangle (eg not 180deg)
            f = np.intersect1d(np.argwhere(ma.D['ma_theta'] < math.pi*0.9), ma_idx)
            # sort on z component of bisector, take 10% lowest and highest, compute ratio of means of that
            l = len(f)
            bz = ma.D['ma_bisec'][f][:,2]
            bz_sort = np.sort(bz )
            b_bot = np.mean(bz_sort[:l//10]) 
            b_top = np.mean(bz_sort[l-l//10:])
            b_ratio = abs(b_top)/abs(b_bot)

            # for v in g.vs:
            #     zsum += v['ma_bisec_mean'][2]#*len(v['ma_idx'])
            # zsum = np.prod(np.sum(np.array(g.vs['ma_bisec_mean'])[:,2], np.array([len(idx) for idx in g.vs['ma_idx']]) ), axis=1)
            # use ration to decide if int or ext cluster. Now
            classname = "exterior"
            if b_ratio > 1.05 or b_bot > 0: # completely `closed` or bounded clusters should have ratio closed to one (only occurs in artifical datasets, since DSM never closed)
                color = (0.,1.,0.)
                classname = "interior"
                if np.mean(g.vs['ma_theta_mean']) > math.pi/4: # artificial/building structures typically have a large sepangle (compared to terrain features)
                    color = (1.,1.,0.)
                    classname += " (building)"

            elif b_ratio > 0.95:
                color = (0.5,0.5,0.5)
                classname = ""

            g['classification'] = classname



def analyse_cluster(ma, g):
    def cluster_spokes(ma,ma_idx):
        """Determine orientation of spoke-pair for each ma point"""
        cross = np.cross(ma.D['ma_f1'][ma_idx], ma.D['ma_f2'][ma_idx])
        l = Line([Point(v) for v in cross])
        # l.t is a unit vector in one direction of the line
        x = np.empty(len(cross))
        for i, v in enumerate(cross):
            x[i] = np.dot(l.t,v)

        return x > 0, cross

    for v in g.vs:
        ma_idx = v['ma_idx']
        r = np.mean(ma.D['ma_radii'][ma_idx], axis=0)
        theta = np.mean(ma.D['ma_theta'][ma_idx], axis=0)
        b = np.mean(ma.D['ma_bisec'][ma_idx], axis=0)
        c = np.mean(ma.D['ma_coords'][ma_idx], axis=0)

        # define x-axis; it opints in the direction of convergence of the 2 planes this sheet is supporting and is 0 at that ridge
        xc = r/np.cos(theta/2)
        x = np.empty(len(ma_idx))
        for i in range(len(ma_idx)):
            x[i] = xc + (np.dot(c-ma.D['ma_coords'][ma_idx[i]], b))

        ## collect al measures
        Y = {}

        ## radius
        Y['radius'] = ma.D['ma_radii'][ma_idx]
        
        ## separation angle
        Y['ma_theta'] = ma.D['ma_theta'][ma_idx]

        ## bisector
        Y['ma_bisec_diff'] = np.empty(len(ma_idx))
        for i in range(len(ma_idx)):
            Y['ma_bisec_diff'][i] = angle(b, ma.D['ma_bisec'][ma_idx[i]])
        
        ## predicted separation angle based on radius and position along x-axis
        Y['ma_thetap'] = np.empty(len(ma_idx))
        for i in range(len(ma_idx)):
            Y['ma_thetap'][i] = 2*np.arccos(ma.D['ma_radii'][ma_idx[i]]/ x[i])

        ## compute 'normal' direction for sheet
        # cross product of spokes is perpendicular to bisector and tangent to sheet
        one_side, cross = cluster_spokes(ma, ma_idx)
        # align al crosses and compute average
        cross_align = cross 
        cross_align[~one_side] *= -1  
        # np.concatenate([cross[one_side], -1*cross[~one_side]])
        vec_coplanar = np.mean(cross_align, axis=0)
        # now compute this cross product to find a vector in the normal direction of the plane that we want to reconstruct
        n = np.cross(vec_coplanar, b)
        n = n / np.linalg.norm(n)
        # plane = Plane(pc, Line(pc, pn))

        # plane fit error
        Y['ma_planfit'] = np.empty(len(ma_idx))
        for i in range(len(ma_idx)):
            q = ma.D['ma_coords'][ma_idx[i]] - c
            q_on_n = np.dot(q,n)
            Y['ma_planfit'][i] = q_on_n

        # diff in cross
        Y['biseco_diff'] = np.empty(len(ma_idx))
        for i in range(len(ma_idx)):
            Y['biseco_diff'][i] = angle(vec_coplanar, cross_align[i])

        # this ratio should say something about the relative number of primary feature points on each side of the sheet
        n_oneside = one_side.sum()
        n_total = len(one_side)
        n_otherside = n_total-n_oneside
        if n_oneside > n_otherside:
            regularity_ratio = n_oneside/n_total
        else:
            regularity_ratio = n_otherside/n_total

        result = {}
        for key, y in Y.items():
            slope, intercept, r_value, p_value, std_err = linregress(x,y)
            result[key] = {}
            result[key]['slope']=slope
            result[key]['intercept']=intercept
            result[key]['r_value']=r_value
            result[key]['p_value']=p_value
            result[key]['std_err']=std_err
            result[key]['regularity_ratio']=regularity_ratio
        v['sheet_analysis'] = result


        

def grow_cluster(g, ma, min_flipcount = 10):
    # grow 1 flip sheet around this cluster
    to_add = []
    for v in g.vs:
        v['fliplinks'] = {}
        s_idx = np.concatenate([np.mod(v['ma_idx'], ma.m), ma.D['ma_qidx'][v['ma_idx']]])
        for s_id in s_idx:
            s_in = ma.D['ma_segment'][s_id]                        
            s_out = ma.D['ma_segment'][s_id+ma.m]

            s_other = s_in
            if s_in == v['s_id']:
                s_other = s_out
            
            if s_other == 0:
                continue

            if s_other in v['fliplinks']:
                v['fliplinks'][s_other] += 1
            else:
                v['fliplinks'][s_other] = 1
        for s_id, count in v['fliplinks'].items():
            if count > min_flipcount:
                to_add.append((v.index, s_id, count))
    for source, target, count in to_add:
        g.add_vertex(**master_g.vs[target].attributes())
        g.add_edge(source, g.vcount()-1, flip_count=count)