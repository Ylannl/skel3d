import numpy as np
from graph import *
import random
from geom3d import *

from sklearn.cluster import KMeans

def ransac_plane_fit(point_array, point_indices, threshold=0.05, n_planes=2, max_iterations=10):
    # choose 3 random points

    idx = list(point_indices)

    # we are searching for n_planes solutions, but we may only find a lower number of planes    
    solutions = []
    for i in range(n_planes):
        if len(idx) < 3:
            break
        iteration_cnt = 0
        candidate_solutions = []
        while iteration_cnt < max_iterations:
            random.shuffle(idx)
            seeds = idx[:3]
            # import ipdb; ipdb.set_trace()
            try: 
                candidate_plane = Plane( [Point( point_array[pi]) for pi in seeds] )
            except ValueError:
                continue

            inliers = seeds
            for pi in idx[3:]:
                d = candidate_plane.distance_to( Point(point_array[pi]) )
                if d < threshold:
                    inliers.append(pi)
            candidate_solutions.append( (candidate_plane, inliers) )

            iteration_cnt+=1
        # determine best Plane
        candidate_solutions.sort(key=lambda tup: len(tup[1]), reverse=True)
        plane, inliers = candidate_solutions[0]

        idx = list(set(idx) - set(inliers))
        solutions.append(plane)
    
    return solutions

def derive_tetra(plane1, plane2, point_indices, ma):
    line = plane1.intersection(plane2)

    line_range = []
    for c in ma.D['ma_coords'][point_indices]:
        p = Point(c)
        s = np.inner(line.t,p.r - line.r)
        line_range.append(s)
    s_mi, s_ma = min(line_range), max(line_range)
    v0,v1 = line.t*s_mi+line.r, line.t*s_ma+line.r

    pi_maxr =  np.argmax(ma.D['ma_radii'][point_indices])
    p_maxr = ma.D['ma_coords'][point_indices][pi_maxr]

    v2 = plane1.project(Point(p_maxr)).r
    v3 = plane2.project(Point(p_maxr)).r

    # center v2,v3 nicely to be halfway across v0,v1
    s = np.inner(line.t,v2 - line.r)
    move_to_center = (s_mi + (s_ma-s_mi)/2 - s)*line.t
    v2 += move_to_center
    v3 += move_to_center

    return np.array([v0, v1, v2, v3], dtype=np.float32)

# def get_graph_library():
#     graph_library = {}
#     g = igraph.Graph()
#     g.add_vertices(4)
#     g.add_edges([(0,1), (1,2), (2,3), (3,0), (0,2)])
#     graph_library['flatcube_top'] = g

#     return graph_library

def gf_flatcube_top(master_graph, mapping, ma):
    """compute cube geometry for matched graph based on the mapping with library graph
         v1
        /|\
       / | \
     v2  |  v3
       \ | /
        \|/
         v0
        each vertex is a sheet
        each sheet contributes to the top plane
        each sheet has one other side plane
    """
    g_flatcube_top = get_graph_library()['flatcube_top']
    g = master_graph
    v0, v1, v2, v3 = mapping
    
    # for each sheet find two sets of surface points, each corresponding to a distict plane. Use cross-product of spoke vectors with avg bisector of sheet to decide
    pointsets = {}
    for vid in mapping:
        # obtain set of surface points (regardless of what is in/out)
        bisec_mean = g.vs[vid]['ma_bisec_mean']
        ma_idx = g.vs[vid]['ma_idx']

        # collect normals and flip outward if necessary
        normals = []
        for ma_id, s_id in zip(ma_idx, np.mod(ma_idx, ma.m)):
            c = ma.D['ma_coords'][ma_id]
            # b = ma.D['ma_bisec'][ma_id]
            n1 = ma.D['normals'][s_id]
            # f1 = c-ma.D['coords'][s_id]
            # if np.inner(f1,n1) < 0:
                # n1*=-1
            n2 = ma.D['normals'][ma.D['ma_qidx'][ma_id]]
            # f2 = c-ma.D['coords'][ma.D['ma_qidx'][ma_id]]
            # if np.inner(f2,n2) < 0:
                # n2*=-1
            normals.extend([n1,n2])

        km = KMeans(n_clusters = 2)
        labels = km.fit_predict(normals)
        # import ipdb; ipdb.set_trace()
        # print labels
        pointsets[vid] = km.cluster_centers_, labels
        # return ma_idx, normals, labels
            

    # aggregate points for common planes based on graph topology

    ## top Plane
    pointsets[v0]

    # compute for each aggregate set of surface points a Plane

    # intersect planes using graph topology, compute the vertices of the roof plane, extrude down to z=0 
