import numpy as np
from graph import *
from  map import *
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

def line_intersect(l1,l2):
    # assuming there is an intersection and the lines are not parallel
    l2m = (l1.t[0] * (l2.r[1]-l1.r[1]) + l1.t[1]*l1.r[0] - l2.r[0]*l1.t[1]) / (l2.t[0]*l1.t[1] - l1.t[0]*l2.t[1])
    return l2.r + l2m*l2.t

def angle(a, b):
    a = a/np.linalg.norm(a)
    b = b/np.linalg.norm(b)
    return np.arccos(np.inner(a, b))

def gf_flatcube_top(master_graph, mapping, ma, ground_level=0):
    """compute cube geometry for matched graph based on the mapping with library graph
         v0
        /|\
       / | \
     v1  |  v3
       \ | /
        \|/
         v2
        each vertex is a sheet
        each sheet contributes to the top plane
        each sheet has one other side plane
    """
    g_flatcube_top = get_graph_library()['flatcube_top']
    g = master_graph
    v0, v1, v2, v3 = mapping

    # for each sheet find two sets of surface points, each corresponding to a distict plane. Use cross-product of spoke vectors with avg bisector of sheet to decide
    pointsets = {}
    coord_idx = []
    for vid in mapping:
        # obtain set of surface points (regardless of what is in/out)
        bisec_mean = g.vs[vid]['ma_bisec_mean']
        ma_idx = g.vs[vid]['ma_idx']
        s_idx = np.mod(ma_idx, ma.m)

        # collect normals and flip outward if necessary
        normals_1 = []
        normals_2 = []
        for ma_id, s_id in zip(ma_idx, s_idx):
            c = ma.D['ma_coords'][ma_id]
            # b = ma.D['ma_bisec'][ma_id]
            # n1 = ma.D['normals'][s_id]
            f1 = c-ma.D['coords'][s_id]
            # if np.inner(f1,n1) > 0:
            #     n1*=-1
            # n2 = ma.D['normals'][ma.D['ma_qidx'][ma_id]]
            f2 = c-ma.D['coords'][ma.D['ma_qidx'][ma_id]]
            # if np.inner(f2,n2) > 0:
            #     n2*=-1
            normals_1.append(f1/np.linalg.norm(f1))
            normals_2.append(f2/np.linalg.norm(f2))

        idx = np.concatenate([s_idx, ma.D['ma_qidx'][ma_idx]])

        km = KMeans(n_clusters = 2)
        labels = km.fit_predict(normals_1+normals_2)
        pts_0 = set(idx[labels==0])
        pts_1 = set(idx[labels==1])
        coord_idx.extend(idx)

        # print labels
        pointsets[vid] = [x/np.linalg.norm(x) for x in km.cluster_centers_], (pts_0, pts_1)
        # return ma_idx, normals, labels

    # return pointsets
    # aggregate points for common planes based on graph topology

    ## top Plane


    # print [len( x[1][0] | x[1][1] ) for x in pointsets.values()]

    c0,c2 = pointsets[v0][0], pointsets[v2][0]
    angle_thres = 20
    print "c0[0], c2[0]",
    print math.degrees(angle(c0[0], c2[0]))
    print "c0[0], c2[1]",
    print math.degrees(angle(c0[0], c2[1]))
    print "c0[1], c2[0]",
    print math.degrees(angle(c0[1], c2[0]))
    print "c0[1], c2[1]",
    print math.degrees(angle(c0[1], c2[1]))
    if angle(c0[0], c2[0]) < math.radians(angle_thres):
        up = c0[0]
        plane_top_pts = pointsets[v0][1][0] | pointsets[v2][1][0]
        plane_0_pts = pointsets[v0][1][1]
        plane_2_pts = pointsets[v2][1][1]
    elif angle(c0[0], c2[1]) < math.radians(angle_thres):
        up = c0[0]
        plane_top_pts = pointsets[v0][1][0] | pointsets[v2][1][1]
        plane_0_pts = pointsets[v0][1][1]
        plane_2_pts = pointsets[v2][1][0]
    elif angle(c0[1], c2[0]) < math.radians(angle_thres):
        up = c0[1]
        plane_top_pts = pointsets[v0][1][1] | pointsets[v2][1][0]
        plane_0_pts = pointsets[v0][1][0]
        plane_2_pts = pointsets[v2][1][1]
    else:
        up = c0[1]
        plane_top_pts = pointsets[v0][1][1] | pointsets[v2][1][1]
        plane_0_pts = pointsets[v0][1][0]
        plane_2_pts = pointsets[v2][1][0]

    if angle(up, pointsets[v1][0][0]) < math.radians(angle_thres):
        plane_top_pts |= pointsets[v1][1][0]
        plane_1_pts = pointsets[v1][1][1]
    else:
        plane_top_pts |= pointsets[v1][1][1]
        plane_1_pts = pointsets[v1][1][0]


    if angle(up, pointsets[v3][0][0]) < math.radians(angle_thres):
        plane_top_pts |= pointsets[v3][1][0]
        plane_3_pts = pointsets[v3][1][1]
    else:
        plane_top_pts |= pointsets[v3][1][1]
        plane_3_pts = pointsets[v3][1][0]

    # compute for each aggregate set of surface points a Plane
    plane_top = Plane( [Point( ma.D['coords'][pi]) for pi in plane_top_pts] )
    plane_0 = Plane( [Point( ma.D['coords'][pi]) for pi in plane_0_pts] )
    plane_1 = Plane( [Point( ma.D['coords'][pi]) for pi in plane_1_pts] )
    plane_2 = Plane( [Point( ma.D['coords'][pi]) for pi in plane_2_pts] )
    plane_3 = Plane( [Point( ma.D['coords'][pi]) for pi in plane_3_pts] )

    # intersect planes using graph topology, compute the vertices of the roof plane, extrude down to z=0
    line_0 = plane_top.intersection(plane_0)
    line_1 = plane_top.intersection(plane_1)
    line_2 = plane_top.intersection(plane_2)
    line_3 = plane_top.intersection(plane_3)

    line_01 = plane_0.intersection(plane_1)
    line_12 = plane_1.intersection(plane_2)
    line_23 = plane_2.intersection(plane_3)
    line_30 = plane_3.intersection(plane_0)

    q0 = line_intersect(line_0, line_01)
    q1 = line_intersect(line_1, line_12)
    q2 = line_intersect(line_2, line_23)
    q3 = line_intersect(line_3, line_30)

    ground_level = ma.D['coords'][coord_idx][:,2].min()
    q0_ = q0.copy()
    q0_[2] = ground_level
    q1_ = q1.copy()
    q1_[2] = ground_level
    q2_ = q2.copy()
    q2_[2] = ground_level
    q3_ = q3.copy()
    q3_[2] = ground_level
    # print q0
    # print q1
    # print q2
    # print q3

    # coords = np.empty((6,3), dtype=np.float32)
    normals = np.empty((30,3), dtype=np.float32)

    coords = np.array([     q0,q1,q2, q0,q2,q3,
                            q0,q0_,q1_, q0,q1_,q1,
                            q1,q1_,q2, q2,q1_,q2_,
                            q2,q2_,q3_, q3,q2,q3_,
                            q0,q3,q3_, q0,q3_,q0_ ], dtype=np.float32)
    normals[0:6] = plane_top.n
    normals[6:12] = plane_0.n
    normals[12:18] = plane_1.n
    normals[18:24] = plane_2.n
    normals[24:30] = plane_3.n

    return coords, normals, (plane_top_pts, plane_0_pts, plane_1_pts, plane_2_pts, plane_3_pts)

def build_map(g, ma):
    p_angle_match = 5
    p_angle_parallel = 165
    p_angle_converge = 170

    # create MAP datastructure
    m = Map()

    # for each sheet find two sets of surface points, each corresponding to a distinct plane. Store results in sheet
    vids_coplanar = [] # list to keep sheets that we ignore because they support parallel planes (separation angle=180deg)
    vid_map = {}
    nid_cnt = 0
    for i,v in enumerate(g.vs):
        ma_idx = v['ma_idx']
        if v['ma_theta_mean'] > math.radians(p_angle_converge):
            vids_coplanar.append(v.index)
            continue
        # elif ma.D['ma_radii'][ma_idx].min() > 0.8: # ignore sheets that are not fully going into the building edges
        #     vids_coplanar.append(v.index)
        #     continue
        
        # obtain set of surface points (regardless of what is in/out)
        s_idx = np.mod(ma_idx, ma.m)

        # collect spokes
        c = ma.D['ma_coords'][ma_idx]
        f1 = c-ma.D['coords'][s_idx]
        f2 = c-ma.D['coords'][ma.D['ma_qidx'][ma_idx]]
        
        f1_norm = np.linalg.norm(f1, axis=1)[:,None]
        f2_norm = np.linalg.norm(f2, axis=1)[:,None]
        spokes = np.concatenate([f1/f1_norm, f2/f2_norm])

        idx = np.concatenate([s_idx, ma.D['ma_qidx'][ma_idx]])

        km = KMeans(n_clusters = 2)
        labels = km.fit_predict(spokes)

        locked = ma.D['ma_radii'][ma_idx].min() > 0.8

        m.add_node(
            {'s_idx':idx[labels==0], 'spoke_cluster_center':km.cluster_centers_[0], 'locked':locked, 'vid':i}, 
            {'s_idx':idx[labels==1], 'spoke_cluster_center':km.cluster_centers_[1], 'locked':locked, 'vid':i}
        )
        vid_map[i]=nid_cnt
        nid_cnt += 1

    
    # label each edge with corresponding spoke sets from incident vertices (ie. an edge can be linked to exactly one plane)
    for ge in list(g.es):
        # convert source/target idx to halfnode idx in map
        source, target = ge.tuple
        # skip edges that connect to sheets supporting parallel planes
        if source in vids_coplanar or target in vids_coplanar:
            continue
        source, target = [vid_map[x] for x in ge.tuple]
        hs0, ht0 = source*2+0, target*2+0
        hs1, ht1 = source*2+1, target*2+1

        # compute all angles between spokesets of source and target sheet, also find minimum angle and corresponding label pair
        angles = {}
        angle_min = 10 #radians
        label_pair_min = None
        for sn in (hs0, hs1):
            for tn in (ht0,ht1):
                a = angle(m.ns[sn]['spoke_cluster_center'], m.ns[tn]['spoke_cluster_center'])
                angles[(sn,tn)] = a
                if a < angle_min: 
                    angle_min = a
                    label_pair_min = (sn,tn)

        # make sure we find a good plane correspondence
        if not angle_min < math.radians(p_angle_match):
            continue

        # store on this edge a dict with for each vertex_id the corresponding spoke_cluster_label
        s_min, t_min = label_pair_min 
        other_pair = other(s_min, [hs0,hs1]), other(t_min, [ht0,ht1])
        s_other, t_other = other_pair

        if angles[other_pair] > math.radians(p_angle_parallel): # what threshold is reasonable here? 175 degrees seems to be too high
            continue # skip parallel edges for now
            kind='parallel'
        else:
            kind='intersect'

        m.add_edge(s_min, t_min, kind='match')
        # m.add_edge(s_other, t_other, kind) # we don't need these it seems

    # remove locked sheets and rewire local topology
    to_delete = []
    for hn in m.ns:        
        if hn['locked']:
            to_delete.append(hn)
    for hn in to_delete:
        m.contract_halfnode(hn)

    # find halfnodes with degree > 2 and try to get rid of some edges
    # first get rid of dangling edges,
    #...
    # then remove edges that split loops
    for hn in m.ns:
        if 'match' in hn.edges:
            if len(hn.edges['match']) > 2:
                for e, hno in hn.neighbours():
                    if len(hno.edges['match']) > 2:
                        m.delete_edge(e)


    # orient all plane cycles ccw

    return m


def polyhedral_reconstruct(m, ma):
    """Reconstruct polyhedral model for this mat-component/cluster assuming it represents some suitable object from the inside"""

    # m = build_map(g, ma)
    # Find hnodes that are linked to the same plane, using connecting edges of 'match' kind
    # planes = {}
    # plane_id = 0
    N = set(m.ns)
    while len(N) > 0:
        n = N.pop()
        f = m.add_face(halfnode=n, cycle=True)
        next = None
        for next in n.cycle(kind='match', direction=0): # assuming here we only have duplicate edges of intersect kind
            next.face = f
            N.discard(next)
            # print next, next.edges['match']
        if next != n: # see if there are some in the other direction
            f['cycle'] = False
            n.face = f
            # note: following statement may crash the function
            # for next in n.cycle(kind='match', direction=1):
            #     next.face = f

    # introduce virtual vertices and edges for missing planes
    # we assume that all present edges are properly connected and labeled now
    # set all current vertices as not virtual:
    # g.vs['is_virtual'] = False
    # for e in g.es:
    #     # if there is an edge now there must be a match in the spoke_cluster_map, so we take that for granted
    #     if 'intersect' in e['spoke_cluster_map']:
    #         # check if the intersecting spoke sets of source target sheets both reference two planes 
    #         # if not we should introduce virtual vertex + 2 virtual edges and create a new plane/update references on an existing plane_id
    #         (source_vid, source_spokecluster), (target_vid, target_spokecluster) = e['spoke_cluster_map']['intersect'].items()
    #         source_v = g.vs[source_vid]
    #         target_v = g.vs[target_vid]
    #         # source_spokecluster, target_spokecluster = inverse_spokeset_pair(source_spokeclusterm, target_spokecluster)

    #         source_plane_list = source_v['spoke_cluster_planeids'][source_spokecluster]
    #         target_plane_list = target_v['spoke_cluster_planeids'][target_spokecluster]

    #         if len(source_plane_list) == 2 or len(target_plane_list) == 2:
    #             continue # no problem here, go to the next edge

    #         source_s_idx = get_surface_points(source_v, ma)[source_v['spoke_cluster_labels']==source_spokecluster]
    #         target_s_idx = get_surface_points(target_v, ma)[target_v['spoke_cluster_labels']==target_spokecluster]
    #         s_idx = np.concatenate([source_s_idx, target_s_idx])#.tolist()
    #         labels = np.empty(len(s_idx), dtype=int)
    #         labels[:len(source_s_idx)] = 0
    #         labels[len(source_s_idx):] = 1

    #         source_eid = g.ecount()
    #         target_eid = source_eid+1

    #         print len(source_plane_list), len(target_plane_list)
    #         if len(source_plane_list) == 0 and len(target_plane_list) == 0:
    #             source_plane = plane_id
    #             target_plane = plane_id+1
    #             plane_id += 2

    #             # create two new planes
    #             planes[source_plane] = [source_eid]
    #             planes[target_plane] = [target_eid]
    #             source_plane_list.append(source_plane)
    #             target_plane_list.append(target_plane)
    #             # add references from source_v/target_v
    #         elif len(source_plane_list) == 1 and len(target_plane_list) == 0:
    #             source_plane = source_plane_list[0]
    #             target_plane = plane_id
    #             plane_id += 1
    #             planes[source_plane].append(source_eid)
    #             planes[target_plane] = [target_eid]
    #             target_plane_list.append(target_plane)
    #         elif len(source_plane_list) == 0 and len(target_plane_list) == 1:
    #             target_plane = target_plane_list[0]
    #             source_plane = plane_id
    #             plane_id += 1
    #             planes[source_plane] = [source_eid]
    #             planes[target_plane].append(target_eid)
    #             source_plane_list.append(source_plane)
    #         else: # both are length 1
    #             source_plane = source_plane_list[0]
    #             target_plane = target_plane_list[0]
    #             planes[source_plane].append(source_eid)
    #             planes[target_plane].append(target_eid)
                
    #         virtual_vid = g.vcount()
    #         g.add_vertex(
    #             is_virtual = True,
    #             ma_idx = s_idx, # NOTE: for virtual vertices we actually store surfacepoint ids here! 
    #             spoke_cluster_planeids = ([source_plane],[target_plane]),
    #             # spoke_cluster_centers=, # not really needed now 
    #             spoke_cluster_labels=labels
    #         )

    #         # here we should also add the intersections...? maybe not required
            
    #         g.add_edge(virtual_vid, source_vid, spoke_cluster_map={
    #             'match':{virtual_vid:0, source_vid:source_spokecluster}
    #         })
            
    #         g.add_edge(virtual_vid, target_vid, spoke_cluster_map={
    #             'match':{virtual_vid:1, target_vid:target_spokecluster}
    #         })

    # fit each plane through its supporting points
    for f in m.fs:
        if f['cycle']:
            s_idx = set()
            for n in f.halfnode.cycle('match'):
                s_idx |= set(n['s_idx'])
            f['plane'] = Plane( [Point(x) for x in ma.D['coords'][list(s_idx)]] )
        else:
            f['plane'] = Plane( [Point(x) for x in ma.D['coords'][list(f.halfnode['s_idx'])]] )

    # orient all plane cycles ccw

    # compute line of intersection for each halfnode pair

    # traverse vertex cycles and compute vertices

    # reconstruct planes from vertices
    # import ipdb; ipdb.set_trace()
    planes = []
    for f in m.fs:
        if f['cycle']:
            p = f['plane']
            lines = []
            for hn in f.halfnode.cycle('match'):
                lines.append(p.intersection(hn.twin.face['plane']))
            vertices = []
            vertices.append( line_intersect(lines[-1], lines[0]) )
            for i in range(len(lines)-1):
                vertices.append(line_intersect(lines[i], lines[i+1]))

            
            # fan-like triangulation
            coords = np.empty((3*(len(vertices)-2),3))
            normals = np.empty((3*(len(vertices)-2),3))
            for i in range(len(vertices)-2):
                coords[3*i] = vertices[0]
                coords[3*i+1] = vertices[i+1]
                coords[3*i+2] = vertices[i+2]
                normals[3*i:3*i+3] = p.n

            if len(vertices) ==5 :
                print vertices
                import ipdb; ipdb.set_trace()

            
            planes.append((coords, normals))
    # print planes
        # print vertices, coords
    return planes
