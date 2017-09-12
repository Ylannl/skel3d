from __future__ import division

allowed_kinds = 'intersect', 'match', 'parallel'

class Map(object):

    def __init__(self):
        self.ns = []
        self.es = []
        self.fs = []
        self.efs = []

    def add_nodes(self, n):
        for i in xrange(n):
            self.add_node()

    def add_node(self, kwargs1={}, kwargs2={}):
        shared = dict() # shared dictionary between pair of halfnodes
        hn1 = HalfNode(shared, **kwargs1)
        hn2 = HalfNode(shared, **kwargs2)
        hn1.twin = hn2
        hn2.twin = hn1
        self.ns.extend([hn1, hn2])
        return hn1, hn2

    def delete_node(self, node):

        # node = self.ns[nid]
        self.delete_halfnode(node)
        self.delete_halfnode(node.twin)

    def contract_halfnode(self, node):
        # only cares about match edges
        # add new edges between neighbours of this node if they connect through this node
        for e1, hn1 in node.neighbours():
            for e2, hn2 in node.neighbours():
                if hn1 != hn2:
                    # check if an edge exists between hn1 and hn2
                    if len(set(hn1.edges['match']) & set(hn2.edges['match'])) == 0:
                        # add that Edge
                        self.add_edge(hn1, hn2, 'match', 'unsure')

        # remove incident edges to this node
        for e, hn in node.neighbours():
            # remove references to e in hn
            hn.edges['match'].remove(e)
            # remove e from map
            self.es.remove(e)
            # delete e itself
            del e
        self.ns.remove(node)
        del node

    def delete_halfnode(self, hn):
        # note take care of refs from edges yourself!
        self.ns.remove(hn)
        
    def delete_edge(self, e):
        for hn in e.nodes:
            hn.edges['match'].remove(e)
        self.es.remove(e)
        del e

    def add_edges(self, node_pairs, kind):
        for s,t in node_pairs:
            self.add_edge(s,t, kind=kind)

    def add_edge(self, s_idx, t_idx, kind, other_kind):
        if type(t_idx) == int:
            s, t = self.ns[s_idx], self.ns[t_idx]
        else:
            s, t = s_idx, t_idx
        e = Edge(s, t, kind, other_kind)
        s.edges.setdefault(kind,[]).append(e)
        t.edges.setdefault(kind,[]).append(e)
        self.es.append(e)
        return e
        
    def add_face(self, link, **kwargs):
        f = Face(**kwargs)
        if type(link)==HalfNode:
            f.halfnode=link
            self.fs.append(f)
        elif type(link)==Edge:
            f.edge=link
            self.efs.append(f)
        return f

def other(this, pair):
    """assume pair is an indexable iterable of length 2 and this is an element in pair, this function will return the other element"""
    return pair[not pair.index(this)]

class HalfNode(object):

    def __init__(self, shared, face=None, **kwargs):
        self.twin = None
        self.shared = shared
        self.edges = {} # 
        self.face = face
        self.attributes = {}
        self.attributes.update(kwargs)

    def __getitem__(self, key):
        return self.attributes[key]
    
    def __setitem__(self, key, value):
        self.attributes[key] = value

    def cycle(self, kind): # direction = [0,1] eg start with first or second edge from this node
        """Does not care about edge direction during the cycle"""
        if not kind in self.edges.keys(): return
        next = self.edges[kind][1].nodes[1]
        while next!=self:
            yield next
            next = next.edges[kind][1].nodes[1]
        yield self
    
    def cycle_undirected(self, kind, direction=0): # direction = [0,1] eg start with first or second edge from this node
        """Does not care about edge direction during the cycle"""
        if not kind in self.edges.keys(): return
        follow_edge = self.edges[kind][direction]
        next = other(self, follow_edge.nodes)
        while next!=self:
            yield next
            if len(next.edges[kind])==1: return # there is no next edge
            if len(next.edges[kind])>2: return # there are too many edges to choose from
            follow_edge = other(follow_edge, next.edges[kind])
            next = other(next, follow_edge.nodes)
        yield self

    def flip(self, kind):
        self.edges[kind] = [self.edges[kind][1], self.edges[kind][0]]

    def neighbours(self, kind='match'):
        if not kind in self.edges.keys(): return
        for e in self.edges[kind]:
            yield e, other(self, e.nodes)

    # def __repr__(self):
    #     return "HalfNode edges:{}".format(self.edges)

class Edge(object):
    
    def __init__(self, source, target, kind, other_kind, face=None):
        self.nodes = (source, target)
        self.kind = kind # 'intersect' 'parallel' 'match'
        self.other_kind = other_kind # 'intersect' 'parallel' 'match'
        self.face = face # i.e. the face that is a vertex cycle 

    def flip(self):
        self.nodes = self.nodes[1], self.nodes[0] 

    # def __repr__(self):
    #     return "Edge [{}] nodes:{}".format(self.kind, self.nodes)

class Face(object):

    def __init__(self, **kwargs):
        self.attributes = {}
        self.attributes.update(kwargs)

    def __getitem__(self, key):
        return self.attributes[key]

    def __setitem__(self, key, value):
        self.attributes[key] = value

    # def __repr__(self):
    #     return "Face node:{}".format(self.halfnodes)

def test():
    m = Map()
    m.add_nodes(3)
    m.add_edges([(0,2),(2,4),(4,0)], kind='intersect')
    n = m.ns[0]
    for x in n.cycle(kind='intersect'):
        print(x)


import math
import numpy as np
from .util import angle
from .geom3d import Plane, Point
def build_map(g, ma):
    def surface_match(n1, n2):
        s_idx = np.concatenate([n1['s_idx'], n2['s_idx']])
        l1 = len(n1['s_idx'])
        l2 = len(n2['s_idx'])
        if l1 > l2:
            bigone = n1
            smallone = n2
        else:
            bigone = n2
            smallone = n1

        plane = Plane( [Point(p) for p in ma.D['coords'][bigone['s_idx']]] )
        threshold = 0.10
        n_inliers = 0
        for pi in smallone['s_idx']:
            d = plane.distance_to( Point(ma.D['coords'][pi]) )
            if d < threshold:
                n_inliers += 1 
        
        # import ipdb; ipdb.set_trace()
        print (n_inliers/len(smallone['s_idx']))
        return n_inliers/len(smallone['s_idx']) > 0.9


    p_angle_match = 5
    p_angle_parallel = 165
    p_angle_converge = 170

    # create MAP datastructure
    m = Map()
    print("Creating nodes...")
    try:
        # for each sheet find two sets of surface points, each corresponding to a distinct plane. Store results in sheet
        vids_coplanar = [] # list to keep sheets that we ignore because they support parallel planes (separation angle=180deg)
        vid_map = {} # mapping we need to keep because we are skipping some nodes in the graph and indexing is based on position in list 
        nid_cnt = 0
        for i,v in enumerate(g.vs):
            ma_idx = v['ma_idx']
            if v['ma_theta_mean'] > math.radians(p_angle_converge):
                vids_coplanar.append(v.index)
                continue
            # if len(ma_idx) < 30: # ignore small sheets, eg. often caused by oversegmentation, should actually be merged
            #     vids_coplanar.append(v.index)
            #     print('skipping small sheet...')
            #     continue
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


            locked = ma.D['ma_radii'][ma_idx].min() > 0.8

            m.add_node(
                {'s_idx':v['s_idx1'], 'spoke_cluster_center':v['spoke_cluster_center1'], 'locked':locked, 'vid':i}, 
                {'s_idx':v['s_idx2'], 'spoke_cluster_center':v['spoke_cluster_center2'], 'locked':locked, 'vid':i}
            )

            # km = KMeans(n_clusters = 2)
            # labels = km.fit_predict(spokes)
            # m.add_node(
            #     {'s_idx':idx[labels==0], 'spoke_cluster_center':km.cluster_centers_[0], 'locked':locked, 'vid':i}, 
            #     {'s_idx':idx[labels==1], 'spoke_cluster_center':km.cluster_centers_[1], 'locked':locked, 'vid':i}
            # )
            vid_map[i]=nid_cnt
            nid_cnt += 1
        # import ipdb; ipdb.set_trace()

        yield m
        print("Creating edges...")

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

            # import ipdb; ipdb.set_trace()
            if angles[other_pair] > math.radians(p_angle_parallel): # what threshold is reasonable here? 175 degrees seems to be too high
                continue # skip parallel edges for now
                other_kind='parallel'
            # elif angles[other_pair] < math.radians(30):
            elif surface_match(m.ns[s_other], m.ns[t_other]):
                other_kind='match'
            else:
                other_kind='intersect'

            m.add_edge(s_min, t_min, kind='match', other_kind=other_kind)
            # m.add_edge(s_other, t_other, other_kind) # we don't need these it seems
        
        yield m
        print("Merging matching sheets ...")

        to_delete = []
        for e in m.es:        
            if e.other_kind=='match':
                to_delete.append(e.nodes[0])
        for hn in to_delete:
            m.contract_halfnode(hn)
        
        yield m
        print("Removing locked sheets...")

        # remove locked sheets and rewire local topology
        to_delete = []
        for hn in m.ns:        
            if hn['locked']:
                to_delete.append(hn)
        for hn in to_delete:
            m.contract_halfnode(hn)

        # find halfnodes with degree > 2 and try to get rid of some edges
        # first get rid of dangling edges,
        # to_delete = []
        # for hn in m.ns:
        #     if 'match' in hn.edges:
        #         if len(hn.edges['match']) < 2:
        #             for e in hn.edges['match']:
        #                 print e
        #                 m.delete_edge(e)
        #             to_delete.append(hn)

        # for hn in to_delete:
        #     m.delete_halfnode(hn)

        
        yield m
        print("Removing edges that split loops...")
        # then remove edges that split loops
        for hn in m.ns:
            if 'match' in hn.edges:
                if len(hn.edges['match']) > 2:
                    for e, hno in hn.neighbours():
                        if len(hno.edges['match']) > 2:
                            m.delete_edge(e)

        yield m
        print("Orienting cycles...")
        
        # attempt to find cycles and consistently direct all edges/halfnodes in a cycle. Incomplete cycles are also assigned faces.
        N = set(m.ns)
        while len(N) > 0:
            n = N.pop()
            f = m.add_face(link=n, cycle=True)
            n.face = f

            next_edge = n.edges['match'][1] # assume id 0 is the previous edge, and id 1 the next edge
            if next_edge.nodes[1] == n:
                next_edge.flip()
            next_node = next_edge.nodes[1]
            while next_node!=n:
                N.discard(next_node)
                next_node.face = f
                if len(next_node.edges['match'])!=2: # there is no next edge or there are too many edges to choose from
                    f['cycle'] = False 
                    break
                # ensure the current next_edge is the current next_nodes' previous edge 
                if next_node.edges['match'][0] != next_edge:
                    next_node.flip('match')
                next_edge = next_node.edges['match'][1]
                # flip the next edge if not consistent with direction of current edge (eg current node should be the source of the next edge)
                if next_edge.nodes[1] == next_node:
                    next_edge.flip()
                next_node = next_edge.nodes[1]
            
        # orient all cycles ccw using cross product of consecutive differences of bisector, Then see make that face the same direction as spoke vector
        # thus from hereon we are dealing with directed edges
        for f in m.fs:
            if not f['cycle']: continue
            e_previous, e_next = f.halfnode.edges['match']

            hnodes = [hn for hn in f.halfnode.cycle('match')]
            spoke = hnodes[0]['spoke_cluster_center']
            b0, b1, b2 = [g.vs[hn['vid']]['ma_bisec_mean'] for hn in hnodes[:3]] 
            dob1 = b1-b0
            dob2 = b2-b1
            # import ipdb; ipdb.set_trace()
            if np.dot(np.cross(dob1,dob2), spoke) > 0:
                for hn in hnodes:
                    hn.flip('match')
                    hn.edges['match'][0].flip()

        yield m
    except Exception as e:
        raise e
        # import ipdb; ipdb.set_trace()

class MAPCreator(object):
    p_angle_match = 5
    p_angle_parallel = 165
    p_angle_converge = 170

    def __init__(self, g, ma):
        self.g = g
        self.ma = ma
        self.m = Map()

    def create_nodes(self):
        # for each sheet find two sets of surface points, each corresponding to a distinct plane. Store results in sheet
        vids_coplanar = [] # list to keep sheets that we ignore because they support parallel planes (separation angle=180deg)
        vid_map = {} # mapping we need to keep because we are skipping some nodes in the graph and indexing is based on position in list 
        nid_cnt = 0
        for i,v in enumerate(self.g.vs):
            ma_idx = v['ma_idx']
            if v['ma_theta_mean'] > math.radians(self.p_angle_converge):
                vids_coplanar.append(v.index)
                continue
            if len(ma_idx) < 30: # ignore small sheets, eg. often caused by oversegmentation, should actually be merged
                vids_coplanar.append(v.index)
                print('skipping small sheet...')
                continue
            # elif ma.D['ma_radii'][ma_idx].min() > 0.8: # ignore sheets that are not fully going into the building edges
            #     vids_coplanar.append(v.index)
            #     continue
            
            # obtain set of surface points (regardless of what is in/out)
            s_idx = np.mod(ma_idx, ma.m)

            # collect spokes
            c = self.ma.D['ma_coords'][ma_idx]
            f1 = c-self.ma.D['coords'][s_idx]
            f2 = c-self.ma.D['coords'][self.ma.D['ma_qidx'][ma_idx]]
            
            f1_norm = np.linalg.norm(f1, axis=1)[:,None]
            f2_norm = np.linalg.norm(f2, axis=1)[:,None]
            spokes = np.concatenate([f1/f1_norm, f2/f2_norm])

            idx = np.concatenate([s_idx, self.ma.D['ma_qidx'][ma_idx]])

            km = KMeans(n_clusters = 2)
            labels = km.fit_predict(spokes)

            locked = self.ma.D['ma_radii'][ma_idx].min() > 0.8

            self.m.add_node(
                {'s_idx':idx[labels==0], 'spoke_cluster_center':km.cluster_centers_[0], 'locked':locked, 'vid':i}, 
                {'s_idx':idx[labels==1], 'spoke_cluster_center':km.cluster_centers_[1], 'locked':locked, 'vid':i}
            )
            vid_map[i]=nid_cnt
            nid_cnt += 1
