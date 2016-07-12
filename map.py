
allowed_kinds = 'intersect', 'match', 'parallel'

class Map(object):

    def __init__(self):
        self.ns = []
        self.es = []
        self.fs = []

    def add_nodes(self, n):
        for i in xrange(n):
            self.add_node()

    def add_node(self, **kwargs):
        hn1 = HalfNode(**kwargs)
        hn2 = HalfNode(**kwargs)
        hn1.twin = hn2
        hn2.twin = hn1
        self.ns.extend([hn1, hn2])
        return hn1, hn2

    def add_edges(self, node_pairs, kind):
        for s,t in node_pairs:
            self.add_edge(s,t, kind=kind)

    def add_edge(self, s_idx, t_idx, kind):
        s, t = self.ns[s_idx], self.ns[t_idx]
        e = Edge(s, t, kind)
        s.edges.setdefault(kind,[]).append(e)
        t.edges.setdefault(kind,[]).append(e)
        self.es.append(e)
        
    def add_face(self, hn_idx):
        f = Face(self.ns[hn_idx])
        self.fs.append(f)

def other(this, pair):
    """assume pair is an indexable iterable of lenth 2 and this is an element in pair, this function will return the other element"""
    return pair[not pair.index(this)]

class HalfNode(object):

    def __init__(self, **kwargs):
        self.twin = None
        self.edges = {} # 
        if 'face' in kwargs:
            self.face = kwargs['face']
        else:
            self.face = None
        if 'pointset' in kwargs:
            self.pointset = kwargs['pointset']
        else:
            self.pointset = None

    def cycle(self, kind):
        yield self
        follow_edge = self.edges[kind][0]
        next = other(self, follow_edge.nodes)
        while next!=self:
            yield next
            if len(next.edges[kind])==1: yield None
            follow_edge = other(follow_edge, next.edges[kind])
            next = other(next, follow_edge.nodes)

    # def __repr__(self):
    #     return "HalfNode edges:{}".format(self.edges)

class Edge(object):
    
    def __init__(self, source, target, kind):
        self.nodes = [source, target]
        self.kind = kind # 'intersect' 'parallel' 'match'

    # def __repr__(self):
    #     return "Edge [{}] nodes:{}".format(self.kind, self.nodes)

class Face(object):

    def __init__(halfnode):
        self.halfnode = halfnode

    # def __repr__(self):
    #     return "Face node:{}".format(self.halfnodes)


def test():
    m = Map()
    m.add_nodes(3)
    m.add_edges([(0,2),(2,4),(4,0)], kind='intersect')
    n = m.ns[0]
    for x in n.cycle(kind='intersect'):
        print x

