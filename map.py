
allowed_kinds = 'intersect', 'match', 'parallel'

class Map(object):

    def __init__(self):
        self.ns = []
        self.es = []
        self.fs = []

    def add_nodes(self, n):
        for i in xrange(n):
            self.add_node()

    def add_node(self, kwargs1={}, kwargs2={}):
        hn1 = HalfNode(**kwargs1)
        hn2 = HalfNode(**kwargs2)
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
        return e
        
    def add_face(self, halfnode, **kwargs):
        f = Face(halfnode, **kwargs)
        self.fs.append(f)
        return f

def other(this, pair):
    """assume pair is an indexable iterable of lenth 2 and this is an element in pair, this function will return the other element"""
    return pair[not pair.index(this)]

class HalfNode(object):

    def __init__(self, face=None, **kwargs):
        self.twin = None
        self.edges = {} # 
        self.face = face
        self.attributes = {}
        self.attributes.update(kwargs)

    def __getitem__(self, key):
        return self.attributes[key]
    
    def __setitem__(self, key, value):
        self.attribute[key] = value

    def cycle(self, kind, direction=0): # direction = [0,1] eg start with first or second edge from this node
        if not kind in self.edges.keys(): return
        follow_edge = self.edges[kind][direction]
        next = other(self, follow_edge.nodes)
        while next!=self:
            yield next
            if len(next.edges[kind])==1: return
            follow_edge = other(follow_edge, next.edges[kind])
            next = other(next, follow_edge.nodes)
        yield self

    # def __repr__(self):
    #     return "HalfNode edges:{}".format(self.edges)

class Edge(object):
    
    def __init__(self, source, target, kind):
        self.nodes = [source, target]
        self.kind = kind # 'intersect' 'parallel' 'match'

    # def twin(self):
    #     return self.nodes[0]

    # def __repr__(self):
    #     return "Edge [{}] nodes:{}".format(self.kind, self.nodes)

class Face(object):

    def __init__(self, halfnode, **kwargs):
        self.halfnode = halfnode
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
        print x

