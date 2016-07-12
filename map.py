

class Map(object):

    def __init__(self):
        self.hnodes = []
        self.edges
        self.faces

    def add_node(self, **kwargs):
        hn1 = HalfNode()
        hn2 = HalfNode()
        hn1.twin = hn2
        hn2.twin = hn1
        hnodes.extend([hn1, hn2])
        return hn1, hn2

    def add_edge(self, s_idx, t_idx, kind):
        s, t = self.hnodes[s_idx], self.hnodes[t_idx]
        e = Edge(s, t, kind)
        s.es[kind] = [e]
        t.es[kind] = [e]
        
    def add_face(self, hn_idx):
        f = Face(self.hnodes[hn_idx])
        self.faces.append(f)

class HalfNode(object):

    def __init__(self):
        self.twin = None
        self.es = {} # 'intersect' 'parallel' 'match'
        self.face = None
        self.pointset = None

    def iterate(self, kind)
        follow_edge = self.es[kind][0]
        next = set(follow_edge.nodes)-self
        while next!=self:
            yield next
            if len(next.es[kind])==1: yield None
            follow_edge = set(next.es[kind])-follow_edge
            next = set(follow_edge.nodes)-next

class Edge(object):
    
    def __init__(self, source, target, kind):
        self.nodes = [source, target]
        self.kind = kind # 'intersect' 'parallel' 'match'

class Face(object):

    def __init__(halfnode):
        self.halfnode = halfnode