from .util import MAHelper
import igraph
import math
import numpy as np
from itertools import chain

# def construct_graph(datadict, min_count):
#     """return connected components in adjacency graph with an adjacency count of at least min_count"""
#     ma_segment = datadict['ma_segment']
#     # datadict['seg_link_flip']
#     seg_link_adj = datadict['seg_link_adj']

#     # build graph datastructure
#     g = igraph.Graph(directed=False)
#     # igraph works with vertex ids that are depend on their order in the graph so we need to map the segment_ids
#     vertex_dict = {}
#     vertex_cnt = 0
#     for start_id, end_id, count in seg_link_adj:
        
#         for index in [start_id, end_id]:
#             if index not in vertex_dict:
#                 vertex_dict[index] = vertex_cnt
#                 g.add_vertex(segment_id=index, ma_idx=np.argwhere(ma_segment==index)[:,0].tolist() )
#                 vertex_cnt += 1
#         if count >= min_count:
#             g.add_edge(vertex_dict[start_id], vertex_dict[end_id], adj_count=count)
    
#     #g.delete_edges(g.es.select(adj_count_lt=min_count)
#     return g
    
    
# def get_graphs(datadict, min_count):
#     g = construct_graph(datadict, min_count)
#     # g.delete_edges(g.es.select(adj_count_lt=min_count)
#     return g.clusters(igraph.WEAK).subgraphs()
def get_graph_library():
    graph_library = {}
    g = igraph.Graph()
    g.add_vertices(4)
    g.add_edges([(0,1), (1,2), (2,3), (3,0), (0,2)])
    graph_library['flatcube_top'] = g

    return graph_library
        
def contract_edges(g, threshold=10):
    threshold = math.radians(threshold)
    for e in g.es:
        s = g.vs[e.source]
        t = g.vs[e.target]
        e['delta_ma_bisec'] = np.arccos(np.dot(s['ma_bisec_mean'], t['ma_bisec_mean']))
        
    tg = g.copy()
    tg.delete_edges(tg.es.select(delta_ma_bisec_gt=threshold))
    
    # note: some attributes may get lost here!
    g.contract_vertices(tg.clusters(igraph.WEAK).membership, 
        combine_attrs={ 'ma_bisec_mean': lambda x: np.mean(x, axis=0), 
                        'ma_theta_mean': lambda x: np.mean(x, axis=0), 
                        'ma_coords_mean': lambda x: np.mean(x, axis=0), 
                        'ma_idx': lambda x: [e for e in chain(*x)]
        })
    g.simplify(combine_edges='sum')

def update_ma_segment(g, datadict):
    for v in g.vs:
        datadict['ma_segment'][v['ma_idx']] = v.index