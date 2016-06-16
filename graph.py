from ma_util import MAHelper
from pointio import io_npy 
import igraph
import math
import numpy as np
from itertools import chain

def construct_graph(datadict, min_count):
    """return connected components in adjacency graph with an adjacency count of at least min_count"""
    ma_segment = datadict['ma_segment']
    # datadict['seg_link_flip']
    seg_link_adj = datadict['seg_link_adj']

    # build graph datastructure
    g = igraph.Graph(directed=False)
    # igraph works with vertex ids that are depend on their order in the graph so we need to map the segment_ids
    vertex_dict = {}
    vertex_cnt = 0
    for start_id, end_id, count in seg_link_adj:
        
        for index in [start_id, end_id]:
            if index not in vertex_dict:
                vertex_dict[index] = vertex_cnt
                g.add_vertex(segment_id=index, ma_idx=np.argwhere(ma_segment==index)[:,0].tolist() )
                vertex_cnt += 1
        if count >= min_count:
            g.add_edge(vertex_dict[start_id], vertex_dict[end_id], adj_count=count)
    
    #g.delete_edges(g.es.select(adj_count_lt=min_count)
    return g
    
    
def get_graphs(datadict, min_count):
    g = construct_graph(datadict, min_count)
    # g.delete_edges(g.es.select(adj_count_lt=min_count)
    return g.clusters(igraph.WEAK).subgraphs()
    
def assign_from_aggregate_dict(g, dic):
    # segment_ids = g.vs['segment_id']
    for v in g.vs:
        v['count'], v['bisec_avg'] = dic[v['segment_id']]
        
def contract_edges(g, threshold=math.radians(5)):
    
    for e in g.es:
        s = g.vs[e.source]
        t = g.vs[e.target]
        e['avg_bisec_angle'] = np.arccos(np.dot(s['bisec_avg'], t['bisec_avg']))
        
    tg = g.copy()
    tg.delete_edges(tg.es.select(avg_bisec_angle_gt=threshold))
    
    g.contract_vertices(tg.clusters(igraph.WEAK).membership, 
        combine_attrs={ 'segment_id': lambda x: x, 
                        'count':'sum', 
                        'ma_idx': lambda x: [e for e in chain(*x)]
        })
    g.simplify(combine_edges='sum')

if __name__ == '__main__':
    from region_growing import compute_segment_aggregate
    
    D = io_npy.read_npy("/Users/ravi/git/masbcpp/rdam_blokken_npy")
    g = construct_graph(D, 6)
    
    bisec_aggregate = compute_segment_aggregate(D, 'ma_bisec')
    assign_from_aggregate_dict(g, bisec_aggregate)
    contract_edges(g)
    # update_points(g, D['ma_segment'])
    import ipdb; ipdb.set_trace()
    # mah = MAHelper(D)
    