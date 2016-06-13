from ma_util import MAHelper
from pointio import io_npy 
import igraph

def get_graphs(datadict, min_count):
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

        if start_id not in vertex_dict:
            vertex_dict[start_id] = vertex_cnt
            g.add_vertex(segment_id=start_id)
            vertex_cnt += 1
        if end_id not in vertex_dict:
            vertex_dict[end_id] = vertex_cnt
            g.add_vertex(segment_id=end_id)
            vertex_cnt += 1
        if count >= min_count:
            g.add_edge(vertex_dict[start_id], vertex_dict[end_id], adj_count=count)
    
    #g.delete_edges(g.es.select(adj_count_lt=min_count)
    
    return g.clusters(igraph.WEAK).subgraphs()
    
if __name__ == '__main__':
    D = io_npy.read_npy("/Users/ravi/git/masbcpp/rdam_blokken_npy")
    mah = MAHelper(D)
    g = get_graphs(D, 6)