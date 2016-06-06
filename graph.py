from pointio import io_npy
import numpy as np
import math

class Graph:
	def __init__(self, edges=set(), nodes=set()):
		self.edges = edges
		self.nodes = nodes

	def addEdge(self, edge):
		self.edges.add(edge)
	
	def addNode(self, node):
		self.nodes.add(node)

	def __str__(self):
		return str(self.nodes)

class Node:
	def __init__(self, segment_id):
		self.segment_id = segment_id
		self.incident_edges = []
		
		self.avg_bisector = None

	def iterate_neighbours(self):
		for e in self.incident_edges:
			yield e, e.get_neighbour_node(self)

	def __str__(self):
		return "node <{}>".format(self.segment_id)

class Edge:
	def __init__(self, start, end, count):
		self.start = start
		self.end = end
		
		self.count = count

	def get_neighbour_node(self, node):
		if node is self.start:
			return self.end
		else:
			return self.start
			

			
def build_graph(datadict, bisec_avg_dict=None):
	ma_segment = datadict['ma_segment']
	# datadict['seg_link_flip']
	seg_link_adj = datadict['seg_link_adj']

	# build graph datastructure
	node_dict = {}
	edge_list = []
	for start_id, end_id, count in seg_link_adj:

		if start_id not in node_dict:
			n = Node(start_id)
			node_dict[start_id] = n
			# n.avg_bisector = bisec_avg_dict[start_id][1]
		if end_id not in node_dict:
			n = Node(end_id)
			node_dict[end_id] = n
			# n.avg_bisector = bisec_avg_dict[end_id][1]

		start = node_dict[start_id]
		end = node_dict[end_id]

		edge = Edge(start, end, count)
		edge_list.append(edge)
		
		start.incident_edges.append(edge)
		end.incident_edges.append(edge)
		
	g = Graph(set(edge_list), set(node_dict.values()))
	# import ipdb; ipdb.set_trace()
	
	return g

def get_graphs(datadict, min_count=4, avg_bisector_threshold=math.cos(math.radians(15))):
	
	# bisec_avg_dict = compute_segment_aggregate(datadict, key_to_aggregate='ma_bisec')
	
	# build graph
	g = build_graph(datadict)		
	# Find connected components using standard graph traversal algorithm
	graph_list = traverse_repeat(g.nodes, f_conditional=f_min_count, f_arg=min_count)
	
	# for g in graph_list:
	# 	traverse_repeat(g.nodes, f_conditional=f_avg_bisector_threshold, f_arg=avg_bisector_threshold) 
		
	return graph_list
	

def f_min_count(e, arg):
	min_count = arg
	return e.count >= min_count

# def f_avg_bisector_threshold(e, arg):
# 	threshold = arg
# 	return np.dot(e.start.avg_bisector, e.end.avg_bisector) <= threshold

def traverse_repeat(node_set, f_conditional=None, f_arg=None):
	# find connected components in node_set, grow components based on function that is evaluated for its incident edges
	graph_list = []
	while len(node_set) != 0:
		Q = [node_set.pop()]
		V = set()
		E = set()
		while len(Q) != 0:
			node = Q.pop(0)
			V.add(node)

			for e in node.incident_edges:
				# import ipdb; ipdb.set_trace()
				if f_conditional(e, f_arg):
					E.add(e)
					adjacent_node = e.get_neighbour_node(node)
					if (not adjacent_node in V) and (not adjacent_node in Q):
						Q.append(adjacent_node)

		node_set -= set(V)
		g = Graph()
		# import ipdb; ipdb.set_trace()
		g.edges = set(E)
		g.nodes = set(V)
		graph_list.append(g)
	return graph_list


def merge_nodes(datadict, nodes):
	pass
	