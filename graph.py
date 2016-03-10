from pointio import io_npy
import numpy as np

class Graph:
	def __init__(self):
		self.edges = set()
		self.nodes = set()

	def addEdge(self, edge):
		self.edges.add(edge)
	
	def addNode(self, node):
		self.nodes.add(node)

class Node:
	def __init__(self, segment_id):
		self.segment_id = segment_id
		self.incident_edges = []

	def iterate_neighbours(self):
		for e in self.incident_edges:
			yield e, e.get_neighbour_node(self)

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


def get_graphs(datadict, min_count=40):
	ma_segment = datadict['ma_segment']
	# datadict['seg_link_flip']
	seg_link_adj = datadict['seg_link_adj']

	node_dict = {}
	edge_list = []
	for start_id, end_id, count in seg_link_adj:

		if not node_dict.has_key(start_id):
			node_dict[start_id] = Node(start_id)
		if not node_dict.has_key(end_id):
			node_dict[end_id] = Node(end_id)

		start = node_dict[start_id]
		end = node_dict[end_id]

		edge = Edge(start, end, count)
		edge_list.append(edge)
		
		start.incident_edges.append(edge)
		end.incident_edges.append(edge)
		
	# Find connected components
	node_set = set(node_dict.values())
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
				if e.count >= min_count:
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

