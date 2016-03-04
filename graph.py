from pointio import io_npy
import numpy as np

class Graph:
	def __init__(self):
		self.edge_list = set()
		self.node_list = set()

	def addEdge(self, edge):
		self.edge_list.add(edge)
	
	def addNode(self, node):
		self.node_list.add(node)

class Node:
	def __init__(self, segment_id):
		self.segment_id = segment_id
		self.incident_edges = []

class Edge:
	def __init__(self, start, end, count):
		self.start = start
		self.end = end
		self.count = count

infile = "/Users/ravi/git/masbcpp/rdam_blokken_npy"
datadict = io_npy.read_npy(infile)

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
	