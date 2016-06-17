import math, sys
import numpy as np
from time import time
from pointio import io_npy
from ma_util import MAHelper
import argparse
import igraph
from pykdtree.kdtree import KDTree

# INFILE = 'data/scan_npy'
INFILE = "/Users/ravi/git/masbcpp/rdam_blokken_npy"
# INFILE = "/Volumes/Data/Data/pointcloud/AHN2_matahn_samples/ringdijk_opmeer_npy"
# INFILE = "/Volumes/Data/Data/pointcloud/AHN2_matahn_samples/denhaag_a12_npy"

def get_neighbours_ma(data, k=15):
    kdt_ma = KDTree(data)
    return kdt_ma.query(data, k)

class RegionGrower(object):
	"""Segmentation based on region growing. Segment '0' is reserved for unsegmented points. Note that interior & exterior MAT points are concatenated and not treated separately."""

	def __init__(self, mah, bisec_thres=5.0, k=10, only_interior=False, method = 'bisec'):
		if method == 'bisec':
			self.valid_candidate = self.valid_candidate_bisec # or 'normal'
		else:
			self.valid_candidate = self.valid_candidate_normal

		self.p_bisecthres = math.cos((bisec_thres / 180.0) * math.pi)
		self.p_normalthres = math.cos((5.0 / 180.0) * math.pi)
		self.p_thetathres_1 = (50.0 / 180.0) * math.pi # during bisect growing
		self.p_thetathres_2 = (3.0 / 180.0) * math.pi # during theta growing
		self.p_k = k
		self.p_mincount = 10

		self.mah = mah
		# self.filt = self.mah.D['ma_radii'] < 190.

		if only_interior:
			self.ma_coords = self.mah.D['ma_coords_in']
			# self.mah.D['m']
			self.m = self.mah.m
			self.ma_bisec = self.mah.D['ma_bisec_in']
			self.ma_theta = self.mah.D['ma_theta_in']
		else:
			self.ma_coords = self.mah.D['ma_coords']
			# self.mah.D['m']
			self.m = self.mah.m*2
			self.ma_bisec = self.mah.D['ma_bisec']
			self.ma_theta = self.mah.D['ma_theta']

		self.neighbours_dist, self.neighbours_idx = get_neighbours_ma(self.ma_coords, self.p_k)
		# self.estimate_normals()

		self.ma_segment = np.zeros(self.m, dtype=np.int64)
		
		self.region_nr = 1
		self.overwrite_regions = False

	# def estimate_normals(self):
	#	from sklearn.decomposition import PCA
	# 	def compute_normal(neighbours):
	# 		pca = PCA(n_components=3)
	# 		pca.fit(neighbours)
	# 		plane_normal = pca.components_[-1] # this is a normalized normal
	# 		# # make all normals point upwards:
	# 		# if plane_normal[-1] < 0:
	# 		# 	plane_normal *= -1
	# 		return plane_normal

	# 	neighbours = self.ma_coords[self.neighbours_idx]
	# 	t1 = time()
	# 	self.ma_normals = np.empty((self.m,3), dtype=np.float32)
	# 	for i, neighbourlist in enumerate(neighbours):
	# 		self.ma_normals[i] = compute_normal(neighbourlist)
	# 	t2 = time()
	# 	print "finished normal computation in {} s".format(t2-t1)

	def apply_region_growing_algorithm(self, seedpoints):
		"""pop seedpoints and try to grow regions until no more seedpoints are left"""
		
		while len(seedpoints) > 0:
			seed = seedpoints.pop()
			if self.ma_segment[seed] == 0:
				self.grow_region(seed)
				self.region_nr += 1

	def grow_region(self, initial_seed):
		"""Use initial_seed to grow a region by testing if its neighbours are valid candidates. Valid candidates are added to the current region/segment and _its_ neighbours are also tested. Stop when we run out of valid candidates."""
		candidate_stack = [initial_seed]
		point_count = 1
		self.ma_segment[initial_seed] = self.region_nr
		while len(candidate_stack) > 0:
			seed = candidate_stack.pop()
			for neighbour in self.neighbours_idx[seed][1:]:
				if not self.overwrite_regions:
					if self.ma_segment[neighbour] != 0:
						continue
				if self.valid_candidate(seed, neighbour):
					self.ma_segment[neighbour] = self.region_nr
					candidate_stack.append(neighbour)
					point_count += 1
		print("found region nr %d with %d points" % (self.region_nr, point_count))
		return point_count

	def valid_candidate_normal(self, seed, candidate):
		"""candidate is valid if angle between normals of seed and candidate is below preset threshold"""
		if math.fabs(np.dot(self.ma_normals[seed], self.ma_normals[candidate])) > self.p_normalthres:
			return True
		else:
			return False
	
	def valid_candidate_bisec(self, seed, candidate):
		"""candidate is valid if angle between bisectors of seed and candidate is below preset threshold"""
		return np.dot(self.ma_bisec[seed], self.ma_bisec[candidate]) > self.p_bisecthres
		# if np.dot(self.ma_bisec[seed], self.ma_bisec[candidate]) > self.p_bisecthres and math.fabs(self.ma_theta[seed]-self.ma_theta[candidate]) < self.p_thetathres_1:

	def valid_candidate_theta(self, seed, candidate):
		"""candidate is valid if difference between separation angles of seed and candidate is below preset threshold"""
		if math.fabs(self.ma_theta[seed]-self.ma_theta[candidate]) < self.p_thetathres_2:
			return True
		else:
			return False

	def unmark_small_clusters(self):
		"""find all segments that are too small and set their segment to 0"""
		# find cluster ids and sizes
		region_numbers, region_counts = np.unique(self.ma_segment, return_counts=True)

		# find small cluster ids
		to_unmark = region_numbers[region_counts < self.p_mincount]
		# find corresponding indices in ma_segment
		to_unmark = np.in1d(self.ma_segment, to_unmark)
		# set those points as unsegmented
		self.ma_segment[to_unmark] = 0


	def assign_unsegmented_points(self):
		"""experimental stuff to distribute unsegmented points among segments."""
		points = np.where(self.ma_segment==0)[0]

		for p in points:
			neighbours = self.neighbours_idx[p][1:]

			# neighbours = neighbours.where(self.ma_segment != 0)
			neighbour_vecs = self.ma_coords[neighbours] - self.ma_coords[p]
			neighbour_vecs = neighbour_vecs/np.linalg.norm(neighbour_vecs, axis=1)[:,None]
			angles = np.arccos(np.sum(self.ma_bisec[p]*neighbour_vecs,axis=1))
			print(self.ma_segment[neighbours])
			print(angles/math.pi * 180)
			print(min(angles/math.pi * 180))
			# import ipdb; ipdb.set_trace()

def perform_segmentation_bisec(mah, bisec_thres, k, infile=INFILE, **args):
	# find segments based on similiraty in bisector orientation
	R = RegionGrower(mah, bisec_thres=bisec_thres, k=k, method='bisec', **args)
	seedpoints = list( np.random.permutation(R.m) )
	R.apply_region_growing_algorithm(seedpoints)
	R.unmark_small_clusters()
	# print(np.unique(R.ma_segment, return_counts=True))
	
	# now try to find segments that have a large separation angle (and unstable bisector orientation)
	seedpoints = list(np.where(np.logical_and(R.ma_segment==0, R.ma_theta > (175.0/180)*math.pi ))[0])
	# R.overwrite_regions = True
	R.valid_candidate = R.valid_candidate_theta
	R.apply_region_growing_algorithm(seedpoints)
	R.unmark_small_clusters()
	# R.assign_unsegmented_points()
	# print(np.unique(R.ma_segment, return_counts=True))
	
	# build graph
	g = igraph.Graph(directed=False)
	
	ma_segment_dict = {}
	for i, seg_id in enumerate(R.ma_segment):
		if seg_id in ma_segment_dict:
			ma_segment_dict[seg_id].append(i)
		else:
			ma_segment_dict[seg_id]=[]
		
	for k,v in ma_segment_dict.items():
		g.add_vertex(ma_idx=v)

	# import ipdb; ipdb.set_trace()
	
	
	ma_segment = np.zeros(R.mah.m*2, dtype=np.int64)
	graph2segmentlist(g, ma_segment)
	
	find_relations(mah, infile)
	
	for start_id, end_id, count in mah.D['seg_link_adj']:
		g.add_edge(start_id, end_id, adj_count=count)

	g.write_pickle(infile+'/ma_segment.pickle')

	D['ma_segment'] = ma_segment
	io_npy.write_npy(infile, D, ['ma_segment'])
	
	return g
	
def graph2segmentlist(g, ma_segment):
	for v in g.vs():
		ma_segment[ v['ma_idx'] ] = v.index	

# def perform_segmentation_normal(mah):	
# 	R = RegionGrower(mah, method='normal')
# 	seedpoints = list( np.random.permutation(R.m) )
# 	R.apply_region_growing_algorithm(seedpoints)
# 	R.unmark_small_clusters()
# 	# import ipdb; ipdb.set_trace()
# 	seedpoints = list(np.where(np.logical_and(R.ma_segment==0, R.ma_theta > (175.0/180)*math.pi ))[0])

# 	print R.region_counts
# 	D['ma_segment'] = R.ma_segment
# 	io_npy.write_npy(INFILE, D, ['ma_segment'])

def find_relations(ma, infile=INFILE, only_interior=False):
	"""
	Find topological relations between segments. Output for each relation: 
		(segment_1, segment_2, count)
	the higher count the stronger the relation.
	"""

	def find_flip_relations():
		"""Find for each pair of segments how many times they are connected by a shared feature point.
			In a pair of segments (tuple) the lowest segmend_id is always put first
		"""

		pdict = {}
		for i in np.arange(ma.m):
			coord_id = i# % ma.m
			s_in = ma.D['ma_segment'][i]
			s_out = ma.D['ma_segment'][i + ma.m]

			if not (s_in == 0 or s_out == 0):
				if s_in < s_out:
					pair = s_in, s_out
				else:
					pair = s_out, s_in

				if pair in pdict:
					pdict[pair]+= 1
				else:
					pdict[pair] = 1

		return pdict


	def find_adjacency_relations():
		"""find pairs of adjacent segments
		"""
		if only_interior:
			neighbours_dist, neighbours_idx = get_neighbours_ma(ma.D['ma_coords_in'])
			m=ma.m
		else:
			neighbours_dist, neighbours_idx = get_neighbours_ma(ma.D['ma_coords'])
			m=ma.m*2
		pdict = {}

		for i in np.arange(m):
			seg_id = ma.D['ma_segment'][i]

			neighbours = neighbours_idx[i][1:]
			n_seg = ma.D['ma_segment'][neighbours]

			for n_seg_id in n_seg:

				if not (seg_id == n_seg_id) :
					if seg_id < n_seg_id:
						pair = seg_id, n_seg_id
					else:
						pair = n_seg_id, seg_id

					if pair[0] == 0: continue

					if pair in pdict:
						pdict[pair]+= 1
					else:
						pdict[pair] = 1

		# import ipdb; ipdb.set_trace()
		return pdict


	if not only_interior:
		flip_relations = find_flip_relations()
		ma.D['seg_link_flip'] = np.zeros(len(flip_relations), dtype = "3int32")
		i=0
		for (s, e), cnt in flip_relations.items():
			ma.D['seg_link_flip'][i] = [s,e,cnt]
			i+=1
		io_npy.write_npy(infile, ma.D, ['seg_link_flip'])

	adj_relations = find_adjacency_relations()
	ma.D['seg_link_adj'] = np.zeros(len(adj_relations), dtype = "3int32")
	i=0
	for (s, e), cnt in adj_relations.items():
		ma.D['seg_link_adj'][i] = [s,e,cnt]
		i+=1
	io_npy.write_npy(infile, ma.D, ['seg_link_adj'])

def compute_segment_aggregate(datadict, key_to_aggregate='ma_coords'):
    """Compute avarage coordinate for each segment"""

    segment_dict = {}
    # segment_point_sums = {}
    for i, segment in enumerate(datadict['ma_segment']):
        # slicing is not copying!
        if segment in segment_dict:
            segment_dict[segment][0] += 1
            segment_dict[segment][1] += datadict[key_to_aggregate][i]
        else:
            segment_dict[segment] = [1, np.copy(datadict[key_to_aggregate][i])]

    for key, value in segment_dict.items():
        segment_dict[key][1] = value[1]/value[0]

    return segment_dict

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Region growing for MAT sheets')
	parser.add_argument('infile', help='npy file', default=INFILE)
	parser.add_argument('-k', help='Number of neighbours to use', default=10, type=int)
	parser.add_argument('-b', help='Bisec threshold in degrees', default=5.0, type=float)
	# parser.add_argument('--topo', help='Also compute topological links', dest='topo', action='store_true')
	parser.add_argument('--interior', help='Only use interior MAT', dest='interior', action='store_true')
	args = parser.parse_args()

	D = io_npy.read_npy(args.infile)
	mah = MAHelper(D)
	g = perform_segmentation_bisec(mah, bisec_thres=args.b, k=args.k, infile=args.infile, only_interior=args.interior)
	# if args.topo:
		# find_relations(mah, infile=args.infile, only_interior=args.interior)
		

