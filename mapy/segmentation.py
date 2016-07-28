import math, sys
import numpy as np
from time import time
import igraph
from pykdtree.kdtree import KDTree

from mapy.io import npy 

# INFILE = 'data/scan_npy'
INFILE = "/Users/ravi/git/masbcpp/rdam_blokken_npy"
# INFILE = "/Volumes/Data/Data/pointcloud/AHN2_matahn_samples/ringdijk_opmeer_npy"
# INFILE = "/Volumes/Data/Data/pointcloud/AHN2_matahn_samples/denhaag_a12_npy"

# Print iterations progress
def printProgress (iteration, total, prefix = '', suffix = '', decimals = 2, barLength = 100):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : number of decimals in percent complete (Int) 
        barLength   - Optional  : character length of bar (Int) 
    """
    if total==0: return
    filledLength    = int(round(barLength * iteration / float(total)))
    percents        = '{:>6.2f}'.format(round(100.00 * (iteration / float(total)), decimals))
    bar             = '#' * filledLength + '-' * (barLength - filledLength)
    sys.stdout.write('%s [%s] %s%s %s\r' % (prefix, bar, percents, '%', suffix)),
    sys.stdout.flush()
    if iteration == total:
        print("\n")

def get_neighbours_ma(data, k=15):
    kdt_ma = KDTree(data)
    return kdt_ma.query(data, k)

class RegionGrower(object):
	"""Segmentation based on region growing. Segment '0' is reserved for unsegmented points. Note that interior & exterior MAT points are concatenated and not treated separately."""

	def __init__(self, mah, **kwargs):
		self.p = {
			'bisec_thres':10.0,
			'theta_thres':10.0,
			'k':10,
			'only_interior':False,
			'method':'bisec',
			'mincount':10
		}
		self.p.update(kwargs)
		# import ipdb;ipdb.set_trace()

		if self.p['method'] == 'bisec':
			self.valid_candidate = self.valid_candidate_bisectheta # or 'normal'
		else:
			self.valid_candidate = self.valid_candidate_normal

		self.p_bisecthres = math.cos((self.p['bisec_thres'] / 180.0) * math.pi)
		self.p_normalthres = math.cos((5.0 / 180.0) * math.pi)
		self.p_thetathres_1 = (self.p['theta_thres'] / 180.0) * math.pi # during bisect growing
		self.p_thetathres_2 = (self.p['theta_thres'] / 180.0) * math.pi # during theta growing
		self.p_k = self.p['k']
		self.p_mincount = self.p['mincount']

		self.mah = mah
		# self.filt = self.mah.D['ma_radii'] < 190.

		if self.p['only_interior']:
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
		totalcount = len(seedpoints)
		pointcount = 0
		seedpoints = set(seedpoints)
		printProgress(0, totalcount, prefix='Segmentation progress from {} seeds:'.format(totalcount), barLength=20)
		while len(seedpoints) > 0:
			seed = seedpoints.pop()
			seedpoints -= self.grow_region(seed)
			printProgress(totalcount-len(seedpoints), totalcount, prefix='Segmentation progress from {} seeds:'.format(totalcount), suffix='({} regions)'.format(self.region_nr), barLength=20)
			self.region_nr += 1	

	def grow_region(self, initial_seed):
		"""Use initial_seed to grow a region by testing if its neighbours are valid candidates. Valid candidates are added to the current region/segment and _its_ neighbours are also tested. Stop when we run out of valid candidates."""
		candidate_stack = [initial_seed]
		point_count = 1
		self.ma_segment[initial_seed] = self.region_nr
		neighbours_in_region = set() 
		while len(candidate_stack) > 0:
			seed = candidate_stack.pop()
			for neighbour in self.neighbours_idx[seed][1:]:
				if not self.overwrite_regions:
					if self.ma_segment[neighbour] != 0:
						continue
				if self.valid_candidate(seed, neighbour):
					self.ma_segment[neighbour] = self.region_nr
					candidate_stack.append(neighbour)
					neighbours_in_region.add(neighbour)
					point_count += 1
		# print("found region nr %d with %d points" % (self.region_nr, point_count))
		return neighbours_in_region

	def valid_candidate_normal(self, seed, candidate):
		"""candidate is valid if angle between normals of seed and candidate is below preset threshold"""
		if math.fabs(np.dot(self.ma_normals[seed], self.ma_normals[candidate])) > self.p_normalthres:
			return True
		else:
			return False
	
	def valid_candidate_bisectheta(self, seed, candidate):
		"""candidate is valid if angle between bisectors of seed and candidate is below preset threshold"""
		return self.valid_candidate_bisec(seed, candidate) and self.valid_candidate_theta(seed, candidate)

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

def perform_segmentation_bisec(mah, **kwargs):
	# find segments based on similiraty in bisector orientation
	print("Initiating region grower...")
	R = RegionGrower(mah, **kwargs)
	# seedorder = np.argsort(mah.D['ma_radii'])[::-1].tolist() # reverse
	seedorder = list( np.random.permutation(R.m) )
	print("\nPerforming bisector-based region growing...")
	R.apply_region_growing_algorithm(seedorder)
	R.unmark_small_clusters()
	# print(np.unique(R.ma_segment, return_counts=True))
	
	# now try to find segments that have a large separation angle (and unstable bisector orientation)
	seedpoints = list(np.where(np.logical_and(R.ma_segment==0, R.ma_theta > (175.0/180)*math.pi ))[0])
	# R.overwrite_regions = True
	R.valid_candidate = R.valid_candidate_theta
	print("Performing theta-based region growing...")
	R.apply_region_growing_algorithm(seedpoints)
	R.unmark_small_clusters()
	# R.assign_unsegmented_points()
	# print(np.unique(R.ma_segment, return_counts=True))
	
	# build graph
	print("Constructing graph...")
	g = igraph.Graph(directed=False)
	
	ma_segment_dict = {}
	for i, seg_id in enumerate(R.ma_segment):
		if seg_id in ma_segment_dict:
			ma_segment_dict[seg_id].append(i)
		else:
			ma_segment_dict[seg_id]=[]
	
	for k,v in ma_segment_dict.items():
		g.add_vertex(ma_idx=v)
	
	mah.D['ma_segment'] = np.zeros(R.mah.m*2, dtype=np.int64)
	graph2segmentlist(g, mah.D['ma_segment'])
	
	find_relations(mah, kwargs['infile'])
	
	for start_id, end_id, count in mah.D['seg_link_adj']:
		g.add_edge(start_id, end_id, adj_count=count)

	compute_segment_aggregate(g, mah.D, 'ma_coords')
	compute_segment_aggregate(g, mah.D, 'ma_bisec')
	compute_segment_aggregate(g, mah.D, 'ma_theta')
	mah.D['ma_segment_graph'] = g
	npy.write(kwargs['infile'], mah.D, ['ma_segment', 'ma_segment_graph'])
	
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
# 	npy.write(INFILE, D, ['ma_segment'])

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


	def find_adjacency_relations(k=25):
		"""find pairs of adjacent segments
		"""
		if only_interior:
			neighbours_dist, neighbours_idx = get_neighbours_ma(ma.D['ma_coords_in'], k=k)
			m=ma.m
		else:
			neighbours_dist, neighbours_idx = get_neighbours_ma(ma.D['ma_coords'], k=k)
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

					# don't add edges to unsegmented points (stored in segment 0)
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
		npy.write(infile, ma.D, ['seg_link_flip'])

	adj_relations = find_adjacency_relations()
	ma.D['seg_link_adj'] = np.zeros(len(adj_relations), dtype = "3int32")
	i=0
	for (s, e), cnt in adj_relations.items():
		ma.D['seg_link_adj'][i] = [s,e,cnt]
		i+=1
	npy.write(infile, ma.D, ['seg_link_adj'])

def compute_segment_aggregate(g, datadict, key='ma_coords'):
	"""Compute eg. avarage coordinate for each segment"""
	for v in g.vs:
		v[key+'_mean'] = datadict[key][v['ma_idx']].mean(axis=0)