import math
import numpy as np
from time import time
from pykdtree.kdtree import KDTree
from sklearn.decomposition import PCA
from pointio import io_npy
from ma_util import MAHelper

# INFILE = 'data/scan_npy'
INFILE = "/Users/ravi/git/masbcpp/rdam_blokken_npy"


class RegionGrower(object):

	def __init__(self, datadict, method = 'bisec'):
		# self.kdt_surf = KDTree(sel f.D['coords'])
		if method == 'bisec':
			self.valid_candidate = self.valid_candidate_bisec # or 'normal'
		else:
			self.valid_candidate = self.valid_candidate_normal

		self.p_bisecthres = math.cos((10.0 / 180.0) * math.pi)
		self.p_normalthres = math.cos((5.0 / 180.0) * math.pi)
		self.p_thetathres_1 = (50.0 / 180.0) * math.pi # during bisect growing
		self.p_thetathres_2 = (3.0 / 180.0) * math.pi # during theta growing
		self.p_k = 15
		self.p_mincount = 50

		self.mah = MAHelper(datadict)
		self.filt = self.mah.D['ma_radii'] < 190.
		# import ipdb; ipdb.set_trace()

		self.ma_coords = self.mah.D['ma_coords'][self.filt]
		# self.mah.D['m']
		self.m = self.ma_coords.shape[0]
		self.ma_bisec = self.mah.D['ma_bisec'][self.filt]
		self.ma_theta = self.mah.D['ma_theta'][self.filt]

		# elif method == 'normal':
		# self.find_neighbours()
		self.neighbours_dist, self.neighbours_idx = self.mah.get_neighbours_ma(self.filt, self.p_k)
		# self.estimate_normals()

		self.ma_segment = np.zeros(self.m, dtype=np.int64)
		self.region_nr = 1
		self.region_counts = [0]
		self.overwrite_regions = False

	# def find_neighbours(self):
	# 	self.kdt_ma = KDTree(self.ma_coords)
	# 	self.neighbours_dist, self.neighbours_idx = self.kdt_ma.query(
	# 		self.ma_coords, 
	# 		self.p_k 
	# 	)


	def estimate_normals(self):
		def compute_normal(neighbours):
			pca = PCA(n_components=3)
			pca.fit(neighbours)
			plane_normal = pca.components_[-1] # this is a normalized normal
			# # make all normals point upwards:
			# if plane_normal[-1] < 0:
			# 	plane_normal *= -1
			return plane_normal

		neighbours = self.ma_coords[self.neighbours_idx]
		# p = Pool()
		t1 = time()
		self.ma_normals = np.empty((self.m,3), dtype=np.float32)
		for i, neighbourlist in enumerate(neighbours):
			self.ma_normals[i] = compute_normal(neighbourlist)
		t2 = time()
		print "finished normal computation in {} s".format(t2-t1)
		# p.close()

	def apply_region_growing_algorithm(self, seedpoints):
		#pick seedpoints at random for now
		
		while len(seedpoints) > 0:
			seed = seedpoints.pop()
			if self.ma_segment[seed] == 0:
				self.region_counts.append(self.grow_region(seed))
				self.region_nr += 1
		# print region_counts
		# print("found %d regions" % region_nr)

	def grow_region(self, initial_seed):
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
		print "found region with %d points" % point_count
		return point_count

	def valid_candidate_normal(self, seed, candidate):
		if math.fabs(np.dot(self.ma_normals[seed], self.ma_normals[candidate])) > self.p_normalthres:
			return True
		else:
			return False
	
	def valid_candidate_bisec(self, seed, candidate):
		if np.dot(self.ma_bisec[seed], self.ma_bisec[candidate]) > self.p_bisecthres:
		# if np.dot(self.ma_bisec[seed], self.ma_bisec[candidate]) > self.p_bisecthres and math.fabs(self.ma_theta[seed]-self.ma_theta[candidate]) < self.p_thetathres_1:
			return True
		else:
			return False

	def valid_candidate_theta(self, seed, candidate):
		if math.fabs(self.ma_theta[seed]-self.ma_theta[candidate]) < self.p_thetathres_2:
			return True
		else:
			return False


	def unmark_bad_regions(self):
		new_counts = [self.region_counts[0]]
		for i, c in enumerate(self.region_counts):
			if c < self.p_mincount:
				self.ma_segment[self.ma_segment==i] = 0
				new_counts[0] += c
			else:
				new_counts.append(c)
		self.region_counts = new_counts
		self.region_nr = self.ma_segment.max()+1

	# def realign_segments(self):
	# 	self.region_numbers, self.region_counts = np.unique(self.ma_segment, return_counts=True)
	# 	mp = np.arange(len(self.region_numbers))

	# 	replace = { old:new for old, new in zip(self.region_numbers, mp)}
		
	# 	mp[replace.keys()] = replace.values()
	# 	self.ma_segment = mp[self.ma_segment]


	def assign_unsegmented_points(self):
		points = np.where(self.ma_segment==0)[0]

		for p in points:
			neighbours = self.neighbours_idx[p][1:]

			# neighbours = neighbours.where(self.ma_segment != 0)
			neighbour_vecs = self.ma_coords[neighbours] - self.ma_coords[p]
			neighbour_vecs = neighbour_vecs/np.linalg.norm(neighbour_vecs, axis=1)[:,None]
			angles = np.arccos(np.sum(self.ma_bisec[p]*neighbour_vecs,axis=1))
			print self.ma_segment[neighbours]
			print angles/math.pi * 180
			print min(angles/math.pi * 180)
			# import ipdb; ipdb.set_trace()

def do_bisec_based():
	D = io_npy.read_npy(INFILE)
	
	R = RegionGrower(D, method='bisec')
	seedpoints = list( np.random.permutation(R.m) )
	R.apply_region_growing_algorithm(seedpoints)
	R.unmark_bad_regions()
	# import ipdb; ipdb.set_trace()
	

	if 1: # do theta 180deg planes
		seedpoints = list(np.where(np.logical_and(R.ma_segment==0, R.ma_theta > (175.0/180)*math.pi ))[0])
		# R.overwrite_regions = True
		R.valid_candidate = R.valid_candidate_theta
		# import ipdb; ipdb.set_trace()
		R.apply_region_growing_algorithm(seedpoints)
		R.region_counts[0] = np.count_nonzero(R.ma_segment == 0)
		# R.unmark_bad_regions()
		# R.assign_unsegmented_points()
		print R.region_counts

	ma_segment = np.zeros(R.mah.m*2, dtype=np.int64)
	ma_segment[R.filt] = R.ma_segment

	D['ma_segment'] = ma_segment
	io_npy.write_npy(INFILE, D, ['ma_segment'])

def do_normal_based():
	D = io_npy.read_npy(INFILE)
	
	R = RegionGrower(D, method='normal')
	seedpoints = list( np.random.permutation(R.m) )
	R.apply_region_growing_algorithm(seedpoints)
	R.unmark_bad_regions()
	# import ipdb; ipdb.set_trace()
	seedpoints = list(np.where(np.logical_and(R.ma_segment==0, R.ma_theta > (175.0/180)*math.pi ))[0])

	print R.region_counts
	D['ma_segment'] = R.ma_segment
	io_npy.write_npy(INFILE, D, ['ma_segment'])

if __name__ == '__main__':
	# do_normal_based()
	do_bisec_based()

#compute k neighbors for each point
#select seed point