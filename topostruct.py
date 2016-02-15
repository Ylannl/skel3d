import math, sys
from time import time
import numpy as np
from pointio import io_npy
from ma_util import MAHelper
from povi import App

# INFILE = 'data/scan_npy'
INFILE = "/Users/ravi/git/masbcpp/rdam_blokken_npy"
# INFILE = "/Volumes/Data/Data/pointcloud/AHN2_matahn_samples/ringdijk_opmeer_npy"
INFILE = "/Volumes/Data/Data/pointcloud/AHN2_matahn_samples/denhaag_a12_npy"


# import ipdb; ipdb.set_trace()
def timeit(func):
	t0 = time()
	r = func()
	print("Executed function %s took %f s" % (func.__name__, time()-t0))
	return r

def assign_seg_point():
	"""find for each coord the set of segments it is linked to"""
	pdict = {}

	for i, segment in enumerate(ma.D['ma_segment']):
		coord_id = i# % ma.m
		fp = coord_id
		fq = ma.D['ma_qidx'][coord_id]

		for idx in [fp,fq]:
			if pdict.has_key(idx):
				pdict[idx].add(segment)
			else:
				pdict[idx] = set([segment])
	del pdict[-1]

	datadict['segment_count'] = np.array([ len(s) for s in pdict.values() ], dtype=np.int32)
	io_npy.write_npy(INFILE, datadict, ['segment_count'])


def count_refs():
	"""count the number of times each coord is used as feature point for a medial ball"""

	pdict = {}
	for i in np.arange(ma.m*2):
		coord_id = i# % ma.m
		fp = coord_id %ma.m
		fq = ma.D['ma_qidx'][coord_id]

		for idx in [fp,fq]:
			if pdict.has_key(idx):
				pdict[idx] += 1
			else:
				pdict[idx] = 1
	del pdict[-1]

	return np.array(pdict.values(), dtype=np.int32)

	# datadict['ref_count'] = np.array(pdict.values(), dtype=np.int32)
	# io_npy.write_npy(INFILE, datadict, ['ref_count'])
	# import ipdb; ipdb.set_trace()

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

			if pdict.has_key(pair):
				pdict[pair]+= 1
			else:
				pdict[pair] = 1

	return pdict
	
	# datadict['flip_relations'] = np.array(pdict.values(), dtype=np.int32)
	# io_npy.write_npy(INFILE, datadict, ['ref_count'])
	import ipdb; ipdb.set_trace()


def find_adjacency_relations():
	"""find pairs of adjacent segments
	"""
	filt = ma.D['ma_segment'] != -10
	neighbours_dist, neighbours_idx = ma.get_neighbours_ma(filt)
	pdict = {}

	for i in np.arange(ma.m*2):
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

				if pdict.has_key(pair):
					pdict[pair]+= 1
				else:
					pdict[pair] = 1

	# import ipdb; ipdb.set_trace()
	return pdict

def compute_segment_centers():
	"""Compute avarage coordinate for each segment"""

	segment_dict = {}
	# segment_point_sums = {}
	for i, segment in enumerate(ma.D['ma_segment']):
		# slicing is not copying!
		if segment_dict.has_key(segment):
			segment_dict[segment][0] += 1
			segment_dict[segment][1] += ma.D['ma_coords'][i]
		else:
			segment_dict[segment] = [1, np.copy(ma.D['ma_coords'][i])]

	for key, value in segment_dict.iteritems():
		segment_dict[key][1] = value[1]/value[0]

	
	return segment_dict

def find_relations(ma):
	flip_relations = timeit(find_flip_relations)
	ma.D['seg_link_flip'] = np.zeros(len(flip_relations), dtype = "3int32")
	i=0
	for (s, e), cnt in flip_relations.iteritems():
		ma.D['seg_link_flip'][i] = [s,e,cnt]
		i+=1

	adj_relations = timeit(find_adjacency_relations)
	ma.D['seg_link_adj'] = np.zeros(len(adj_relations), dtype = "3int32")
	i=0
	for (s, e), cnt in adj_relations.iteritems():
		ma.D['seg_link_adj'][i] = [s,e,cnt]
		i+=1

	io_npy.write_npy(INFILE, ma.D, ['seg_link_flip', 'seg_link_adj'])


def view(ma):
	# ref_count = timeit(count_refs)
	min_link_adj = 20

	segment_centers_dict = timeit(compute_segment_centers)

	seg_centers = np.array([v[1] for v in segment_centers_dict.values()], dtype=np.float32)
	seg_cnts = np.array([v[0] for v in segment_centers_dict.values()], dtype=np.float32)
	seg_ids = np.array([k for k in segment_centers_dict.keys()], dtype=np.float32)

	flip_rel_start = np.zeros((len(ma.D['seg_link_flip']),3), dtype=np.float32)
	flip_rel_end = np.zeros((len(ma.D['seg_link_flip']),3), dtype=np.float32)
	i=0
	for s,e in ma.D['seg_link_flip'][:,:2]:
		flip_rel_start[i] = segment_centers_dict[s][1]
		flip_rel_end[i] = segment_centers_dict[e][1]
		i+=1

	adj_rel_start = np.zeros((len(ma.D['seg_link_adj']),3), dtype=np.float32)
	adj_rel_end = np.zeros((len(ma.D['seg_link_adj']),3), dtype=np.float32)
	i=0
	f = ma.D['seg_link_adj'][:,2] > min_link_adj
	for s,e in ma.D['seg_link_adj'][:,:2][f]:
		adj_rel_start[i] = segment_centers_dict[s][1]
		adj_rel_end[i] = segment_centers_dict[e][1]
		i+=1

	max_r=190.
	c = App()

	c.add_data_source(
		opts=['splat_disk', 'with_normals'],
		points=ma.D['coords'], normals=ma.D['normals']
	)

	if ma.D.has_key('ma_segment'):
		f = np.logical_and(ma.D['ma_radii'] < max_r, ma.D['ma_segment']>0)
		c.add_data_source(
			opts=['splat_point', 'with_intensity'],
			points=ma.D['ma_coords'][f], 
			category=ma.D['ma_segment'][f].astype(np.float32),
			colormap='random'
		)
	
		f = np.logical_and(ma.D['ma_radii'] < max_r, ma.D['ma_segment']==0)
		c.add_data_source(
			opts = ['splat_point', 'blend'],
			points=ma.D['ma_coords'][f]
		)
	else:
		f = ma.D['ma_radii_in'] < max_r
		c.add_data_source(
			opts = ['splat_point', 'blend'],
			points=ma.D['ma_coords_in'][f]
		)
		f = ma.D['ma_radii_out'] < max_r
		c.add_data_source(
			opts = ['splat_point', 'blend'],
			points=ma.D['ma_coords_out'][f]
		)

	f = seg_cnts!=1
	c.add_data_source(
		opts = ['splat_point'],
		points = seg_centers[f]
	)
	if len(flip_rel_start)>0:
		c.add_data_source_line(
			coords_start = flip_rel_start,
			coords_end = flip_rel_end
		)

	if len(adj_rel_start)>0:
		f = seg_cnts!=1
		c.add_data_source_line(
			coords_start = adj_rel_start,
			coords_end = adj_rel_end,
			color = (0,1,0)
		)

	# f = ref_count > 20
	# c.add_data_source(
	# 	opts = ['splat_point', 'fixed_color'],
	# 	points = ma.D['coords'][f],
	# 	# intensity = np.clip(ref_count,0,15).astype(np.float32)[f]
	# 	color = (1,1,1)
	# )

	f = ma.D['ma_radii'] < max_r
	c.add_data_source_line(
		coords_start = ma.D['ma_coords'][f],
		coords_end = ma.D['ma_bisec'][f]+ma.D['ma_coords'][f]
	)
	c.add_data_source_line(
		coords_start = ma.D['ma_coords'][f],
		coords_end = np.concatenate([ma.D['coords'],ma.D['coords']])[f]
	)
	c.add_data_source_line(
		coords_start = ma.D['ma_coords'][f],
		coords_end = np.concatenate([ma.D['coords'][ma.D['ma_qidx_in']],ma.D['coords'][ma.D['ma_qidx_out']]])[f]
	)
	
	c.run()

if __name__ == '__main__':
	if len(sys.argv)>1:
		INFILE = sys.argv[1]
	datadict = io_npy.read_npy(INFILE)
	ma = MAHelper(datadict, origin=True)

	find_relations(ma)
	view(ma)
