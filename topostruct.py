import math
import numpy as np
from pointio import io_npy
from ma_util import MAHelper

# INFILE = 'data/scan_npy'
INFILE = "/Users/ravi/git/masbcpp/rdam_blokken_npy"

datadict = io_npy.read_npy(INFILE)
ma = MAHelper(datadict, origin=False)

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
	import ipdb; ipdb.set_trace()


def ref_count():
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

	datadict['ref_count'] = np.array(pdict.values(), dtype=np.int32)
	io_npy.write_npy(INFILE, datadict, ['ref_count'])
	import ipdb; ipdb.set_trace()

def find_flip_relations():
	"""find for each pair of segments how many times they are connected by a shared feature point

		in a pair of segments (tuple) the lowest segmend_id is always put first
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


	
	# datadict['flip_relations'] = np.array(pdict.values(), dtype=np.int32)
	# io_npy.write_npy(INFILE, datadict, ['ref_count'])
	import ipdb; ipdb.set_trace()

def compute_segment_centers():

	segment_dict = {}
	# segment_point_sums = {}
	for i, segment in enumerate(ma.D['ma_segment']):
		p = ma.D['ma_coords'][i]
		if segment_dict.has_key(segment):
			segment_dict[segment][0] += 1
			segment_dict[segment][1] += p
		else:
			segment_dict[segment] = [1, p]

	for key, value in segment_dict.iteritems():
		segment_dict[key][1] = value[1]/value[0]

	import ipdb; ipdb.set_trace()

if __name__ == '__main__':
	# ref_count()
	# find_flip_relations()
	compute_segment_centers()
