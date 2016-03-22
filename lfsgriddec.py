import numpy as np
from pointio import io_npy
from ma_util import MAHelper
import sys
infile = '/Volumes/Data/Data/pointcloud/dense_matching_kadaster/60_DSM_Cloud_1594_3826_npy/'
infile = '/Users/ravi/git/masbcpp/rdam_blokken_npy'

def lfsgriddec(infile, cellsize=0.25, epsilon=0.4, D=2):
	#read coords into grid dict, store indices in lists
	d = io_npy.read_npy(infile)
	ma = MAHelper(d)
	ma.compute_lfs(only_interior=True)

	x_max, y_max, z_max = ma.D['coords'].max(axis=0)
	x_min, y_min, z_min = ma.D['coords'].min(axis=0)

	x_d, y_d, z_d = x_max-x_min, y_max-y_min, z_max-z_min

	x_count = int(x_d / cellsize) + 1
	y_count = int(y_d / cellsize) + 1
	z_count = int(z_d / cellsize) + 1

	gridic = {}
	for i in xrange(x_count):
		for j in xrange(y_count):
			gridic[(i,j)] = []

	for i,p in enumerate(((ma.D['coords'] - [x_min, y_min, z_min]) / cellsize).astype(int)):
		x,y,z = p
		gridic[(x,y)].append(i)

	#iterate through grid dict and determine thinning factor for each cell, write out thinned points
	decimation_filter = np.zeros(ma.m).astype(bool)
	A = cellsize**2
	for xy, idx in gridic.iteritems():
		n = len(idx)
		if n:
			lfs_mean = np.mean( ma.D['lfs'][idx] )
			d_mean = np.sqrt( A/n )

			target_n = A/(epsilon*lfs_mean)**2
			mask = np.random.random(n) <= (target_n/n)

			decimation_filter[idx] = mask

	ma.D['decimate_lfs'] = decimation_filter
	io_npy.write_npy(infile, ma.D, ['decimate_lfs', 'lfs'])
# import ipdb;ipdb.set_trace()

if __name__ == '__main__':
	infile = sys.argv[-1]
	lfsgriddec(infile)