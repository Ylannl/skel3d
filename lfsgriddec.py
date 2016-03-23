import numpy as np
from pointio import io_npy
from ma_util import MAHelper
import sys
infile = '/Volumes/Data/Data/pointcloud/dense_matching_kadaster/60_DSM_Cloud_1594_3826_npy/'
infile = '/Users/ravi/git/masbcpp/rdam_blokken_npy'

def lfsgriddec(infile, cellsize=0.5, epsilon=0.4, D=2, edge_preserve=True):
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
	if D==2:
		for i in xrange(x_count):
			for j in xrange(y_count):
				gridic[(i,j)] = []

		for i,p in enumerate(((ma.D['coords'] - [x_min, y_min, z_min]) / cellsize).astype(int)):
			x,y,z = p
			gridic[(x,y)].append(i)
	elif D==3:
		for i in xrange(x_count):
			for j in xrange(y_count):
				for k in xrange(z_count):
					gridic[(i,j,k)] = []

		for i,p in enumerate(((ma.D['coords'] - [x_min, y_min, z_min]) / cellsize).astype(int)):
			x,y,z = p
			gridic[(x,y,z)].append(i)

	#iterate through grid dict and determine thinning factor for each cell, write out thinned points
	decimation_filter = np.zeros(ma.m).astype(bool)
	A = cellsize**2
	for idx in gridic.values():
		n = len(idx)
		if n:
			lfs_mean = np.mean( ma.D['lfs'][idx] )
			
			# d_mean = np.sqrt( A/n )
			if edge_preserve:
				elevation_diff = np.ptp(ma.D['coords'][idx][:,2],0)
				if elevation_diff > 0.5:
					lfs_mean/=10

			target_n = A/(epsilon*lfs_mean)**D
			mask = np.random.random(n) <= (target_n/n)

			decimation_filter[idx] = mask

			if edge_preserve:
				ma.D['lfs'][idx]/=10

	
	if edge_preserve:
		ma.D['decimate_lfs_edge'] = decimation_filter
		ma.D['lfs_edge'] = ma.D['lfs']
		io_npy.write_npy(infile, ma.D, ['decimate_lfs_edge', 'lfs_edge'])
	else:
		ma.D['decimate_lfs'] = decimation_filter
		io_npy.write_npy(infile, ma.D, ['decimate_lfs', 'lfs'])
# import ipdb;ipdb.set_trace()

if __name__ == '__main__':
	infile = sys.argv[-1]
	lfsgriddec(infile)