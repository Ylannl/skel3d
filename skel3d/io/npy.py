# This file is part of pointio.

# pointio is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# pointio is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with pointio.  If not, see <http://www.gnu.org/licenses/>.

# Copyright 2015 Ravi Peters

import os
from glob import glob
import numpy as np
import igraph

def write(dir, datadict, keys=[]):
	if not os.path.exists(dir):
	    os.makedirs(dir)

	for key,val in list(datadict.items()):
		if key == 'ma_segment_graph':
			if type(datadict[key]) is igraph.Graph:
				datadict[key].write_pickle(os.path.join(dir, key+'.pickle'))
		elif key == 'ma_clusters':
			clusterdir = os.path.join(dir, 'ma_clusters')
			if not os.path.exists(clusterdir):
				os.makedirs(clusterdir)
			else:
				for cluster_file in glob(os.path.join(clusterdir, '*.pickle')):
					os.remove(cluster_file)
			for i, cluster in enumerate(datadict[key]):
				cluster.write_pickle(os.path.join(clusterdir, 'ma_cluster_'+str(i)+'.pickle'))
		elif key in keys or len(keys)==0:
			fname = os.path.join(dir,key)
			np.save(fname, val)

def read(dir, keys=[]):
	assert os.path.exists(dir)

	if len(keys)==0:
		keys = inspect(dir)

	datadict = {}	
	for key in keys:
		fname = os.path.join(dir,key+'.npy')
		if os.path.exists(fname):
			datadict[key] = np.load(fname)

	fname = os.path.join(dir,'ma_segment_graph.pickle')
	if os.path.exists(fname):
		datadict['ma_segment_graph'] = igraph.read(fname)

	fname = os.path.join(dir, 'ma_clusters')
	if os.path.exists(fname):
		datadict['ma_clusters'] = []
		for f in glob(os.path.join(fname, '*.pickle')):
			datadict['ma_clusters'].append(igraph.read(f))

	return datadict

def inspect(dir):
	dir = os.path.join(dir,'*.npy')
	return [os.path.split(f)[-1][:-4] for f in glob(dir)]