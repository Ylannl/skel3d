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
import numpy as np
try:
	import pypcd
except ImportError:
	print("Cannot read pcd files without pypcd module")
	raise

def write(dir, datadict, keys=[]):
	# limited to 2d float arrays now
	if not os.path.exists(dir):
	    os.makedirs(dir)

	for key,val in list(datadict.items()):
		if key in ['coords', 'ma_coords_in', 'ma_coords_out']:
			if key in keys or len(keys)==0:
				fname = os.path.join(dir,key+'.pcd')
				pc = pypcd.make_xyz_point_cloud(datadict[key])
				pc.save_pcd(fname)

def read(dir, keys=[]):
	assert os.path.exists(dir)

	if len(keys)==0:
		keys = inspect(dir)

	datadict = {}	
	for key in keys:
		fname = os.path.join(dir,key+'.pcd')
		if os.path.exists(fname) and key in ['coords', 'ma_coords_in', 'ma_coords_out']:
			pc = pypcd.PointCloud.from_path(fname)
			datadict[key] = np.empty((pc.points,3),dtype=np.float32)
			datadict[key][:,0] = pc.pc_data['x']
			datadict[key][:,1] = pc.pc_data['y']
			datadict[key][:,2] = pc.pc_data['z']
	return datadict

def inspect(dir):
	from glob import glob
	dir = os.path.join(dir,'*.pcd')
	return [os.path.split(f)[-1][:-4] for f in glob(dir)]
