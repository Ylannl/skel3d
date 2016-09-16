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

import numpy as np
import os

def read(infile, move_to_origin=False, limit_points=0, delimiter=' '):
	"""collect vertex coordinates from ascii input file"""
	ox,oy,oz = (0,0,0)
	datadict = {}

	linecount=0
	with open(infile) as f:
		for l in f:
			linecount+=1
    
	with open(infile) as f:
		datadict['coords'] = np.empty((linecount,3), dtype=np.float32)
		for i, line in enumerate(f):
			columns = line.split(delimiter)
			x = float(columns[0])
			y = float(columns[1])
			z = float(columns[2])

			if limit_points != 0 and i > limit_points:
				break

			if move_to_origin and i==0:
				ox,oy,oz = float(x), float(y), float(z)
				datadict['offset'] = np.zeros(3, dtype=np.double)
				datadict['offset'][0] = ox
				datadict['offset'][1] = oy
				datadict['offset'][2] = oz

			datadict['coords'][i] = x-ox,y-oy,z-oz

	return datadict

def write(dir, datadict, keys=[]):
	if not os.path.exists(dir):
	    os.makedirs(dir)

	for key,val in list(datadict.items()):
		if key == 'ma_segment_graph':
			continue
		elif key in keys or len(keys)==0:
			fname = os.path.join(dir,key) + '.xyz'
			with open(fname, 'w') as fo:
				for item in val:
					line = ' '.join([str(i) for i in item]) + '\n'
					fo.write(line)


