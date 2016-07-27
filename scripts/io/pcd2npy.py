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

import sys, argparse
from pointio import io_npy, io_pcd

def main(args):
	keys = []
	if args.key:
		keys.append(args.key)
	datadict = io_pcd.read_pcd(args.infile, keys)
	io_npy.write_npy(args.outfile, datadict, keys)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='convert pcl .pcd files to numpy binary files')
    parser.add_argument('infile', help='input pcd dir')
    parser.add_argument('outfile', help='npy output directory')
    parser.add_argument('-k', '--key', help='read only this key', dest='key', default='')

    args = parser.parse_args()
    main(args)