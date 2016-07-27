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
from pointio import io_npy, io_las

def main(args):
	datadict = io_npy.read_npy(args.infile)
	io_las.write_las(args.outfile,datadict)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='convert numpy binary file to LAS')
    parser.add_argument('infile', help='input npy directory')
    parser.add_argument('outfile', help='LAS output path')

    args = parser.parse_args()
    main(args)