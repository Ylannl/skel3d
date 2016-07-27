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
	datadict = io_las.read_las(args.infile, move_to_origin=args.move_to_origin)
	io_npy.write_npy(args.outfile, datadict, ['coords', 'offset'])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='convert LAS to numpy binary file')
    parser.add_argument('infile', help='input .las')
    parser.add_argument('outfile', help='npy output directory')
    parser.add_argument('-d', '--dont-move_to_origin', help='Don\'t move points to origin', dest='move_to_origin', action='store_false')

    args = parser.parse_args()
    main(args)