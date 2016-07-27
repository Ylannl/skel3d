from mapy.segmentation import perform_segmentation_bisec
from mapy.io import npy
from mapy.util import MAHelper
import argparse

# INFILE = 'data/scan_npy'
INFILE = "/Users/ravi/git/masbcpp/rdam_blokken_npy"
# INFILE = "/Volumes/Data/Data/pointcloud/AHN2_matahn_samples/ringdijk_opmeer_npy"
# INFILE = "/Volumes/Data/Data/pointcloud/AHN2_matahn_samples/denhaag_a12_npy"

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Region growing for MAT sheets')
	parser.add_argument('infile', help='npy file', default=INFILE)
	parser.add_argument('-k', help='Number of neighbours to use', default=10, type=int)
	parser.add_argument('-b', '--bisec_thres', help='Bisec threshold in degrees', default=10.0, type=float)
	parser.add_argument('-t', '--theta_thres', help='Theta threshold in degrees', default=10.0, type=float)
	parser.add_argument('-m', '--mincount', help='Minimum number of points in a segment', default=10, type=int)
	# parser.add_argument('--topo', help='Also compute topological links', dest='topo', action='store_true')
	parser.add_argument('--only_interior', help='Only use interior MAT', dest='interior', action='store_true')
	args = parser.parse_args()

	D = npy.read(args.infile)
	mah = MAHelper(D)
	# import ipdb; ipdb.set_trace()
	g = perform_segmentation_bisec(mah, **args.__dict__)
	# if args.topo:
		# find_relations(mah, infile=args.infile, only_interior=args.interior)