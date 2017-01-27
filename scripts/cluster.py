from skel3d.io import npy
from skel3d.graph import *
from skel3d.util import MAHelper
from skel3d import clustering

import argparse

# INFILE = 'data/scan_npy'
INFILE = "/Users/ravi/git/masbcpp/rdam_blokken_npy"
# INFILE = "/Volumes/Data/Data/pointcloud/AHN2_matahn_samples/ringdijk_opmeer_npy"
# INFILE = "/Volumes/Data/Data/pointcloud/AHN2_matahn_samples/denhaag_a12_npy"

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Clusterer of MAT sheets')
    parser.add_argument('infile', help='npy file', default=INFILE)
    parser.add_argument('-m', '--mincount', help='Minimum edge count used during connected compenent analysis', default=20, type=int)
    # parser.add_argument('-c', '--contract', help='contract edges', default=None, type=float)
    # parser.add_argument('-a', '--analyse', help='Also compute statistics for each sheet', dest='analyse', action='store_true')
    args = parser.parse_args()

    D = npy.read(args.infile)
    mah = MAHelper(D)

    # if not args.contract is None:
    #     contract_edges(mah.D['ma_segment_graph'], args.contract)
    clustering.get_clusters(mah, args.mincount)
    clustering.classify_clusters(mah)
    for g in mah.D['ma_clusters']:
        if g['classification'] == 4: # building class
            clustering.analyse_cluster(mah, g)

    npy.write(args.infile,mah.D, ['ma_clusters'])