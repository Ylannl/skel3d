from skel3d.segmentation import perform_segmentation_bisec
from skel3d.io import npy
from skel3d.util import MAHelper
from skel3d.segmentation import RegionGrower
import numpy as np
import argparse

# INFILE = 'data/scan_npy'
INFILE = "/Users/ravi/git/masbcpp/rdam_blokken_npy"
# INFILE = "/Volumes/Data/Data/pointcloud/AHN2_matahn_samples/ringdijk_opmeer_npy"
# INFILE = "/Volumes/Data/Data/pointcloud/AHN2_matahn_samples/denhaag_a12_npy"


def list_region_idx(mah):
    region_numbers = np.unique(mah.D['ma_segment'])

    L = []
    for region_number in region_numbers:
        idx = np.where(mah.D['ma_segment']==region_number)
        L.append(idx[0])

    return L

def filter_regions(mah, max_r, minavg_theta):
    region_numbers = np.unique(mah.D['ma_segment'])

    ma_idx = np.zeros(mah.m*2, dtype=np.bool)
    for region_number in region_numbers:
        idx = np.where(mah.D['ma_segment']==region_number)
        if np.nanmax(mah.D['ma_radii'][idx]) < max_r and np.nanmean(mah.D['ma_theta'][idx] > minavg_theta):
            ma_idx[idx] = 1

    return ma_idx

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Region growing for MAT sheets')
    parser.add_argument('infile', help='npy file', default=INFILE)
    parser.add_argument('-k', help='Number of neighbours to use', default=10, type=int)
    parser.add_argument('-b', '--balloverlap_thres', help='Balloverlap threshold', default=2.0, type=float)
    parser.add_argument('-r', '--max_r', help='Filter regions with a maximum radius greater than this value in model coordinates', default=100, type=float)
    parser.add_argument('-s', '--minavg_theta', help='Filter regions with a mean separation angle smaller than this value in radians', default=2, type=float)
    args = parser.parse_args()

    D = npy.read(args.infile)
    mah = MAHelper(D)
    
    R = RegionGrower(mah, method='balloverlap', **args.__dict__)
    seedorder = list( np.random.permutation(R.m) )
    R.apply_region_growing_algorithm(seedorder)
    # R.unmark_small_clusters()
    mah.D['ma_segment']= R.ma_segment
    mah.D['ma_segment_lidx']= list_region_idx(mah)
    npy.write(args.infile, mah.D, ['ma_segment', 'ma_segment_lidx'])
    

    # ma_idx = filter_regions(mah, args.max_r, args.minavg_theta)
    # print(mah.D['ma_coords'][ma_idx])
    # print(mah.D['ma_radii'][ma_idx])


