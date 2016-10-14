# from pointio import io_npy
from numpy import array, zeros, empty, invert, concatenate, nanmax, isnan
import numpy as np

def angle(a, b):
    a = a/np.linalg.norm(a)
    b = b/np.linalg.norm(b)
    return np.arccos(np.inner(a, b))

class MAHelper(object):

    def __init__(self, datadict, origin=True):
        
        self.D=datadict
        self.D['coords'] = datadict['coords']
        self.D['ma_coords_in'] = datadict['ma_coords_in']
        self.D['ma_coords_out'] = datadict['ma_coords_out']
        if origin==True:
            self.mean = np.mean(datadict['coords'], axis=0, dtype=np.float32)
            self.D['coords'] -= self.mean
            self.D['ma_coords_in'] -= self.mean
            self.D['ma_coords_out'] -= self.mean
        self.D['normals'] = datadict['normals']
        if 'ma_segment' in datadict:
            self.D['ma_segment'] = datadict['ma_segment']
        self.m, self.n = self.D['coords'].shape
        self.D['ma_qidx_in'] = datadict['ma_qidx_in']
        self.D['ma_qidx_out'] = datadict['ma_qidx_out']
        self.D['ma_radii_in'] = np.linalg.norm(self.D['coords'] - self.D['ma_coords_in'], axis=1)
        self.D['ma_radii_out'] = np.linalg.norm(self.D['coords'] - self.D['ma_coords_out'], axis=1)

        f1_in = self.D['coords']-self.D['ma_coords_in']
        f2_in = self.D['coords'][self.D['ma_qidx_in']]-self.D['ma_coords_in']
        self.D['ma_f1_in'] = f1_in.copy() 
        self.D['ma_f2_in'] = f2_in.copy()
        f1_in = f1_in/np.linalg.norm(f1_in, axis=1)[:,None]
        f2_in = f2_in/np.linalg.norm(f2_in, axis=1)[:,None]
        self.D['ma_bisec_in'] = (f1_in+f2_in)
        self.D['ma_bisec_in'] = self.D['ma_bisec_in']/np.linalg.norm(self.D['ma_bisec_in'], axis=1)[:,None]
        self.D['ma_theta_in'] = np.arccos(np.sum(f1_in*f2_in,axis=1))

        f1_out = self.D['coords']-self.D['ma_coords_out']
        f2_out = self.D['coords'][self.D['ma_qidx_out']]-self.D['ma_coords_out']
        self.D['ma_f1_out'] = f1_out.copy() 
        self.D['ma_f2_out'] = f2_out.copy()
        f1_out = f1_out/np.linalg.norm(f1_out, axis=1)[:,None]
        f2_out = f2_out/np.linalg.norm(f2_out, axis=1)[:,None]
        self.D['ma_bisec_out'] = (f1_out+f2_out)
        self.D['ma_bisec_out'] = self.D['ma_bisec_out']/np.linalg.norm(self.D['ma_bisec_out'], axis=1)[:,None]
        self.D['ma_theta_out'] = np.arccos(np.sum(f1_out*f2_out,axis=1))

        self.D['ma_coords'] = np.concatenate([self.D['ma_coords_in'], self.D['ma_coords_out']])
        self.D['ma_bisec'] = np.concatenate([self.D['ma_bisec_in'], self.D['ma_bisec_out']])
        self.D['ma_theta'] = np.concatenate([self.D['ma_theta_in'], self.D['ma_theta_out']])
        self.D['ma_radii'] = np.concatenate([self.D['ma_radii_in'], self.D['ma_radii_out']])
        self.D['ma_qidx'] = np.concatenate([self.D['ma_qidx_in'], self.D['ma_qidx_out']])
        self.D['ma_f1'] = np.concatenate([self.D['ma_f1_in'], self.D['ma_f1_out']])
        self.D['ma_f2'] = np.concatenate([self.D['ma_f2_in'], self.D['ma_f2_out']])

        self.D['spoke_cnt'] = zeros(self.m, dtype=int)
        s_idx, cnts = np.unique(self.D['ma_qidx'], return_counts=True)
        self.D['spoke_cnt'][s_idx] = cnts

        if 'ma_segment_graph' in datadict:
            self.g = datadict['ma_segment_graph']

        self.filtered = {}
        self.reset_filter()

    def reset_filter(self):
        self.filtered['in'] = zeros(self.m) == True
        self.filtered['out'] = zeros(self.m) == True
        
    def compute_lfs(self, k=10, only_interior=False):
        from pykdtree.kdtree import KDTree
        # collect all ma_coords that are not NaN
        if only_interior:
            ma_coords = self.D['ma_coords_in'][invert(self.filtered['in'])]
        else:
            ma_coords = concatenate([self.D['ma_coords_in'][invert(self.filtered['in'])], self.D['ma_coords_out'][invert(self.filtered['out'])]])
        ma_coords = ma_coords[~np.isnan(ma_coords).any(axis=1)]

        kd_tree = KDTree(ma_coords)
        # we can get *squared* distances for free, so take the square root
        if k > 1:
            self.D['lfs'] = np.sqrt(np.median(kd_tree.query(self.D['coords'], k)[0], axis=1))
        else:
            self.D['lfs'] = np.sqrt(kd_tree.query(self.D['coords'], k)[0])

    # def decimate_lfs(self, m, max_pointspacing=None, scramble = False, sort = False, squared=False):

    #     if not hasattr(self, 'flann_tree'):
    #         self.flann_tree = FLANN()
    #         self.flann_tree.build_index(self.D['coords'])#, algorithm='linear'
    #     self.D['decimate_lfs'] = zeros(self.m) == True

    #     order = np.arange(self.m)
    #     plfs = list(zip(order, self.D['coords'], self.D['lfs']))

    #     if scramble: 
    #         from random import shuffle
    #         shuffle( plfs )
    #     if sort:
    #         plfs.sort(key = lambda item: item[2])
    #         plfs.reverse()

    #     for i, p, lfs in plfs:
    #         if squared:
    #             lfs = lfs**2
    #         if max_pointspacing is None:
    #             rho = lfs*m
    #         else:
    #             rho = min(lfs*m, max_pointspacing)

    #         qts = self.flann_tree.nn_radius(p, rho**2)[0][1:]
    #         qts = order[qts]
            
    #         # qts = self.flann_tree.nn_index(p, 6)[0][1:]
            
    #         iqts = invert(self.D['decimate_lfs'][qts])
    #         if iqts.any():
    #             self.D['decimate_lfs'][i] = True