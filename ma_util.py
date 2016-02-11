# from pointio import io_npy
import numpy as np
from pykdtree.kdtree import KDTree

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
        if datadict.has_key('ma_segment'):
            self.D['ma_segment'] = datadict['ma_segment']
        self.m, self.n = self.D['coords'].shape
        self.D['ma_qidx_in'] = datadict['ma_qidx_in']
        self.D['ma_qidx_out'] = datadict['ma_qidx_out']
        self.D['ma_radii_in'] = np.linalg.norm(self.D['coords'] - self.D['ma_coords_in'], axis=1)
        self.D['ma_radii_out'] = np.linalg.norm(self.D['coords'] - self.D['ma_coords_out'], axis=1)

        f1_in = self.D['coords']-self.D['ma_coords_in']
        f2_in = self.D['coords'][self.D['ma_qidx_in']]-self.D['ma_coords_in']
        f1_in = f1_in/np.linalg.norm(f1_in, axis=1)[:,None]
        f2_in = f2_in/np.linalg.norm(f2_in, axis=1)[:,None]
        self.D['ma_bisec_in'] = (f1_in+f2_in)
        self.D['ma_bisec_in'] = self.D['ma_bisec_in']/np.linalg.norm(self.D['ma_bisec_in'], axis=1)[:,None]
        self.D['ma_theta_in'] = np.arccos(np.sum(f1_in*f2_in,axis=1))

        f1_out = self.D['coords']-self.D['ma_coords_out']
        f2_out = self.D['coords'][self.D['ma_qidx_out']]-self.D['ma_coords_out']
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

    def get_neighbours_ma(self, filt, k=15):
        self.kdt_ma = KDTree(self.D['ma_coords'][filt])
        return self.kdt_ma.query(self.D['ma_coords'][filt], k)