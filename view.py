from time import time
import sys

import numpy as np
from povi import App
from pointio import io_npy
from ma_util import MAHelper

import click

@click.command(help='Visualiser for the MAT')
@click.argument('input', type=click.Path(exists=True))
@click.option('-r', '--max_r', default=190., type=float, help='Only show MAT with a radius lower than specified with this value.')
def mat(input, max_r):
    c = App()
    
    t0=time()
    datadict = io_npy.read_npy(input)
    ma = MAHelper(datadict)
    print("{} points loaded from file in {} s".format(ma.m, time()-t0))

    # import ipdb; ipdb.set_trace()
    c.add_data_source(
        opts=['splat_disk', 'with_normals'],
        points=ma.D['coords'], normals=ma.D['normals']
    )

    if ma.D.has_key('ma_segment'):
        f = np.logical_and(ma.D['ma_radii'] < max_r, ma.D['ma_segment']>0)
        c.add_data_source(
            opts=['splat_point', 'with_intensity'],
            points=ma.D['ma_coords'][f], 
            category=ma.D['ma_segment'][f].astype(np.float32),
            colormap='random'
        )
    
        f = np.logical_and(ma.D['ma_radii'] < max_r, ma.D['ma_segment']==0)
        c.add_data_source(
            opts = ['splat_point', 'blend'],
            points=ma.D['ma_coords'][f]
        )
    else:
        f = ma.D['ma_radii_in'] < max_r
        c.add_data_source(
            opts = ['splat_point', 'blend'],
            points=ma.D['ma_coords_in'][f]
        )
        f = ma.D['ma_radii_out'] < max_r
        c.add_data_source(
            opts = ['splat_point', 'blend'],
            points=ma.D['ma_coords_out'][f]
        )

    f = ma.D['ref_count'] > 10
    c.add_data_source(
        opts = ['splat_point', 'with_intensity'],
        points = ma.D['coords'][f],
        intensity = np.clip(ma.D['ref_count'],0,15).astype(np.float32)[f]
    )

    f = ma.D['ma_radii'] < max_r
    c.add_data_source_line(
        coords_start = ma.D['ma_coords'][f],
        coords_end = ma.D['ma_bisec'][f]+ma.D['ma_coords'][f]
    )
    c.add_data_source_line(
        coords_start = ma.D['ma_coords'][f],
        coords_end = np.concatenate([ma.D['coords'],ma.D['coords']])[f]
    )
    c.add_data_source_line(
        coords_start = ma.D['ma_coords'][f],
        coords_end = np.concatenate([ma.D['coords'][ma.D['ma_qidx_in']],ma.D['coords'][ma.D['ma_qidx_out']]])[f]
    )
    
    c.run()

if __name__ == '__main__':
    mat()