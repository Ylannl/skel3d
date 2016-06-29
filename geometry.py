import numpy as np
import random
from geom3d import *

def ransac_plane_fit(point_array, point_indices, threshold=0.05, n_planes=2, max_iterations=10):
    # choose 3 random points

    idx = list(point_indices)

    # we are searching for n_planes solutions, but we may only find a lower number of planes    
    solutions = []
    for i in range(n_planes):
        if len(idx) < 3:
            break
        iteration_cnt = 0
        candidate_solutions = []
        while iteration_cnt < max_iterations:
            random.shuffle(idx)
            seeds = idx[:3]
            # import ipdb; ipdb.set_trace()
            try: 
                candidate_plane = Plane( [Point( point_array[pi]) for pi in seeds] )
            except ValueError:
                continue

            inliers = seeds
            for pi in idx[3:]:
                d = candidate_plane.distance_to( Point(point_array[pi]) )
                if d < threshold:
                    inliers.append(pi)
            candidate_solutions.append( (candidate_plane, inliers) )

            iteration_cnt+=1
        # determine best Plane
        candidate_solutions.sort(key=lambda tup: len(tup[1]), reverse=True)
        plane, inliers = candidate_solutions[0]

        idx = list(set(idx) - set(inliers))
        solutions.append(plane)
    
    return solutions

def derive_tetra(plane1, plane2, point_indices, ma):
    line = plane1.intersection(plane2)

    line_range = []
    for c in ma.D['ma_coords'][point_indices]:
        p = Point(c)
        s = np.inner(line.t,p.r - line.r)
        line_range.append(s)
    s_mi, s_ma = min(line_range), max(line_range)
    v0,v1 = line.t*s_mi+line.r, line.t*s_ma+line.r

    pi_maxr =  np.argmax(ma.D['ma_radii'][point_indices])
    p_maxr = ma.D['ma_coords'][point_indices][pi_maxr]

    v2 = plane1.project(Point(p_maxr)).r
    v3 = plane2.project(Point(p_maxr)).r

    # center v2,v3 nicely to be halfway across v0,v1
    s = np.inner(line.t,v2 - line.r)
    move_to_center = (s_mi + (s_ma-s_mi)/2 - s)*line.t
    v2 += move_to_center
    v3 += move_to_center

    return np.array([v0, v1, v2, v3], dtype=np.float32) 