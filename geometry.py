import numpy as np
import random
from geom3d import *

def ransac_plane_fit(point_array, point_indices, threshold=0.05, n_planes=2, max_iterations=30):
    # choose 3 random points

    idx = list(point_indices)

    # we are searching for n_planes solutions    
    solutions = []
    for i in range(n_planes):
        iteration_cnt = 0
        candidate_solutions = []
        while iteration_cnt < max_iterations:
            random.shuffle(idx)
            seeds = idx[:3]
            candidate_plane = Plane( [Point( point_array(pi) for pi in seeds )] )

            inliers = seeds
            for pi in idx[3:]:
                if candidate_plane.distance_to( Point(point_array[pi]) ) < threshold:
                    inliers.append(pi)
            candidate_solutions.append( (candidate_plane, inliers) )

            iteration_cnt+=1
        # determine best Plane
        candidate_solutions.sort(key=lambda tup: len(tup[1]), reverse=True)
        plane, inliers = candidate_solutions[0]

        idx = list(set(idx) - set(inliers))
        solutions.append(plane)
    
    return solutions

def 