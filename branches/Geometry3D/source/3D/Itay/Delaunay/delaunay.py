#!/bin/python3
"""
Calculates the 3D DT of the points in points.pickle
"""

import pickle
from scipy.spatial import Delaunay

__author__ = 'zmbq'

points = pickle.load(open('points.pickle', 'rb'))
print(points[:7])

dt = Delaunay(points[:7])
for simplex in dt.simplices:
    print("%d,%d,%d,%d" % tuple(simplex))