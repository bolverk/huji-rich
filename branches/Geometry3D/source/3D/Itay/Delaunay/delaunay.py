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
for i in range(len(dt.simplices)):
    print(i)
    print("\t%d,%d,%d,%d" % tuple(dt.simplices[i]))
    print("\t%d,%d,%d,%d" % tuple(dt.neighbors[i]))
    print()