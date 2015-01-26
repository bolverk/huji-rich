#!/bin/python
"""
Generate 10 random points between (0,0,0) and (100,100,100), save data as both a Python file and a TetGen .node file
"""
import random
import pickle

__author__ = 'zmbq'

def rand_point():
    x = random.randrange(0, 100)
    y = random.randrange(0, 100)
    z = random.randrange(0, 100)
    return x,y,z

points = [rand_point() for i in range(10)]
pickle.dump(points, open('points.pickle', 'wb'))

with open('points.node', 'w') as f:
    f.write('# 10 randomally generated points\n')
    f.write('10 3 0 0\n')
    for i in range(len(points)):
        f.write('%2d %3d %3d %3d\n' % (i+1, points[i][0], points[i][1], points[i][2]))

print("10 points generated and saved in points.pickle and points.node")