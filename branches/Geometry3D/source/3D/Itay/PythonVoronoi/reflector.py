from TetGenReader import read_tetgen_points
from boundary import Boundary, reflect_points

__author__ = 'zmbq'

import os.path

def main_reflect_points(folder, name):
    filename = os.path.join(folder, name)
    all_pts = read_tetgen_points(filename)
    big_tetrahedron = all_pts[-4:]
    pts = all_pts[:-4]
    boundary = Boundary((-500, 1000, 300), (-1000, 0, -300))
    offsets = ('-  ', ' - ', '  -', '  +', ' + ', '+  ')
    reflected = reflect_points(pts, boundary, offsets)

    num = 501
    for pt in reflected:
        print("%4d: %s" % (num, pt))
        num += 1

if __name__=='__main__':
    main_reflect_points('d:/chelem/rich/temp/stress', 'orig.node')
