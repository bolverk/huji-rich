def parse_delaunay(fname):

    import h5py
    import numpy

    res = {}
    with h5py.File(fname,'r+') as f:
        for field in f:
            res[field] = numpy.array(f[field])
    res['facets'] = numpy.reshape(res['triangles'],(-1,3))
    return res

def locate_point(x_list, y_list, xp, yp, tol = 1e-10):

    import numpy

    d_list = (x_list-xp)**2 + (y_list-yp)**2
    assert(numpy.min(d_list)<tol)
    return numpy.argmin(d_list)

def is_boxed(x, y, offset=0):

    return 1-offset>=x>=0+offset and 1-offset>=y>=0+offset

def main():

    import os
    #import pylab
    import numpy
    import re
    import glob

    if os.path.isfile('serial_ignore.txt'):
        return True

    whole = parse_delaunay('whole.h5')
    parts = [parse_delaunay(fname) for fname in 
             sorted(glob.glob('part_*.h5'),
                    key=
                    lambda x: int(re.search(r'part_(\d+).h5',x).group(1)))]

    whole['sorted facets'] = [sorted(itm) for itm in whole['facets']]
    for n, part in enumerate(parts):
        index_table = [
            locate_point(
                whole['x_coordinate'],
                whole['y_coordinate'],
                x,y) if is_boxed(x,y)
            else -1
            for x,y in 
            zip(part['x_coordinate'],part['y_coordinate'])]
        for facet in part['facets']:
            if numpy.all([itm<part['point number'] for itm in facet]):
                t_facet = sorted([index_table[itm] for itm in facet])
                if not t_facet in whole['sorted facets']:
                    return False
            elif numpy.all(
                [is_boxed(part['x_coordinate'][itm],
                          part['y_coordinate'][itm],0.04)
                 for itm in facet]) and numpy.any([itm<part['point number'] for itm in facet]):
                t_facet = sorted([index_table[itm] for itm in facet])
                if not t_facet in whole['sorted facets']:
                    return False                

    return True

if __name__ == '__main__':

    import os

    os.system('rm *.res;')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
