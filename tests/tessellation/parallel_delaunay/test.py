def parse_delaunay(fname):

    import h5py
    import numpy

    res = {}
    with h5py.File(fname,'r+') as f:
        for field in f:
            res[field] = numpy.array(f[field])
    #res['facets'] = numpy.reshape(res['triangles'],(3,-1)).T
    res['facets'] = numpy.reshape(res['triangles'],(-1,3))
    return res

def main():

    import os
    import pylab
    import numpy
    import glob

    if os.path.isfile('serial_ignore.txt'):
        return True

    whole = parse_delaunay('whole.h5')

    pylab.subplot(211)
    pylab.plot(whole['x_coordinate'][range(whole['point number'])],
               whole['y_coordinate'][range(whole['point number'])],
               '.')
    for facet in whole['facets']:
        closed = numpy.concatenate((facet,[facet[0]]))
        pylab.plot(whole['x_coordinate'][closed],
                   whole['y_coordinate'][closed],
                   color='g')

    pylab.subplot(212)
    for fname in glob.glob('part_*.h5'):
        data = parse_delaunay(fname)
        pylab.plot(
            data['x_coordinate'][range(data['point number'])],
            data['y_coordinate'][range(data['point number'])],
            '.')
        #for facet in data['facets']:
        #    closed = numpy.concatenate((facet,[facet[0]]))
        #    pylab.plot(data['x_coordinate'][closed],
        #               data['y_coordinate'][closed],
        #               color='g') 
    pylab.show()

    return False

if __name__ == '__main__':

    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
