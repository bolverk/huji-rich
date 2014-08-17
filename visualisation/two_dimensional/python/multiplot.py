def unstructured_contour(x_list, y_list, z_list):

    from matplotlib.tri import Triangulation
    from scipy.spatial import Delaunay
    import matplotlib.pyplot as plt
    import numpy

    tri = Delaunay(numpy.array([[x,y] for x,y
                                in zip(x_list, y_list)]))
    mtri = Triangulation(x_list, y_list, tri.vertices.copy())
    res = plt.tricontourf(mtri, z_list)
    return res

def consolidate_data(fname):

    import h5py
    import numpy

    f = h5py.File(fname)
    res = {}
    for field in f:
        res[field] = numpy.array(f[field])
    return res

def collect_snapshot_data(file_list):

    import numpy

    part_list = [consolidate_data(fname) for fname
                 in file_list]
    res = {}
    for field in part_list[0]:
        res[field] = numpy.concatenate([itm[field] for itm in part_list])
    return res

def main():

    import argparse
    import glob
    import pylab

    parser = argparse.ArgumentParser(description='Plots RICH snapshots')
    parser.add_argument('pattern',help='Pattern to match all output files')
    parser.add_argument('varname',help='name of variable to be plotted (z axis)')
    args = parser.parse_args()

    data = collect_snapshot_data(glob.glob(args.pattern))
    unstructured_contour(data['x_coordinate'],
                         data['y_coordinate'],
                         data[args.varname])
    pylab.show()

if __name__ == '__main__':

    main()

    
