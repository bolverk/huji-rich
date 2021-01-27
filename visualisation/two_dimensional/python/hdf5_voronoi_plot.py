def main():

    import h5py
    import pylab
    import numpy
    import matplotlib
    import sys
    import argparse
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection

    parser = argparse.ArgumentParser(description='Displays snapshot of the hydrodynamic simualation')
    parser.add_argument("file_name",
                        help="path to snapshot file")
    parser.add_argument("field",
                        help="Name of hydrodynamic variable")
    args = parser.parse_args()

    with h5py.File(args.file_name, 'r') as f:
        x_list = numpy.array(f['geometry']['x_coordinate'])
        y_list = numpy.array(f['geometry']['y_coordinate'])
        z_list = numpy.array(f['hydrodynamic'][args.field])
        vert_x_raw = numpy.array(f['geometry']['x_vertices'])
        vert_y_raw = numpy.array(f['geometry']['y_vertices'])
        vert_n_raw = numpy.array(f['geometry']['n_vertices'])
        vert_idx_list = numpy.concatenate(([0], numpy.cumsum(vert_n_raw))).astype(int)
        time = numpy.array(f['time'])[0]

    polygon_list = [Polygon(
        numpy.vstack((vert_x_raw[low:high],
                      vert_y_raw[low:high])).T)
                    for low, high
                    in zip(vert_idx_list[:-1],
                           vert_idx_list[1:])]
    patch_collection = PatchCollection(polygon_list)
    patch_collection.set_array(z_list)

    fig, ax = pylab.subplots()
    ax.add_collection(patch_collection)
    pylab.suptitle('t = %.4f' % time)
    pylab.axis('scaled')
    pylab.xlim((numpy.min(x_list), numpy.max(x_list)))
    pylab.ylim((numpy.min(y_list), numpy.max(y_list)))

    pylab.show()

if __name__ == '__main__':

    main()
