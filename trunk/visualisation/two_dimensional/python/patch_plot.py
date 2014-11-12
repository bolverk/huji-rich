def consolidate(fname,g=5./3.):

    import h5py
    import numpy

    f = h5py.File(fname)
    res = {}
    for field in f:
        res[field] = numpy.array(f[field])
    res['sound_speed'] = numpy.sqrt(g*res['pressure']/res['density'])
    res['abs_velocity'] = numpy.sqrt(res['x_velocity']**2+res['y_velocity']**2)
    res['mach_number'] = res['abs_velocity']/res['sound_speed']
    f.close()
    return res

def get_user_args():

    from argparse import ArgumentParser

    parser = ArgumentParser(description='Displays a snapshot of a hydrodynamic simulation')
    parser.add_argument('file_name',
                        help='path to output file')
    parser.add_argument('field',
                        help='name of hydrodynamic variable')    
    return parser.parse_args()

def display_snapshot(fname,field):

    import matplotlib
    import numpy
    import pylab
    from matplotlib.collections import PolyCollection
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    raw = consolidate(fname)
    vert_idx_list = numpy.concatenate(([0],numpy.cumsum(raw['Number of vertices in cell'])))

    verts = []
    for i in range(len(raw[field])):
        lowbound = vert_idx_list[i]
        upbound = vert_idx_list[i+1]
        verts.append([[x,y] 
                      for x,y 
                      in zip(raw['x position of vertices'][lowbound:upbound],
                             raw['y position of vertices'][lowbound:upbound])])

    coll = PolyCollection(verts,array=raw[field],cmap=mpl.cm.jet,
                          edgecolors='none')
    fig, ax = plt.subplots()
    fig.suptitle(field+' @ t = '+str(raw['time'][0]))
    ax.add_collection(coll)
    ax.autoscale_view()
    ax.set_aspect('equal')

    fig.colorbar(coll, ax=ax)
    plt.show()

def main():

    args = get_user_args()

    display_snapshot(args.file_name, args.field)

if __name__ == '__main__':

    main()
