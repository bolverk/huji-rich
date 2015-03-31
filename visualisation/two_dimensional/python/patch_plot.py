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
    res['log_abs_velocity'] = numpy.log(1e-30+res['abs_velocity'])
    res['log_pressure'] = numpy.log(res['pressure'])
    res['log_density'] = numpy.log(res['density'])
    res['radius'] = numpy.sqrt(res['x_coordinate']**2+res['y_coordinate']**2)
    res['r_velocity'] = (res['x_velocity']*res['x_coordinate']+
                         res['y_velocity']*res['y_coordinate'])/res['radius']
    if 'entropy' in f:
        res['log_entropy'] = numpy.log(res['entropy'])
    f.close()
    return res

def h5_unite(f_list, g=5./3.):

    import numpy

    temp = [consolidate(fname) for fname in f_list]
    res = {}
    for field in temp[0]:
        res[field] = numpy.concatenate([itm[field] for itm in temp])

    return res

def get_user_args():

    from argparse import ArgumentParser

    parser = ArgumentParser(description='Displays a snapshot of a hydrodynamic simulation')
    #parser.add_argument('file_name',
    #                    help='path to output file')
    parser.add_argument('file_names', metavar='fn', nargs='+',
                        help='input files')
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

    #raw = consolidate(fname)
    raw = h5_unite(fname)
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

    display_snapshot(args.file_names, args.field)

if __name__ == '__main__':

    main()
