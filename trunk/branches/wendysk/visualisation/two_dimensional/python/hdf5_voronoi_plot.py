#def voronoi_plot(x_list,
#                 y_list,
#                 z_list,
#                 vert_x_raw,
#                 vert_y_raw,
#                 vert_n_raw):    

def main():

    import h5py
    import pylab
    import numpy
    import matplotlib
    import sys
    import argparse

    parser = argparse.ArgumentParser(description='Displays snapshot of the hydrodynamic simualation')
    parser.add_argument("file_name",
                        help="path to snapshot file")
    parser.add_argument("field",
                        help="Name of hydrodynamic variable")
    args = parser.parse_args()

    h5f = h5py.File(args.file_name)
    x_list = numpy.array(h5f['x_coordinate'])
    y_list = h5f['y_coordinate']
    if 'log_'==args.field[:4]:
        z_list = numpy.log(h5f[args.field.replace('log_','')])
    elif 'log10_'==args.field[:6]:
        z_list = numpy.log10(h5f[args.field.replace('log10_','')])
        z_list = [max(z,-3) for z in z_list]
    elif args.field=='r_velocity':
        z_list = numpy.array([(vx*x+vy*y)/numpy.sqrt(x**2+y**2) for x,y,vx, vy
                              in zip(x_list,y_list,h5f['x_velocity'],h5f['y_velocity'])])
    elif args.field=='q_velocity':
        z_list = numpy.array([(vx*y-vy*x)/numpy.sqrt(x**2+y**2) for x,y,vx, vy
                              in zip(x_list,y_list,h5f['x_velocity'],h5f['y_velocity'])])
    else:
        z_list = h5f[args.field]
    vert_x_raw = numpy.array(h5f['x position of vertices'])
    vert_y_raw = numpy.array(h5f['y position of vertices'])
    vert_n_raw = h5f['Number of vertices in cell']
    vert_idx_list = numpy.cumsum(vert_n_raw)

    z_range = max(z_list) - min(z_list) + 1e-15
    z_min = min(z_list)
    z_scaled = [(z-z_min)/z_range for z in z_list]

    for i in range(len(vert_n_raw)):
        upbound = vert_idx_list[i]
        if i==0:
            lowbound = 0
        else:
            lowbound = vert_idx_list[i-1]
        pylab.fill(vert_x_raw[lowbound:upbound],
                   vert_y_raw[lowbound:upbound],
                   fc=matplotlib.cm.jet(z_scaled[i]),
                   ec='none')

    time = numpy.array(h5f['time'])[0]
    pylab.title(args.field+' @ t = '+str(time))
    ax, _ = matplotlib.colorbar.make_axes(pylab.gca())
    matplotlib.colorbar.ColorbarBase(ax,norm=matplotlib.colors.Normalize(vmin=min(z_list),vmax=max(z_list)))
    pylab.show()

if __name__ == '__main__':

    main()
