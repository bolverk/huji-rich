def demo2():

    import numpy as np
    from tvtk.api import tvtk
    from mayavi import mlab

    points = np.random.normal(0, 1, (1000, 3))
    ug = tvtk.UnstructuredGrid(points=points)
    ug.point_data.scalars = np.sqrt(np.sum(points**2, axis=1))
    ug.point_data.scalars.name = "value"
    ds = mlab.pipeline.add_dataset(ug)
    delaunay = mlab.pipeline.delaunay3d(ds)
    iso = mlab.pipeline.iso_surface(delaunay)
    iso.actor.property.opacity = 0.1
    iso.contour.number_of_contours = 10
    mlab.show()

def demo():

    import numpy as np
    from tvtk.api import tvtk
    from mayavi import mlab

    data = np.loadtxt('./final.txt')

    points = np.vstack((data.T[0], data.T[1], data.T[2])).T
    ug = tvtk.UnstructuredGrid(points=points)
    ug.point_data.scalars = data.T[3]
    ug.point_data.scalars.name = "density"
    ds = mlab.pipeline.add_dataset(ug)
    delaunay = mlab.pipeline.delaunay3d(ds)
    iso = mlab.pipeline.iso_surface(delaunay)
    iso.actor.property.opacity = 0.1
    iso.contour.number_of_contours = 10
    mlab.show()

def plot2d():

    import numpy
    import pylab

    init = numpy.loadtxt('./initial.txt')
    r_init = numpy.sqrt(init.T[0]**2+
                        init.T[1]**2+
                        init.T[2]**2)                        
    final = numpy.loadtxt('./final.txt')
    r_final = numpy.sqrt(final.T[0]**2+
                         final.T[1]**2+
                         final.T[2]**2)
    pylab.subplot(211)
    pylab.semilogy(r_init, init.T[3], '.')
    pylab.semilogy(r_final, final.T[3], '.')
    pylab.subplot(212)
    pylab.semilogy(r_init, init.T[4], '.')
    pylab.semilogy(r_final, final.T[4], '.')
    pylab.show()

def visualise():

    import numpy
    import pylab
    from tvtk.api import tvtk

    data = numpy.loadtxt('./initial.txt')
    r_list = numpy.sqrt(data.T[0]**2+
                        data.T[1]**2+
                        data.T[2]**2)
    pylab.subplot(211)
    pylab.plot(r_list, data.T[3], '.')
    pylab.subplot(212)
    pylab.plot(r_list, data.T[4], '.')
    pylab.show()

def plot_trajectory():

    import numpy
    import pylab

    data = numpy.loadtxt('shock_trajectory.txt')
    pylab.loglog(data.T[0],
                 data.T[1],
                 '.')
    mask = data.T[0]>1e-2
    fitd = numpy.polyfit(numpy.log(data.T[0])[mask],
                         numpy.log(data.T[1])[mask],
                         1)
    x_list = sorted(data.T[0])
    pylab.loglog(x_list,
                 numpy.exp(numpy.polyval(fitd, numpy.log(x_list))))
    pylab.show()

def main():

    import numpy
    import os

    data = numpy.loadtxt('shock_trajectory.txt')
    mask = data.T[1]>0.05
    fitd = numpy.polyfit(numpy.log(data.T[0])[mask],
                         numpy.log(data.T[1])[mask],
                         1)
    if abs(fitd[0]-0.4)<0.1:
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')

if __name__ == '__main__':

    main()
