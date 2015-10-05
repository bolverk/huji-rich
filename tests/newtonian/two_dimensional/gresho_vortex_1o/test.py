def l1_fit(a1, a2):

    return sum([abs(x-y) for x,y in zip(a1,a2)])/len(a1)

def pressure_dist(r):

    import math

    if r<0.2:
        return 5+(25./2.)*r**2
    elif r>0.4:
        return 3+4*math.log(2.)
    else:
        return 9+(25./2.)*r**2-20*r+4*math.log(r/0.2)

def azimuthal_velocity_dist(r):

    import math

    if r<0.2:
        return 5*r
    elif r>0.4:
        return 0
    else:
        return 2-5*r

def main():

    import numpy
    import math
    import h5py
    import glob

    np = len(glob.glob('process_*_final.h5'))
    if np>0:
        rx_list = []
        ry_list = []
        d_list = []
        p_list = []
        vx_list = []
        vy_list = []
        for fname in glob.glob('process_*_final.h5'):
            f = h5py.File(fname)
            rx_list.extend(f['x_coordinate'])
            ry_list.extend(f['y_coordinate'])
            d_list.extend(f['density'])
            p_list.extend(f['pressure'])
            vx_list.extend(f['x_velocity'])
            vy_list.extend(f['y_velocity'])
    else:
        h5f = h5py.File('final.h5')
        rx_list = h5f['geometry']['x_coordinate']
        ry_list = h5f['geometry']['y_coordinate']
        d_list = h5f['hydrodynamic']['density']
        p_list = h5f['hydrodynamic']['pressure']
        vx_list = h5f['hydrodynamic']['x_velocity']
        vy_list = h5f['hydrodynamic']['y_velocity']

    r_list = [math.sqrt(x**2+y**2)
              for x,y in zip(rx_list,ry_list)]
    r_sorted = numpy.sort(r_list)
    vf_list = [(vy*x-vx*y)/math.sqrt(x**2+y**2)
               for x,y,vx,vy in
               zip(rx_list,
                   ry_list,
                   vx_list,
                   vy_list)] 

    vr_list = [(vy*y+vx*x)/math.sqrt(x**2+y**2)
               for x,y,vx,vy in
               zip(rx_list,
                   ry_list,
                   vx_list,
                   vy_list)] 

    if False:
        import pylab
        pylab.subplot(221)
        pylab.plot(r_list,p_list,'.')
        pylab.plot(r_sorted,[pressure_dist(r) for r in r_sorted])

        pylab.subplot(222)
        pylab.plot(r_list,vf_list,'.')
        pylab.plot(r_sorted,
                   [azimuthal_velocity_dist(r)
                    for r in r_sorted])

        pylab.subplot(223)
        pylab.plot(r_list,vr_list,'.')
        pylab.plot(r_sorted,
                   0*r_sorted)

        pylab.subplot(224)
        pylab.plot(r_list,d_list,'.')
        pylab.plot(r_sorted,
                   [1 for r in r_sorted])

        pylab.show()

    l1_density = l1_fit(d_list, [1 for r in r_list])
    l1_pressure = l1_fit(p_list,
                         [pressure_dist(r) for r in r_list])
    l1_radial_velocity = l1_fit(
        vr_list,[0 for r in r_list])
    l1_azimuthal_velocity = l1_fit(
        vf_list,[azimuthal_velocity_dist(r) for r in r_list])

    f = open('gradesheet.txt','w')
    f.write(str(l1_density)+'\n')
    f.write(str(l1_pressure)+'\n')
    f.write(str(l1_radial_velocity)+'\n')
    f.write(str(l1_azimuthal_velocity)+'\n')
    f.close()

    return l1_density<3.3e-4 and \
        l1_pressure<2.9e-3 and \
        l1_radial_velocity<1.3e-3 and \
        l1_azimuthal_velocity<2.1e-3

if __name__=='__main__':
    
    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')

