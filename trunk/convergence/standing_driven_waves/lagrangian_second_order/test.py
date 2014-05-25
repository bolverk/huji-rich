#! /usr/bin/python

def l1_fit_factor(a1, a2):

    import math

    diff2 = [abs(x-y) for x,y in zip(a1,a2)]
    return sum(diff2)/len(a1)

def main():

    import numpy
    import os
    #import imp
    import math

    x_list, d_list, p_list, vx_list, vy_list = \
        numpy.loadtxt('prof1d.txt', unpack=True)

    rho_0 = numpy.loadtxt('mean_density.txt')
    amp = numpy.loadtxt('amplitude.txt')
    k = 2*numpy.pi/numpy.loadtxt('wavelength.txt')
    vp = numpy.loadtxt('phase_velocity.txt')
    p_0 = numpy.loadtxt('mean_pressure.txt')
    gamma = numpy.loadtxt('adiabatic_index.txt')
    c_0 = numpy.sqrt(gamma*p_0/rho_0)
    t = numpy.loadtxt('time.txt')

    def density_exact(x):

        return rho_0-\
            amp*rho_0*math.cos(k*x)*\
            (vp*math.sin(c_0*k*t)-c_0*math.sin(vp*k*t))/\
            (k*c_0*(c_0-vp)*(c_0+vp))

    def pressure_exact(x):

        return p_0-\
            amp*c_0*rho_0*math.cos(k*x)*\
            (vp*math.sin(c_0*k*t)-c_0*math.sin(k*t*vp))/\
            (k*(c_0-vp)*(c_0+vp))

    def velocity_exact(x):

        return -amp*vp*math.sin(k*x)*\
            (math.cos(k*t*vp)-math.cos(k*t*c_0))/\
            (k*(c_0-vp)*(c_0+vp))

    d_exact = [density_exact(x)
               for x in x_list]
    p_exact = [pressure_exact(x)
               for x in x_list]
    v_exact = [velocity_exact(x)
               for x in x_list]

    if False:
        import pylab
        pylab.subplot(311)
        pylab.plot(x_list, d_list, '.')
        pylab.plot(x_list, d_exact, '.')
        pylab.subplot(312)
        pylab.plot(x_list, p_list, '.')
        pylab.plot(x_list, p_exact, '.')
        pylab.subplot(313)
        pylab.plot(x_list, vx_list, '.')
        pylab.plot(x_list, v_exact, '.')
        pylab.show()

    l1_density = l1_fit_factor(d_list, d_exact)
    l1_pressure = l1_fit_factor(p_list, p_exact)
    l1_velocity = l1_fit_factor(vx_list, v_exact)

    f = open('gradesheet.txt','w')
    f.write(str(l1_density)+'\n')
    f.write(str(l1_pressure)+'\n')
    f.write(str(l1_velocity)+'\n')
    f.close()
    
    return l1_density < 1.27e-5 and \
        l1_pressure < 2.2e-5 and \
        l1_velocity < 3.4e-5

import sys
if __name__=='__main__':
    print main()
