#! /usr/bin/python

def l1_fit_factor(a1, a2):

    import math

    diff2 = [abs(x-y) for x,y in zip(a1,a2)]
    return sum(diff2)/(max(a1)-min(a1))/len(a1)

def main():

    import numpy
    import os
    import imp
    import h5py

    h5f = h5py.File('final.h5','r+')
    x_list = h5f['geometry']['x_coordinate']
    d_list = h5f['hydrodynamic']['density']
    p_list = h5f['hydrodynamic']['pressure']
    vx_list = h5f['hydrodynamic']['x_velocity']

    rho_0 = 1.0
    amp = 0.001
    k = 2*numpy.pi/1.0
    v_0 = 0.1
    omega = k*v_0
    p_0 = 1.0
    gamma = 5./3.
    c_0 = numpy.sqrt(gamma*p_0/rho_0)
    t = numpy.loadtxt('time.txt')

    def density_exact(x):

        return rho_0-\
            (amp*rho_0*(2*c_0*numpy.sin(k*x) - (c_0 + v_0)*numpy.sin(k*(c_0*t - \
t*v_0 + x)) + (-c_0 + v_0)*numpy.sin(k*(-(t*(c_0 + v_0)) + \
x))))/(2.*c_0*k*(c_0 - v_0)*(c_0 + v_0))

    def pressure_exact(x):

        return p_0-\
            (amp*c_0*rho_0*(2*c_0*numpy.sin(k*x) + (-c_0 + \
v_0)*numpy.sin(k*(c_0*t - t*v_0 + x)) - (c_0 + \
v_0)*numpy.sin(k*(-(t*(c_0 + v_0)) + x))))/(2.*k*(c_0 - v_0)*(c_0 + \
v_0))

    def velocity_exact(x):

        return v_0-\
            (amp*(-2*v_0*numpy.sin(k*x) + (c_0 + v_0)*numpy.sin(k*(c_0*t - t*v_0 \
+ x)) + (-c_0 + v_0)*numpy.sin(k*(-(t*(c_0 + v_0)) + x))))/(2.*k*(c_0 \
- v_0)*(c_0 + v_0))

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
    
    return l1_density < 0.07 and \
        l1_pressure < 0.1 and \
        l1_velocity < 0.23

import sys
if __name__=='__main__':

    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
