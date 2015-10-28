#! /usr/bin/python

def my_sign(x):

    if x>0:
        return 1.0
    elif x<0:
        return -1.0
    elif x==0:
        return 0.0
    else:
        raise NameError('Error in my_sign: x is not numeric')

def afd_single_mode(x,t,k,omega,v,g0,dx,xi):

    from cmath import exp, sin, cos

    return (((-2)*1j)*dx*exp(1j*k*x)*g0*xi*((2*1j)*exp((1j*dx*omega*xi)/(2.*abs(\
v)))*abs(v) + v*xi*sin(dx*k))*(exp(1j*omega*t) - (1 - xi**2/4. + \
(xi**2*cos(2*dx*k))/4. - \
((1j)*xi*abs(v)*sin(dx*k))/v)**((t*abs(v))/(dx*xi))))/(v*(v*(-4 + \
4*exp((1j*dx*omega*xi)/abs(v)) + xi**2 - xi**2*cos(2*dx*k)) + \
(4*1j)*xi*abs(v)*sin(dx*k)))

def afd_standing_wave(x,t,k,omega,v,g0,dx,xi):

    return 0.25*(afd_single_mode(x,t,k,omega,v,g0,dx,xi)-
                 afd_single_mode(x,t,-k,omega,v,g0,dx,xi)-
                 afd_single_mode(x,t,k,-omega,v,g0,dx,xi)+
                 afd_single_mode(x,t,-k,-omega,v,g0,dx,xi))

def afd_standing_wave_2(x,t,k,omega,v,g0,dx,xi):

    import cmath

    return (dx*g0*xi*(-((cmath.exp(2*1j*k*x)*(cmath.exp(-(1j*omega*t)) - (1 + \
(-1 + cmath.exp(-(1j*dx*k)))*xi)**((t*v)/(dx*xi))))/(-1 + \
cmath.exp(-(1j*dx*omega*xi)/v) + xi - xi/cmath.exp(1j*dx*k))) + \
(cmath.exp(2*1j*k*x)*(cmath.exp(1j*omega*t) - (1 + (-1 + \
cmath.exp(-(1j*dx*k)))*xi)**((t*v)/(dx*xi))))/(-1 + \
cmath.exp((1j*dx*omega*xi)/v) + xi - xi/cmath.exp(1j*dx*k)) + \
(cmath.exp(-(1j*omega*t)) - (1 + (-1 + \
cmath.exp(1j*dx*k))*xi)**((t*v)/(dx*xi)))/(-1 + \
cmath.exp(-(1j*dx*omega*xi)/v) + xi - cmath.exp(1j*dx*k)*xi) + \
(-cmath.exp(1j*omega*t) + (1 + (-1 + \
cmath.exp(1j*dx*k))*xi)**((t*v)/(dx*xi)))/(-1 + \
cmath.exp((1j*dx*omega*xi)/v) + xi - \
cmath.exp(1j*dx*k)*xi)))/(4.*cmath.exp(1j*k*x)*v)

def afd_hydro_profiles(x_list,t,k,omega,g0,dx,xi,rho0,p0,gamma):

    import math

    c0 = math.sqrt(gamma*p0/rho0)
    j_p_list = [afd_standing_wave(x,t,k,omega,c0,g0,dx,xi).real for x in x_list]
    j_m_list = [afd_standing_wave(x,t,k,omega,-c0,g0,dx,xi).real for x in x_list]
    v_list = [0.5*(j_p+j_m) for j_p, j_m in zip(j_p_list,j_m_list)]
    rho_list = [rho0+0.5*(j_p-j_m)*rho0/c0
                for j_p, j_m
                in zip(j_p_list, j_m_list)]
    p_list = [p0+gamma*p0*(rho-rho0)/rho0 for rho in rho_list]

    return rho_list, p_list, v_list

def l1_fit_factor(a1, a2):

    import math

    diff2 = [abs(x-y) for x,y in zip(a1,a2)]
    return sum(diff2)/len(a1)

def main():

    import numpy
    import os
    import math
    import h5py

    h5f = h5py.File('final.h5')
    x_list = h5f['geometry']['x_coordinate']
    d_list = h5f['hydrodynamic']['density']
    p_list = h5f['hydrodynamic']['pressure']
    vx_list = h5f['hydrodynamic']['x_velocity']

    rho_0 = numpy.loadtxt('mean_density.txt')
    amp = numpy.loadtxt('amplitude.txt')
    k = 2*numpy.pi/numpy.loadtxt('wavelength.txt')
    vp = numpy.loadtxt('phase_velocity.txt')
    p_0 = numpy.loadtxt('mean_pressure.txt')
    gamma = numpy.loadtxt('adiabatic_index.txt')
    c_0 = numpy.sqrt(gamma*p_0/rho_0)
    t = numpy.loadtxt('time.txt')
    
    s_list = [p/d**gamma for d,p in zip(d_list,p_list)]
    j_p_list = [v+(rho-rho_0)*c_0/rho_0
                for rho,v in zip(d_list,vx_list)]
    j_m_list = [v-(rho-rho_0)*c_0/rho_0
                for rho,v in zip(d_list,vx_list)]

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
    s_exact = [p/d**gamma for d,p in zip(d_exact,p_exact)]
    j_p_exact = [v+(rho-rho_0)*c_0/rho_0
                 for rho,v in zip(d_exact,v_exact)]
    j_m_exact = [v-(rho-rho_0)*c_0/rho_0
                 for rho,v in zip(d_exact,v_exact)]

    dx = 1.0/math.sqrt(len(x_list))
    d_afd, p_afd, v_afd =  afd_hydro_profiles(x_list,t,k,k*vp,amp,dx,0.3,rho_0,p_0,gamma)
    s_afd = [p/d**gamma for d,p in zip(d_afd,p_afd)]
    j_p_afd = [v+(d-rho_0)*c_0/rho_0
               for d,v in zip(d_afd,v_afd)]
    j_m_afd = [v-(d-rho_0)*c_0/rho_0
               for d,v in zip(d_afd,v_afd)]

    if False:

        import pylab
        #pylab.subplot(311)
        #pylab.plot(x_list, s_list, '.')
        #pylab.plot(x_list, s_exact, '.')
        #pylab.plot(x_list, s_afd, '.')
        #pylab.subplot(312)
        #pylab.plot(x_list, j_p_list, '.')
        #pylab.plot(x_list, j_p_exact, '.')
        #pylab.plot(x_list, j_p_afd, '.')
        #pylab.subplot(313)
        #pylab.plot(x_list, j_m_list, '.')
        #pylab.plot(x_list, j_m_exact, '.')
        #pylab.plot(x_list, j_m_afd, '.')
        #pylab.show()

        pylab.subplot(311)
        pylab.plot(x_list, d_list, '.')
        pylab.plot(x_list, d_exact, '.')
        pylab.plot(x_list, d_afd, '.')
        pylab.subplot(312)
        pylab.plot(x_list, p_list, '.')
        pylab.plot(x_list, p_exact, '.')
        pylab.plot(x_list, p_afd, '.')
        pylab.subplot(313)
        pylab.plot(x_list, vx_list, '.')
        pylab.plot(x_list, v_exact, '.')
        pylab.plot(x_list, v_afd, '.')
        pylab.show()

    l1_density = l1_fit_factor(d_list, d_exact)
    l1_pressure = l1_fit_factor(p_list, p_exact)
    l1_velocity = l1_fit_factor(vx_list, v_exact)

    f = open('gradesheet.txt','w')
    f.write(str(l1_density)+'\n')
    f.write(str(l1_pressure)+'\n')
    f.write(str(l1_velocity)+'\n')
    f.close()
    
    return l1_density < 3.45e-5 and \
        l1_pressure < 3.44e-5 and \
        l1_velocity < 8.2e-5

if __name__=='__main__':
    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
