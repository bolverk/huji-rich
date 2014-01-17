#! /usr/bin/python

import math

def chi_squared(a1, a2):

    import math

    diff2 = [(x-y)**2 for x,y in zip(a1,a2)]
    return math.sqrt(sum(diff2)/(max(a1)-min(a1))**2/len(a1))

def L1_error_norm(a1,a2):

    import math

    abs_diff = [abs(x-y) for x,y in zip(a1,a2)]
    return sum(abs_diff)/len(a1)

def error_norm(a1,a2,test_name='L1'):

    if 'L1'==test_name:
        return L1_error_norm(a1,a2)
    elif 'chi2'==test_name:
        return chi_squared(a1,a2)
    else:
        raise NameError('I do not know test '+test_name)

def density_prof(x,t):

    return 1+1e-6*math.sin(2*math.pi*(x-t))

def pressure_prof(x,t):

    return 3./5.+1e-6*math.sin(2*math.pi*(x-t))

def velocity_prof(x,t):

    return 1e-6*math.sin(2*math.pi*(x-t))

def plot_primitives(x_list,
                    d_list,
                    p_list,
                    v_list,
                    d_exact,
                    p_exact,
                    v_exact):
    import pylab

    pylab.subplot(311)
    pylab.plot(x_list,d_list,
               x_list,d_exact)

    pylab.subplot(312)
    pylab.plot(x_list,p_list,
               x_list,p_exact)

    pylab.subplot(313)
    pylab.plot(x_list,v_list,
               x_list,v_exact)

    pylab.show()

def plot_primitives_dfr(x_list,
                        d_list,
                        p_list,
                        v_list,
                        d_exact,
                        p_exact,
                        v_exact):
    import pylab

    pylab.subplot(311)
    pylab.plot(x_list,d_list-d_exact)

    pylab.subplot(312)
    pylab.plot(x_list,p_list-p_exact)

    pylab.subplot(313)
    pylab.plot(x_list,v_list-v_exact)

    pylab.show()

def calc_entropy(d,p,g):

    return p/d**g

def calc_riemann_invariant(d,p,v,g,sgn):

    if abs(sgn)<>1:
        raise NameError('Sing must be either 1 or -1')

    c = (g*p/d)**0.5
    return 2*c/(g-1)+sgn*v

def calc_invariants(d_list,
                    p_list,
                    v_list,
                    g):
    
    entropy = [calc_entropy(d,p,g) for d,p in zip(d_list,p_list)]
    jplus = [calc_riemann_invariant(d,p,v,g,1) for d,p,v in zip(d_list,p_list,v_list)]
    jminus = [calc_riemann_invariant(d,p,v,g,-1) for d,p,v in zip(d_list,p_list,v_list)]
    
    return entropy, jplus, jminus

def plot_invariants(x_list,
                    d_numeric,
                    p_numeric,
                    v_numeric,
                    d_exact,
                    p_exact,
                    v_exact,
                    g):

    s_n, jp_n, jm_n = calc_invariants(d_numeric,
                                      p_numeric,
                                      v_numeric, g)

    s_e, jp_e, jm_e = calc_invariants(d_exact,
                                      p_exact,
                                      v_exact, g)

    import pylab

    pylab.subplot(311)
    pylab.plot(x_list, s_n,
               x_list, s_e)

    pylab.subplot(312)
    pylab.plot(x_list, jp_n,
               x_list, jp_e)

    pylab.subplot(313)
    pylab.plot(x_list, jm_n,
               x_list, jm_e)

    pylab.show()

def main():

    import os
    import h5py
    import numpy

    h5f = h5py.File('final.h5')
    t= numpy.array(h5f['time'])[0]
    x_list = h5f['grid']
    d_list = h5f['density']
    p_list = h5f['pressure']
    v_list = h5f['x_velocity']
    da = [density_prof(x,t) for x in x_list]
    pa = [pressure_prof(x,t) for x in x_list]
    va = [velocity_prof(x,t) for x in x_list]

    if False:
        plot_primitives(x_list,
                        d_list,
                        p_list,
                        v_list,
                        da,pa,va)

    if False:
        plot_primitives_dfr(x_list,
                            d_list,
                            p_list,
                            v_list,
                            da,pa,va)

    if False:
        plot_invariants(x_list,
                        d_list,
                        p_list,
                        v_list,
                        da,pa,va,5./3.)

    test_name = 'L1'
    gof1 = error_norm(d_list,da,test_name)
    gof2 = error_norm(p_list,pa,test_name)
    gof3 = error_norm(v_list,va,test_name)

    f = open('gradesheet.txt','w')
    f.write(str(gof1)+'\n')
    f.write(str(gof2)+'\n')
    f.write(str(gof3)+'\n')
    f.close()

    return gof1<0.06 and gof2<0.02 and gof3 < 0.04

import sys

if __name__=='__main__':

    print main()
