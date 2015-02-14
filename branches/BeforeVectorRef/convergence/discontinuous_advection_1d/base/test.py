#! /usr/bin/python

import math

def chi2_test(a1,a2):

    import math

    diff2 = [(x-y)**2 for x,y in zip(a1,a2)]
    return math.sqrt(sum(diff2)/len(a1))

def L1_test(a1,a2):

    import math

    abs_diff = [abs(x-y) for x,y in zip(a1,a2)]
    return sum(abs_diff)/len(a1)

def goodness_of_fit(a1, a2,test_name='L1'):

    import math

    if 'chi2'==test_name:
        return chi2_test(a1,a2)
    elif 'L1'==test_name:
        return L1_test(a1,a2)
    else:
        raise NameError("I don't know test "+
                        test_name)

def density_init_prof(x):

    if x<=0.3 or x>=0.7:
        return 1
    else:
        return 2

def density_prof(x,t):

    xeq = x-t
    if xeq < 0:
        xeq = xeq + 1
    return density_init_prof(xeq)

def pressure_prof(x):

    return 1

def velocity_prof(x):

    return 1

def main():

    import os

    import numpy

    time = numpy.loadtxt('time.txt')
    x_list = numpy.loadtxt('cell_centres.txt')
    d_list = numpy.loadtxt('densities.txt')
    p_list = numpy.loadtxt('pressures.txt')
    v_list = numpy.loadtxt('velocities.txt')
    da = [density_prof(x,time) for x in x_list]
    pa = [pressure_prof(x) for x in x_list]
    va = [velocity_prof(x) for x in x_list]

    if False:
        import pylab

        pylab.subplot(311)
        pylab.plot(x_list,d_list,x_list,da)
        pylab.ylabel('Density')

        pylab.subplot(312)
        pylab.plot(x_list,p_list,x_list,pa)
        pylab.ylabel('Pressure')

        pylab.subplot(313)
        pylab.plot(x_list,v_list,x_list,va)
        pylab.xlabel('Distance')
        pylab.ylabel('Velocity')

        pylab.show()

    test_name = 'L1'
    gof1 = goodness_of_fit(d_list,da,test_name)
    gof2 = goodness_of_fit(p_list,pa,test_name)
    gof3 = goodness_of_fit(v_list,va,test_name)

    f = open('gradesheet.txt','w')
    f.write(str(gof1)+'\n')
    f.write(str(gof2)+'\n')
    f.write(str(gof3)+'\n')
    f.close()

    return gof1<0.06 and gof2<0.02 and gof3 < 0.04

import sys

if __name__=='__main__':

    print main()
