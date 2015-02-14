#! /usr/bin/python

import math

def goodness_of_fit(a1, a2):

    import math

    diff2 = [(x-y)**2 for x,y in zip(a1,a2)]
    return math.sqrt(sum(diff2)/(max(a1)-min(a1))**2/len(a1))

def density_prof(x):

    return 1+1e-6*math.sin(2*math.pi*x)

def pressure_prof(x):

    return 3./5.+1e-6*math.sin(2*math.pi*x)

def velocity_prof(x):

    return 1e-6*math.sin(2*math.pi*x)

def main():

    import os
    import numpy
    import h5py

    h5f = h5py.File('final.h5')
    x_list = h5f['grid']
    d_list = h5f['density']
    p_list = h5f['pressure']
    v_list = h5f['x_velocity']
    da = [density_prof(x) for x in x_list]
    pa = [pressure_prof(x) for x in x_list]
    va = [velocity_prof(x) for x in x_list]

    if False:
        import pylab

        pylab.subplot(311)
        pylab.plot(x_list,d_list,x_list,da)

        pylab.subplot(312)
        pylab.plot(x_list,p_list,x_list,pa)

        pylab.subplot(313)
        pylab.plot(x_list,v_list,x_list,va)

        pylab.show()

    gof1 = goodness_of_fit(d_list,da)
    gof2 = goodness_of_fit(p_list,pa)
    gof3 = goodness_of_fit(v_list,va)

    f = open('gradesheet.txt','w')
    f.write(str(gof1)+'\n')
    f.write(str(gof2)+'\n')
    f.write(str(gof3)+'\n')
    f.close()

    return gof1<0.06 and gof2<0.02 and gof3 < 0.04

if __name__=='__main__':

    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')

