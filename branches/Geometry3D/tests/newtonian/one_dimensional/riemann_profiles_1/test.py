#! /usr/bin/python

def goodness_of_fit(a1, a2):

    import math

    diff2 = [(x-y)**2 for x,y in zip(a1,a2)]
    return math.sqrt(sum(diff2)/(max(a1)-min(a1))**2/len(a1))

def main():

    import os
    import imp
    enrs = imp.load_source('enrs',os.environ['RICH_ROOT']+'/analytic/enrs.py')

    import numpy
    import h5py

    left = enrs.Primitive(1,10,0);
    right = enrs.Primitive(1,1,0);
    prof = enrs.RiemannProfile(left,right,5./3.)
    offset = 0.5

    h5file = h5py.File('final.h5','r')
    t = h5file['time']
    x = h5file['grid']
    d = h5file['density']
    p = h5file['pressure']
    v = h5file['x_velocity']
    da = [prof.CalcPrim((i-offset)/t).Density for i in x]
    pa = [prof.CalcPrim((i-offset)/t).Pressure for i in x]
    va = [prof.CalcPrim((i-offset)/t).Velocity for i in x]

    gof1 = goodness_of_fit(d,da)
    gof2 = goodness_of_fit(p,pa)
    gof3 = goodness_of_fit(v,va)

    f = open('gradesheet.txt','w')
    f.write(str(gof1)+'\n')
    f.write(str(gof2)+'\n')
    f.write(str(gof3)+'\n')
    f.close()

    return gof1<0.14 and gof2<0.06 and gof3 < 0.13

if __name__=='__main__':
    
    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')

