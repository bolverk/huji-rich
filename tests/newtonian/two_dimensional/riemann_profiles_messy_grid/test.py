#! /usr/bin/python

def goodness_of_fit(a1, a2):

    import math

    diff2 = [(x-y)**2 for x,y in zip(a1,a2)]
    return math.sqrt(sum(diff2)/(max(a1)-min(a1))/len(a1))

def main():

    import numpy
    import os
    import imp
    import h5py
    enrs = imp.load_source('enrs',os.environ['RICH_ROOT']+'/analytic/enrs.py')

    h5f = h5py.File('final.h5')
    x = h5f['x_coordinate']
    d = h5f['density']
    p = h5f['pressure']
    v = h5f['x_velocity']
    t = h5f['time']

    xs = numpy.sort(x)
    ids = numpy.argsort(x)
    ds = [d[i] for i in ids]
    ps = [p[i] for i in ids]
    vs = [v[i] for i in ids]

    left = enrs.Primitive(1,2,0)
    right = enrs.Primitive(1,1,0)
    prof = enrs.RiemannProfile(left,right,5./3.)
    offset = 0.5;
    da = [prof.CalcPrim((i-offset)/t).Density for i in xs]
    pa = [prof.CalcPrim((i-offset)/t).Pressure for i in xs]
    va = [prof.CalcPrim((i-offset)/t).Velocity for i in xs]

    gof1 = goodness_of_fit(ds,da)
    gof2 = goodness_of_fit(ps,pa)
    gof3 = goodness_of_fit(vs,va)
    
    return gof1<0.1 and gof2<0.11 and gof3<0.11

if __name__=='__main__':

    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
