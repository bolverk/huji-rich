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
    import glob
    enrs = imp.load_source('enrs',os.environ['RICH_ROOT']+'/analytic/enrs.py')

    ns = len(glob.glob('process_*_final.h5'))
    if ns>0:
        x = []
        d = []
        p = []
        v = []
        for fname in glob.glob('process_*_final.h5'):
            f = h5py.File(fname,'r+')
            x.extend(f['x_coordinate'])
            d.extend(f['density'])
            p.extend(f['pressure'])
            v.extend(f['x_velocity'])
    else:
        h5f = h5py.File('final.h5')
        x = h5f['geometry']['x_coordinate']
        d = h5f['hydrodynamic']['density']
        p = h5f['hydrodynamic']['pressure']
        v = h5f['hydrodynamic']['x_velocity']

    xs = numpy.sort(x)
    ids = numpy.argsort(x)
    ds = [d[i] for i in ids]
    ps = [p[i] for i in ids]
    vs = [v[i] for i in ids]

    left = enrs.Primitive(1,2,0)
    right = enrs.Primitive(1,1,0)
    prof = enrs.RiemannProfile(left,right,5./3.)
    t = numpy.loadtxt('time.txt')
    offset = 0.5;
    da = [prof.CalcPrim((i-offset)/t).Density for i in xs]
    pa = [prof.CalcPrim((i-offset)/t).Pressure for i in xs]
    va = [prof.CalcPrim((i-offset)/t).Velocity for i in xs]

    gof1 = goodness_of_fit(ds,da)
    gof2 = goodness_of_fit(ps,pa)
    gof3 = goodness_of_fit(vs,va)

    f = open('gradesheet.txt','w')
    f.write(str(gof1)+'\n')
    f.write(str(gof2)+'\n')
    f.write(str(gof3)+'\n')
    f.close()
    
    return gof1<0.067 and gof2<0.053 and gof3<0.062

if __name__=='__main__':

    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
