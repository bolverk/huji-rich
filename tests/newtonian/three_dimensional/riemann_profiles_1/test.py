#! /usr/bin/python

def goodness_of_fit(a1, a2):

    import math

    diff2 = [(x-y)**2 for x,y in zip(a1,a2)]
    return math.sqrt(sum(diff2)/(max(a1)-min(a1))/len(a1))

def consolidate_single(fname):

    import h5py
    import numpy

    res = {}
    with h5py.File(fname, 'r') as f :
        res['x_coordinate'] = numpy.array(f['geometry']['x_coordinate'])
        res['density'] = numpy.array(f['hydrodynamic']['density'])
        res['pressure'] = numpy.array(f['hydrodynamic']['pressure'])
        res['velocity'] = numpy.array(f['hydrodynamic']['x_velocity'])
        res['time'] = numpy.array(f['time'])[0]
    return res

def main():

    import numpy
    import os
    import h5py
    import glob
    from importlib.machinery import SourceFileLoader
    enrs = SourceFileLoader(
        'enrs',
        os.environ['RICH_ROOT']+'/analytic/enrs.py').load_module()

    np = len(glob.glob('process_*_final.h5'))
    if np>0:
        x = []
        d = []
        p = []
        v = []
        for fname in glob.glob('process_*_final.h5'):
            f = h5py.File(fname)
            x.extend(f['x_coordinate'])
            d.extend(f['density'])
            p.extend(f['pressure'])
            v.extend(f['x_velocity'])
            t = f['time']
    else:
        h5f = h5py.File('final.h5')
        x_list = h5f['X']
        d_list = h5f['Density']                     
        p_list = h5f['Pressure']
        v_list = h5f['Vx']
        t = h5f['Time']

    combined = sorted([(x,d,p,v) for x,d,p,v in zip(x_list, d_list, p_list, v_list)], key=lambda tup:tup[0])

    xs = [tup[0] for tup in combined]
    ds = [tup[1] for tup in combined]
    ps = [tup[2] for tup in combined]
    vs = [tup[3] for tup in combined]

    left = enrs.Primitive(1,2,0)
    right = enrs.Primitive(1,1,0)
    prof = enrs.RiemannProfile(left,right,5./3.)
    offset = 0;
    da = [prof.CalcPrim((i-offset)/t).Density for i in xs]
    pa = [prof.CalcPrim((i-offset)/t).Pressure for i in xs]
    va = [prof.CalcPrim((i-offset)/t).Velocity for i in xs]

    if False:
        import pylab
        pylab.subplot(311)
        pylab.plot(xs,ds,'.')
        pylab.plot(xs,da)
        pylab.subplot(312)
        pylab.plot(xs,ps,'.')
        pylab.plot(xs,pa)
        pylab.subplot(313)
        pylab.plot(xs,vs,'.')
        pylab.plot(xs,va)
        pylab.show()

    gof1 = goodness_of_fit(ds,da)
    gof2 = goodness_of_fit(ps,pa)
    gof3 = goodness_of_fit(vs,va)

    with open('gradesheet.txt','w') as f:
        f.write(str(gof1)+'\n')
        f.write(str(gof2)+'\n')
        f.write(str(gof3)+'\n')
    
    return gof1<0.12 and gof2<0.11 and gof3<0.13

if __name__=='__main__':

    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
