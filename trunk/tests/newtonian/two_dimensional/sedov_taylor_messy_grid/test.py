#! /usr/bin/python

def goodness_of_fit(a1, a2):

    import math

    diff2 = [(x-y)**2 for x,y in zip(a1,a2)]
    return math.sqrt(sum(diff2)/(max(a1)-min(a1))/len(a1))

def main():

    import numpy
    import imp
    import h5py
    import os
    import glob

    sedov_taylor = imp.load_source(\
        'sedov_taylor',\
            os.environ['RICH_ROOT']+'/analytic/sedov_taylor.py')

    # Simulation results

    np = len(glob.glob('process_*_final.h5'))

    if np>0:
        x_numeric = []
        y_numeric = []
        vx_numeric = []
        vy_numeric = []
        dr = []
        pr = []
        for idx in range(np):
            f = h5py.File('process_'+str(idx)+'_final.h5')
            x_numeric.extend(f['x_coordinate'])
            y_numeric.extend(f['y_coordinate'])
            vx_numeric.extend(f['x_velocity'])
            vy_numeric.extend(f['y_velocity'])
            dr.extend(f['density'])
            pr.extend(f['pressure'])
    else:
        h5f = h5py.File('final.h5')
        x_numeric = h5f['x_coordinate']
        y_numeric = h5f['y_coordinate']
        vx_numeric = h5f['x_velocity']
        vy_numeric = h5f['y_velocity']
        dr = h5f['density']
        pr = h5f['pressure']
    rr = [numpy.sqrt(x*x+y*y) for x,y in zip(x_numeric, y_numeric)]
    vr = [numpy.sqrt(vx*vx+vy*vy) for vx, vy
          in zip(vx_numeric, vy_numeric)]

    rs = numpy.sort(rr)
    ids = numpy.argsort(rr)
    ds = [dr[i] for i in ids]
    ps = [pr[i] for i in ids]
    vs = [vr[i] for i in ids]

    # Analytic results
    R = rs[numpy.argmax(ps)]
    vatsf = vs[numpy.argmax(ps)]
    patsf = ps[numpy.argmax(ps)]
    datsf = ds[numpy.argmax(ps)]
    g = 5./3. 
    w = 0
    n = 2
    v = numpy.linspace(1e-6+1./g,2/(g+1),num=1000)
    ra = [R*sedov_taylor.vtoz(i,w,g,n) for i in v]
    da = [datsf*sedov_taylor.vtod(i,w,g,n)/((g+1)/(g-1)) for i in v]
    pa = [patsf*sedov_taylor.vtop(i,w,g,n)/(2/(g+1)) for i in v]
    va = [vatsf*sedov_taylor.vtoz(i,w,g,n)*i/(2/(g+1)) for i in v]

    df = numpy.interp(ra,rs,ds)
    pf = numpy.interp(ra,rs,ps)
    vf = numpy.interp(ra,rs,vs)

    gof1 = goodness_of_fit(df,da)
    gof2 = goodness_of_fit(pf,pa)
    gof3 = goodness_of_fit(vf,va)

    f = open('gradesheet.txt','w')
    f.write(str(gof1)+'\n')
    f.write(str(gof2)+'\n')
    f.write(str(gof3)+'\n')
    f.close()

    return gof1<0.09 and \
        gof2<0.8 and gof3<0.11

import sys
if __name__=='__main__':

    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
        
