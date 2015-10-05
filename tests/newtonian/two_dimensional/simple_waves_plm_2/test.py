#! /usr/bin/python

def goodness_of_fit(a1, a2):

    import math

    diff2 = [(x-y)**2 for x,y in zip(a1,a2)]
    return math.sqrt(sum(diff2)/(max(a1)-min(a1))/len(a1))

def main():

    import numpy
    import h5py

    # Simulation results

    g = numpy.loadtxt('adiabatic_index.txt')
    initial_data = h5py.File('initial.h5')
    final_data = h5py.File('final.h5')
    xr0 = initial_data['geometry']['x_coordinate']
    dr0 = initial_data['hydrodynamic']['density']
    pr0 = initial_data['hydrodynamic']['pressure']
    vr0 = initial_data['hydrodynamic']['x_velocity']
    cr0 = [numpy.sqrt(g*p/d) for d,p in zip(dr0,pr0)]
    time = numpy.array(final_data['time'])[0]

    xs0 = numpy.sort(xr0)
    ids0 = numpy.argsort(xr0)
    ds0 = [dr0[i] for i in ids0]
    ps0 = [pr0[i] for i in ids0]
    vs0 = [vr0[i] for i in ids0]
    cs0 = [cr0[i] for i in ids0]

    xr = final_data['geometry']['x_coordinate']
    dr = final_data['hydrodynamic']['density']
    pr = final_data['hydrodynamic']['pressure']
    vr = final_data['hydrodynamic']['x_velocity']
    cr = [numpy.sqrt(g*p/d) for d,p in zip(dr,pr)]

    xs = numpy.sort(xr)
    ids = numpy.argsort(xr)
    ds = [dr[i] for i in ids]
    ps = [pr[i] for i in ids]
    vs = [vr[i] for i in ids]
    cs = [cr[i] for i in ids]
    
    # Analytic results
    xa = [x+time*(c+v) for x,v,c in zip(xs0,vs0,cs0)]

    # Prepare for interpolation
    x_inside = [z for z in xa if z<numpy.max(xs) and 
                z>numpy.min(xs)]
    d_analytic = [ds0[i] for i in range(len(xa)) 
                  if xa[i]<numpy.max(xs) and
                  xa[i]>numpy.min(xs)]
    p_analytic = [ps0[i] for i in range(len(xa)) 
                  if xa[i]<numpy.max(xs) and
                  xa[i]>numpy.min(xs)]
    v_analytic = [vs0[i] for i in range(len(xa)) 
                  if xa[i]<numpy.max(xs) and
                  xa[i]>numpy.min(xs)]

    #d_analytic = numpy.interp(x_inside,xa,ds0)
    d_numeric = numpy.interp(x_inside,xs,ds)
    #p_analytic = numpy.interp(x_inside,xa,ps0)
    p_numeric = numpy.interp(x_inside,xs,ps)
    #v_analytic = numpy.interp(x_inside,xa,vs0)
    v_numeric = numpy.interp(x_inside,xs,vs)

    gof1 = goodness_of_fit(d_numeric,d_analytic)
    gof2 = goodness_of_fit(p_numeric,p_analytic)
    gof3 = goodness_of_fit(v_numeric,v_analytic)

    f = open('gradesheet.txt','w')
    f.write(str(gof1)+'\n')
    f.write(str(gof2)+'\n')
    f.write(str(gof3)+'\n')
    f.close()

    return gof1<0.09 and \
        gof2<0.26 and \
        gof3<0.11

if __name__=='__main__':

    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')

