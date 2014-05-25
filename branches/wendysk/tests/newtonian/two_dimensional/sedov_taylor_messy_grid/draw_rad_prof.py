#! /usr/bin/python

import matplotlib
matplotlib.use('Qt4Agg')
import pylab
import numpy
import sys
import imp
sedov_taylor = imp.load_source('sedov_taylor','../analytic/sedov_taylor.py')

# Simulation results

pref = sys.argv[0].replace('draw_rad_prof.py','')

rawd = numpy.loadtxt(pref+'rad_prof.txt')
rr = rawd[:,0]
dr = rawd[:,1]
pr = rawd[:,2]
vr = rawd[:,3]

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

# Show data
pylab.subplot(311)
pylab.plot(rs,ds,ra,da)
pylab.xlabel('r')
pylab.ylabel('Density')

pylab.subplot(312)
pylab.plot(rs,ps,ra,pa)
pylab.xlabel('r')
pylab.ylabel('Pressure')

pylab.subplot(313)
pylab.plot(rs,vs,ra,va)
pylab.xlabel('r')
pylab.ylabel('Velocity')

pylab.show()
