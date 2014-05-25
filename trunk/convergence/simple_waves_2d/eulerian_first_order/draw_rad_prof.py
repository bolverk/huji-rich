#! /usr/bin/python

import matplotlib
matplotlib.use('Qt4Agg')
import pylab
import numpy
import sys
import imp

# Simulation results

pref = sys.argv[0].replace('draw_rad_prof.py','')

rawd = numpy.loadtxt(pref+'x_prof_initial.txt')
xr0 = rawd[:,0]
dr0 = rawd[:,1]
pr0 = rawd[:,2]
vr0 = rawd[:,3]
cr0 = rawd[:,4]
time = numpy.loadtxt(pref+'time.txt')

xs0 = numpy.sort(xr0)
ids0 = numpy.argsort(xr0)
ds0 = [dr0[i] for i in ids0]
ps0 = [pr0[i] for i in ids0]
vs0 = [vr0[i] for i in ids0]
cs0 = [cr0[i] for i in ids0]

rawd = numpy.loadtxt(pref+'x_prof_final.txt')
xr = rawd[:,0]
dr = rawd[:,1]
pr = rawd[:,2]
vr = rawd[:,3]
cr = rawd[:,4]

xs = numpy.sort(xr)
ids = numpy.argsort(xr)
ds = [dr[i] for i in ids]
ps = [pr[i] for i in ids]
vs = [vr[i] for i in ids]
cs = [cr[i] for i in ids]

# Analytic results
xa = [x+time*(c+v) for x,v,c in zip(xs0,vs0,cs0)]

# Show data
pylab.subplot(311)
pylab.plot(xs,ds,xa,ds0,xs0,ds0)
pylab.xlabel('r')
pylab.ylabel('Density')

pylab.subplot(312)
pylab.plot(xs,ps,xa,ps0,xs0,ps0)
pylab.xlabel('r')
pylab.ylabel('Pressure')

pylab.subplot(313)
pylab.plot(xs,vs,xa,vs0,xs0,vs0)
pylab.xlabel('r')
pylab.ylabel('Velocity')
pylab.legend(('Numeric','Analytic','Initial'))

pylab.show()
