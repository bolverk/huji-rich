#! /usr/bin/python

import matplotlib
matplotlib.use('Qt4Agg')
import pylab
import numpy
import sys
import os
import imp
enrs = imp.load_source('enrs','../analytic/enrs.py')

pref = sys.argv[0].replace('draw_prof_1d.py','')

rawd = numpy.loadtxt(pref+'prof1d.txt')
x = rawd[:,0]
d = rawd[:,1]
p = rawd[:,2]
v = rawd[:,3]

xs = numpy.sort(x)
ids = numpy.argsort(x)
ds = [d[i] for i in ids]
ps = [p[i] for i in ids]
vs = [v[i] for i in ids]

left = enrs.Primitive(1,2,0)
right = enrs.Primitive(1,1,0)
prof = enrs.RiemannProfile(left,right,5./3.)
t = numpy.loadtxt(pref+'time.txt')
offset = 0.5;
da = [prof.CalcPrim((i-offset)/t).Density for i in xs]
pa = [prof.CalcPrim((i-offset)/t).Pressure for i in xs]
va = [prof.CalcPrim((i-offset)/t).Velocity for i in xs]

pylab.subplot(311)
pylab.plot(xs,ds,xs,da)
pylab.xlabel('x')
pylab.ylabel('Density')
pylab.title('t = '+str(t))

pylab.subplot(312)
pylab.plot(xs,ps,xs,pa)
pylab.xlabel('x')
pylab.ylabel('Pressure')

pylab.subplot(313)
pylab.plot(xs,vs,xs,va)
pylab.xlabel('x')
pylab.ylabel('Velocity')

pylab.show()
