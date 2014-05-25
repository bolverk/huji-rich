#! /usr/bin/python

import matplotlib
matplotlib.use('Qt4Agg')
import pylab
import numpy
import sys
import imp
enrs = imp.load_source('enrs','../analytic/enrs.py')

pref = sys.argv[0].replace('draw_profiles.py','')

left = enrs.Primitive(1,10,0);
right = enrs.Primitive(1,1,0);
prof = enrs.RiemannProfile(left,right,5./3.)
t = numpy.loadtxt(pref+'time.txt')
offset = 0.5

x = numpy.loadtxt(pref+'cell_centres.txt')
d = numpy.loadtxt(pref+'densities.txt')
p = numpy.loadtxt(pref+'pressures.txt')
v = numpy.loadtxt(pref+'velocities.txt')
da = [prof.CalcPrim((i-offset)/t).Density for i in x]
pa = [prof.CalcPrim((i-offset)/t).Pressure for i in x]
va = [prof.CalcPrim((i-offset)/t).Velocity for i in x]

pylab.subplot(311)
pylab.plot(x,d,x,da)
pylab.xlabel('x')
pylab.ylabel('Density')

pylab.subplot(312)
pylab.plot(x,p,x,pa)
pylab.xlabel('x')
pylab.ylabel('Pressure')

pylab.subplot(313)
pylab.plot(x,v,x,va)
pylab.xlabel('x')
pylab.ylabel('Velocity')

pylab.show()
