#! /usr/bin/python

import matplotlib
matplotlib.use('Qt4Agg')
import pylab
import numpy
import sys

rfname = sys.argv[1]
vfname = sys.argv[2]

r = numpy.loadtxt(rfname)
v = numpy.loadtxt(vfname)

pylab.plot(r,v)
pylab.show()
