#! /usr/bin/python

import matplotlib
matplotlib.use('Qt4Agg')
import pylab
import numpy
import sys

aux = sys.argv[0]
pref = aux[:aux.find('show_mesh_points.py')]

rawd = numpy.loadtxt(pref+'mesh_generating_points.txt')
pylab.plot(rawd[:,0], rawd[:,1],'x')

rawd2 = numpy.loadtxt(pref+'edges.txt')
for i in range(len(rawd2[:,0])):
    pylab.plot((rawd2[i,0],rawd2[i,2]),(rawd2[i,1],rawd2[i,3]),'r')
pylab.axis('equal')
pylab.xlim((-0.05, 1.05))
pylab.ylim((-0.05, 1.05))
pylab.show()

