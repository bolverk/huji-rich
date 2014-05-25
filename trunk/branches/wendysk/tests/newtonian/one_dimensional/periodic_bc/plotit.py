import pylab
import numpy
import re
import sys
import Tkinter
import glob

pref = re.sub(r'/[^/]*$','',sys.argv[0])

n = len(glob.glob(pref+'/center_list.[0-9]+.txt'))

if len(sys.argv)>1:
    val = int(sys.argv[1])
else:
    val = n

x_list = numpy.loadtxt(pref+'/center_list.'+str(val)+'.txt')
d_list = numpy.loadtxt(pref+'/density_list.'+str(val)+'.txt')
p_list = numpy.loadtxt(pref+'/pressure_list.'+str(val)+'.txt')
v_list = numpy.loadtxt(pref+'/xvelocity_list.'+str(val)+'.txt')

pylab.subplot(311)
pylab.plot(x_list,d_list)
pylab.ylabel('Density')
    
pylab.subplot(312)
pylab.plot(x_list,p_list)
pylab.ylabel('Pressure')

pylab.subplot(313)
pylab.plot(x_list,v_list)
pylab.ylabel('Velocity')
pylab.xlabel('Distance')

pylab.show()
