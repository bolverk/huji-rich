import pylab
import numpy
import re
import sys
import Tkinter
import glob

pref = re.sub(r'/[^/]*$','',sys.argv[0])

n = len(glob.glob(pref+'/center_list.*.txt'))

if len(sys.argv)>1:
    val = int(sys.argv[1])
else:
    val = n

edge_list = numpy.loadtxt(pref+'/edge_list.'+str(val)+'.txt')
mass_flux_list = numpy.loadtxt(pref+'/mass_flux_list.'+str(val)+'.txt')
xmom_flux_list = numpy.loadtxt(pref+'/x_mom_flux_list.'+str(val)+'.txt')
ymom_flux_list = numpy.loadtxt(pref+'/y_mom_flux_list.'+str(val)+'.txt')
enr_flux_list = numpy.loadtxt(pref+'/enr_flux_list.'+str(val)+'.txt')

pylab.subplot(221)
pylab.plot(edge_list,
           mass_flux_list)
pylab.ylabel('Mass flux')
    
pylab.subplot(222)
pylab.plot(edge_list,
           xmom_flux_list)
pylab.ylabel('x momentum flux')

pylab.subplot(223)
pylab.plot(edge_list,
           ymom_flux_list)
pylab.ylabel('y momentum flux')

pylab.subplot(224)
pylab.plot(edge_list,
           enr_flux_list)
pylab.ylabel('energy flux')

pylab.show()
