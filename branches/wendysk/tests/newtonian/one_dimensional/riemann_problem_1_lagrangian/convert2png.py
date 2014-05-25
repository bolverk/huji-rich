import pylab
import numpy
import re
import glob
import sys
import os

if __name__=='__main__':

    pref = re.sub(r'/[^/]*$','',sys.argv[0])

    for txtfile in glob.glob(pref+'/center_list.*.txt'):
        x_list = numpy.loadtxt(txtfile)
        d_list = numpy.loadtxt(txtfile.replace('center','density'))
        p_list = numpy.loadtxt(txtfile.replace('center','pressure'))
        v_list = numpy.loadtxt(txtfile.replace('center','xvelocity'))

        pylab.clf()

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

        pylab.savefig(txtfile.replace('center_list','snapshot').replace('txt','png'))
