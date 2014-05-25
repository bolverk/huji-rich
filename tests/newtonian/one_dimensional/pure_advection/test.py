#! /usr/bin/python

import math

def goodness_of_fit(a1, a2):

    import math

    diff2 = [(x-y)**2 for x,y in zip(a1,a2)]
    return math.sqrt(sum(diff2)/(max(a1)-min(a1))**2/len(a1))

def density_prof(x):

    if x<0.3 or x>0.7:
        return 1
    else:
        return 2

def pressure_prof(x):

    return 1

def velocity_prof(x):

    return 1

def plot_fluxes():

    import numpy
    import pylab

    rawd = numpy.loadtxt('fluxes.txt')
    radius = rawd[:,0]
    mass_flux = rawd[:,1]
    momentum_flux = rawd[:,2]
    energy_flux = rawd[:,3]
    
    pylab.subplot(311)
    pylab.plot(radius,mass_flux)

    pylab.subplot(312)
    pylab.plot(radius,momentum_flux)

    pylab.subplot(313)
    pylab.plot(radius,energy_flux)

    pylab.show()

def plot_hydro_vars():

    import numpy
    import pylab

    for int_name in ['pcm','plm_naive','plm']:
        x_list = numpy.loadtxt(
            int_name+'_cell_centres.txt')
        d_list = numpy.loadtxt(
            int_name+'_densities.txt')
        p_list = numpy.loadtxt(
            int_name+'_pressures.txt')
        v_list = numpy.loadtxt(
            int_name+'_velocities.txt')

        pylab.subplot(311)
        pylab.plot(x_list,d_list)
        pylab.ylabel('Density')
        
        pylab.subplot(312)
        pylab.plot(x_list,p_list)
        pylab.ylabel('Pressure')
        
        pylab.subplot(313)
        pylab.plot(x_list,v_list)
        pylab.ylabel('Velocity')
        pylab.xlabel('Radius')

    da = [density_prof(x) for x in x_list]
    pa = [pressure_prof(x) for x in x_list]
    va = [velocity_prof(x) for x in x_list]

    pylab.subplot(311)
    pylab.plot(x_list,da)

    pylab.subplot(312)
    pylab.plot(x_list,pa)

    pylab.subplot(313)
    pylab.plot(x_list,va)

    pylab.show()
        
def main():

    import os
    import numpy
    import h5py

    if False:
        return False

    if False:
        plot_fluxes()

    if False:
        plot_hydro_vars()

    pcm_data = h5py.File('pcm_final.h5')
    plm_data = h5py.File('plm_final.h5')
    x_list = pcm_data['grid']
    d_list_pcm = pcm_data['density']
    d_list_plm = plm_data['density']
    da = [density_prof(x) for x in x_list]
    gof_pcm = goodness_of_fit(d_list_pcm,da)
    gof_plm = goodness_of_fit(d_list_plm,da)

    f = open('gradesheet.txt','w')
    f.write(str(gof_pcm)+'\n')
    f.write(str(gof_plm)+'\n')
    f.close()

    return gof_plm<0.12 and \
        gof_pcm<0.22 and \
        gof_plm<gof_pcm

if __name__=='__main__':

    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')

