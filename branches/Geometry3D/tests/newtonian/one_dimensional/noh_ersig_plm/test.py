#! /usr/bin/python

def l1_fit(a1,a2):

    abs_dif = [abs(x-y) for x,y in zip(a1,a2)]
    return sum(abs_dif)/len(abs_dif)

import numpy

def main(aux_flag=''):

    import h5py

    graphic_flag = False
    if aux_flag=='show':
        graphic_flag = True

    # Numeric profiles
    h5f = h5py.File('final.h5')
    center_list = h5f['grid']
    density_list = h5f['density']
    pressure_list = h5f['pressure']
    velocity_list = h5f['x_velocity']
    time = numpy.array(h5f['time'])[0]

    # Analyic profiles
    gamma = 5./3.
    upstream_speed = 1
    shock_speed = (gamma-1)*upstream_speed/2
    shock_position = shock_speed*time
    upstream_density = 1
    downstream_density = (gamma+1)*upstream_density/(gamma-1)
    upstream_velocity = -1
    downstream_velocity = 0
    upstream_pressure = 1e-6
    downstream_pressure = ((gamma+1)/2)*upstream_density*\
        (downstream_velocity-upstream_velocity)**2
    
    def density_profile(x):
        if x<shock_position:
            return downstream_density
        else:
            return upstream_density

    def pressure_profile(x):
        if x<shock_position:
            return downstream_pressure
        else:
            return upstream_pressure

    def velocity_profile(x):
        if x<shock_position:
            return downstream_velocity
        else:
            return upstream_velocity

    analytic_density = [density_profile(x)
                        for x in center_list]
    analytic_pressure = [pressure_profile(x)
                        for x in center_list]
    analytic_velocity = [velocity_profile(x)
                        for x in center_list]

    if graphic_flag:
        import pylab

        pylab.subplot(311)
        pylab.plot(center_list,density_list)
        pylab.plot(center_list,analytic_density)

        pylab.subplot(312)
        pylab.plot(center_list,pressure_list)
        pylab.plot(center_list,analytic_pressure)

        pylab.subplot(313)
        pylab.plot(center_list,velocity_list)
        pylab.plot(center_list,analytic_velocity)

        pylab.show()

    l1_density = l1_fit(analytic_density,density_list)
    l1_pressure = l1_fit(analytic_pressure,pressure_list)
    l1_velocity = l1_fit(analytic_velocity,velocity_list)

    cond1 = l1_density<2e-2
    cond2 = l1_pressure<1e-2
    cond3 = l1_velocity<1e-2

    f = open('gradesheet.txt','w')
    f.write(str(l1_density)+'\n')
    f.write(str(l1_pressure)+'\n')
    f.write(str(l1_velocity)+'\n')
    f.close()

    return cond1 and cond2 and cond3

if __name__=='__main__':
    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')

