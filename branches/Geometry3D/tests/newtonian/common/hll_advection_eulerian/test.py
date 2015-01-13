class Primitive:
    """
    Primitive hydrodynamic variables
    """

    def __init__(self,density,pressure,xvelocity,
                 yvelocity,energy,soundspeed):

        self.density = density
        self.pressure = pressure
        self.xvelocity = xvelocity
        self.yvelocity = yvelocity
        self.energy = energy
        self.soundspeed = soundspeed

class Flux:
    """
    Hydrodynamic flux
    """

    def __init__(self,mass,xmomentum,ymomentum,energy):

        self.mass = mass
        self.xmomentum = xmomentum
        self.ymomentum = ymomentum
        self.energy = energy

def read_primitive(fname):

    import numpy

    rawd = numpy.loadtxt(fname)
    return Primitive(rawd[0],
                     rawd[1],
                     rawd[4],
                     rawd[5],
                     rawd[2],
                     rawd[3])

def read_flux(fname):

    import numpy

    rawd = numpy.loadtxt(fname)
    return Flux(rawd[0],
                rawd[1],
                rawd[2],
                rawd[3])

def approx_compare(val1,val2,thres):

    if 0==val1 and 0==val2:
        return True
    else:
        return abs(val1-val2)/abs(val1+val2)<thres

def main():

    left = read_primitive('left.txt')
    right = read_primitive('right.txt')
    numeric_flux = read_flux('flux.txt')

    exact_mass_flux = left.xvelocity*left.density
    exact_xmomentum_flux = left.pressure+\
        left.density*(left.xvelocity)**2
    exact_ymomentum_flux = left.density*\
        left.xvelocity*left.yvelocity
    exact_energy_flux = left.xvelocity*\
        (left.pressure+\
             left.density*left.energy+\
             0.5*left.density*left.xvelocity**2)

    f = open('gradesheet.txt','w')
    f.write('Component Numeric Exact\n')
    f.write('Mass '+str(numeric_flux.mass)+' '+str(exact_mass_flux)+'\n')
    f.write('xMomentum '+str(numeric_flux.xmomentum)+' '+str(exact_xmomentum_flux)+'\n')
    f.write('yMomentum '+str(numeric_flux.ymomentum)+' '+str(exact_ymomentum_flux)+'\n')
    f.write('Energy '+str(numeric_flux.energy)+' '+str(exact_energy_flux)+'\n')

    flag1 = approx_compare(numeric_flux.mass,exact_mass_flux,0.01)
    flag2 = approx_compare(numeric_flux.xmomentum,exact_xmomentum_flux,0.01)
    flag3 = approx_compare(numeric_flux.ymomentum,exact_ymomentum_flux,0.01)
    flag4 =  approx_compare(numeric_flux.energy,exact_energy_flux,0.01)
    return flag1 and flag2 and flag3 and flag4

if __name__=='__main__':

    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
