def calc_error_norm(a1,a2):

    abs_dif = [abs(x-y) for x,y in zip(a1,a2)]
    return sum(abs_dif)/len(a1)

def read_x_prof(fname):

    import numpy
    rawd = numpy.loadtxt(fname)
    return {'coordinate': rawd[:,0],
            'density': rawd[:,1],
            'pressure': rawd[:,2],
            'velocity': rawd[:,3],
            'sound speed': rawd[:,4]}

def calc_sound_speed(d,p,g):

    import math

    return math.sqrt(g*p/d)

class SineWave:
    """
    Sine wave spatial profile
    """

    def __init__(self,
                 amplitude,
                 wavelength,
                 phase,
                 offset):

        import math

        self.amp = amplitude
        self.k = 2*math.pi/float(wavelength)
        self.ph = phase
        self.offset = offset

    def eval(self, x):

        import math

        amp = self.amp
        k = self.k
        ph = self.ph
        offset = self.offset
        return amp*math.sin(k*x+ph)+offset

class SpatProf:

    """
    Acoustic spatial profiles for this test
    """

    def __init__(self, d0, p0, g,
                 dd, l):

        c0 = calc_sound_speed(d0,p0,g)
        self.c0 = c0
        self.density = SineWave(dd,l,0,d0)
        self.pressure = SineWave(dd*c0**2,
                                 l,0,p0)
        self.velocity = SineWave(dd*c0/d0,
                                 l,0,0)

def main():

    graphic_flag = False
    if graphic_flag:
        import pylab

    import numpy
    import h5py

    initial = h5py.File('initial.h5')
    final_data = h5py.File('final.h5')

    time = numpy.array(final_data['time'])[0]
    spat_prof = SpatProf(numpy.loadtxt('ambient_density.txt'),
                         numpy.loadtxt('ambient_pressure.txt'),
                         numpy.loadtxt('adiabatic_index.txt'),
                         numpy.loadtxt('amplitude.txt'),
                         numpy.loadtxt('width.txt'))
    analytic = dict()
    analytic['density'] = [spat_prof.density.eval(x-spat_prof.c0*time) \
                            for x in final_data['geometry']['x_coordinate']]
    analytic['pressure'] = [spat_prof.pressure.eval(x-spat_prof.c0*time) \
                            for x in final_data['geometry']['x_coordinate']]
    analytic['x_velocity'] = [spat_prof.velocity.eval(x-spat_prof.c0*time) \
                            for x in final_data['geometry']['x_coordinate']]

    if graphic_flag:
        for n, field in enumerate(['density','pressure','x_velocity']):
            pylab.subplot(3,1,n+1)
            pylab.plot(final_data['x_coordinate'],
                       final_data[field],'.')
            pylab.plot(final_data['x_coordinate'],
                       analytic[field],'.')
            pylab.plot(initial['x_coordinate'],
                       initial[field], '.')
        pylab.show()

    error_norm = dict()
    for i in analytic.keys():
        error_norm[i] = calc_error_norm(
            final_data['hydrodynamic'][i],analytic[i])

    f = open('gradesheet.txt','w')
    for i in analytic.keys():
        f.write(str(error_norm[i])+'\n')
    f.close()

    return error_norm['density']<5e-8 and \
        error_norm['pressure']<5e-8 and \
        error_norm['x_velocity']<5e-8

if __name__=='__main__':
    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
