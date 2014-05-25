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

    import numpy

    init_cond = \
        read_x_prof('x_prof_initial.txt')
    time = numpy.loadtxt('time.txt')
    final_cond = \
        read_x_prof('x_prof_final.txt')
    spat_prof = SpatProf(1.,3./5.,5./3.,1e-6,1.)
    analytic = dict()
    analytic['density'] = [spat_prof.density.eval(x-spat_prof.c0*time) \
                            for x in final_cond['coordinate']]
    analytic['pressure'] = [spat_prof.pressure.eval(x-spat_prof.c0*time) \
                            for x in final_cond['coordinate']]
    analytic['velocity'] = [spat_prof.velocity.eval(x-spat_prof.c0*time) \
                            for x in final_cond['coordinate']]

    if False:

        import pylab

        for i,v in enumerate(analytic):

            pylab.subplot(3,1,i+1)
            pylab.plot(final_cond['coordinate'],
                       final_cond[v],'.')
            pylab.plot(final_cond['coordinate'],
                       analytic[v],'.')
            pylab.ylabel(v)

        pylab.ylabel('x')
        pylab.show()

    error_norm = dict()
    for i in analytic.keys():
        error_norm[i] = calc_error_norm(final_cond[i],analytic[i])

    f = open('gradesheet.txt','w')
    for field in ['density','pressure','velocity']:
        f.write(str(error_norm[i])+'\n')
    f.close()

    return error_norm['density']<5e-8 and \
        error_norm['pressure']<5e-8 and \
        error_norm['velocity']<5e-8

if __name__=='__main__':
    print main()
