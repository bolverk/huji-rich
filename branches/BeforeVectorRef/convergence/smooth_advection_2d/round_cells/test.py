import math
import logging

logging.basicConfig(filename="test_log.txt",level=logging.INFO)

def l1_test(a1,a2):

    return sum([abs(x-y) for x,y in zip(a1,a2)])/len(a1)

drift_velocity = 1

def density_prof(x,t):

    v = drift_velocity
    return 2+math.sin(2*math.pi*(x-v*t))

def pressure_prof(x):

    return 1

def velocity_prof(x):

    return 1

def main():

    import numpy

    logging.debug('Start of the test')

    rawd = numpy.loadtxt("sim_results.txt")

    logging.debug('Finished reading data from simulation')

    numeric = {'time':numpy.loadtxt('time.txt'),
               'coordinate':rawd[:,0],
               'density':rawd[:,1],
               'pressure':rawd[:,2],
               'velocity':rawd[:,3]}

    logging.debug('Finished putting simulation data in dictionary')

    analytic = {'density':[density_prof(x,numeric['time']) 
                           for x in numeric['coordinate']],
                'pressure':[pressure_prof(x) for x in numeric['coordinate']],
                'velocity':[velocity_prof(x) for x in numeric['coordinate']]}

    logging.debug('Finished calculating analytic profiles')

    if False:
        import pylab

        for i,field in enumerate(analytic):
            pylab.subplot(3,1,i+1)
            pylab.plot(numeric['coordinate'],
                       numeric[field])
            pylab.plot(numeric['coordinate'],
                       analytic[field])
            pylab.ylabel(field)

        pylab.xlabel('Distance')
        pylab.show()

    gof = {}
    for field in analytic:
        gof[field] = l1_test(numeric[field],analytic[field])

    f = open('gradesheet.txt','w')
    for field in ['density','pressure','velocity']:
        f.write(str(gof[field])+'\n')
    f.close()

    return gof['density']<0.1 and gof['pressure']<0.1 and gof['velocity']<0.1

if __name__=='__main__':

    print main()
