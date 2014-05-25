def l1_norm(a1, a2):

    return sum([abs(x-y) for x,y in zip(a1, a2)])/len(a1)

def main(graphic_flag=False):

    import h5py
    import numpy
    import imp
    import os
    afd = imp.load_source('afd',os.environ['RICH_ROOT']+'/analytic/afd.py')
    import argparse

    initial = h5py.File('initial.h5')
    final_numeric = h5py.File('final.h5')

    init_parsed = {}
    init_parsed['grid'] = initial['grid']
    init_parsed['adiabatic index'] = numpy.loadtxt('adiabatic_index.txt')
    init_parsed['ambient'] = {}
    init_parsed['ambient']['density'] = numpy.loadtxt('ambient_density.txt')
    init_parsed['ambient']['pressure'] = numpy.loadtxt('ambient_pressure.txt')
    init_parsed['ambient']['velocity'] = numpy.loadtxt('ambient_velocity.txt')
    init_parsed['pert'] = {}
    init_parsed['pert']['density'] = [d-init_parsed['ambient']['density']
                                      for d in initial['density']]
    init_parsed['pert']['pressure'] = [d-init_parsed['ambient']['pressure']
                                       for d in initial['pressure']]
    init_parsed['pert']['velocity'] = [d-init_parsed['ambient']['velocity']
                                       for d in initial['x_velocity']]
    time = numpy.array(final_numeric['time'])[0]

    final_exact = afd.exact_time_advance(init_parsed, time)
    final_2o = afd.second_order_time_advance(init_parsed, time)

    if graphic_flag:

        import pylab
        import matplotlib

        for i, field in enumerate(['density',
                                   'pressure',
                                   'velocity']):
            pylab.subplot(3,1,i+1)
            pylab.gca().ticklabel_format(axis='y',style='sci',
                                         scilimits=(-2,2))
            pylab.tick_params(axis='x',
                              which='both',
                              bottom='off',
                              labelbottom='off')
            pylab.plot(initial['grid'],
                       init_parsed['pert'][field])
            pylab.plot(final_numeric['grid'],
                       [var - init_parsed['ambient'][field] for var in 
                        final_numeric[field.replace('velocity','x_velocity')]],
                       'k')
            pylab.plot(initial['grid'],
                       final_2o[field], linewidth=3, alpha=0.5)
            pylab.plot(initial['grid'], final_exact[field], 'r')
            pylab.ylabel(field)
        pylab.tick_params(axis='x',
                          which='both',
                          bottom='on',
                          labelbottom='on')
        pylab.xlabel('x')
        pylab.show()

    l1_list = {}
    for field in final_2o:
        l1_list[field] = l1_norm([init_parsed['ambient'][field]+var
                                  for var in final_2o[field]],
                                 final_numeric[field.replace('velocity','x_velocity')])

    f = open('gradehseet.txt','w')
    for field in ['density', 'pressure', 'velocity']:
        f.write(str(l1_list[field])+'\n')
    f.close()
        
    return abs(l1_list['density'])<2e-6 and \
        abs(l1_list['pressure'])<2e-6 and \
        abs(l1_list['velocity'])<2e-6

if __name__=='__main__':

    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')

