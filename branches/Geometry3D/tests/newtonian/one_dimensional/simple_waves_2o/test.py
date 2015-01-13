def L1_error_norm(a1, a2):

    abs_diff = [abs(x-y) for x,y in zip(a1,a2)]
    return sum(abs_diff)/len(a1)

def chi2_error_norm(a1,a2):

    import math

    diff2 = [(x-y)**2 for x,y in zip(a1,a2)]
    return math.sqrt(sum(diff2)/len(a1))

def error_norm(a1, a2, test_name='L1'):

    if 'L1'==test_name:
        return L1_error_norm(a1,a2)
    elif 'chi2'==test_name:
        return chi2_error_norm(a1,a2)
    else:
        raise NameError('Unknown test name '+test_name)

def calc_sound_speed(d,p,g):

    import math

    return math.sqrt(g*p/d)

def simple_waves_propagate(x_list,
                           d_list,
                           p_list,
                           v_list,
                           time, g,
                           direction):

    if 0==direction or 1==direction or -1==direction:
        c_list = [calc_sound_speed(d,p,g) for d,p in zip(d_list,p_list)]
        propag_vel = [v+direction*c for v,c in zip(v_list,c_list)]
        res = [x+time*u for x,u in zip(x_list,propag_vel)]
        return res
    else:
        raise NameError('direction must be either 0,1 or -1')

def main():

    import numpy
    import h5py

    initial_data = h5py.File('initial.h5')
    x_init = initial_data['grid']
    d_init = initial_data['density']
    p_init = initial_data['pressure']
    v_init = initial_data['x_velocity']

    final_data = h5py.File('final.h5')
    time = numpy.array(final_data['time'])[0]
    x_numeric = final_data['grid']
    d_numeric = final_data['density']
    p_numeric = final_data['pressure']
    v_numeric = final_data['x_velocity']

    x_exact = simple_waves_propagate(x_init,
                                     d_init,
                                     p_init,
                                     v_init,
                                     time, 5./3.,1)

    if False:
        
        import pylab

        pylab.subplot(311)
        pylab.plot(x_numeric, d_numeric)
        pylab.plot(x_exact, d_init)
        pylab.plot(x_init,d_init)

        pylab.subplot(312)
        pylab.plot(x_numeric,p_numeric)
        pylab.plot(x_exact,p_init)
        pylab.plot(x_init,p_init)

        pylab.subplot(313)
        pylab.plot(x_numeric,v_numeric)
        pylab.plot(x_exact,v_init)
        pylab.plot(x_init,v_init)

        pylab.show()

    test_name = 'L1'
    xc = [x for x in x_exact if x>min(x_numeric) and x<max(x_numeric)]
    d1 = numpy.interp(xc,x_numeric,d_numeric)
    d2 = numpy.interp(xc,x_exact,d_init)
    errnorm1 = error_norm(d1,d2,test_name)

    p1 = numpy.interp(xc,x_numeric,p_numeric)
    p2 = numpy.interp(xc,x_exact,p_init)
    errnorm2 = error_norm(p1,p2,test_name)

    v1 = numpy.interp(xc,x_numeric,v_numeric)
    v2 = numpy.interp(xc,x_exact,v_init)
    errnorm3 = error_norm(v1,v2,test_name)

    f = open('gradesheet.txt','w')
    f.write(str(errnorm1)+'\n')
    f.write(str(errnorm2)+'\n')
    f.write(str(errnorm3)+'\n')
    f.close()

    return errnorm1<0.3 and errnorm2 < 1.7 and errnorm3 < 0.22
                                     
if __name__=='__main__':
    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')

