def extract_contour(x_list, y_list, z_list, z_value):

    import matplotlib
    matplotlib.use('Agg')
    import pylab
    import numpy

    temp = pylab.tricontour(x_list, y_list, z_list, levels=[z_value])
    pylab.clf()
    return numpy.concatenate([itm.vertices for itm
                              in temp.collections[0].get_paths()])

def main():

    import numpy
    import math
    import imp
    import os
    oblique_shock = imp.load_source(\
        'oblique_shock',\
            os.environ['RICH_ROOT']+'/analytic/oblique_shock.py')
    import h5py

    with h5py.File('final.h5','r+') as f:
        raw = {}
        for field in ['x_coordinate','y_coordinate','y_velocity','wedge','density','pressure']:
            raw[field] = numpy.array(f[field])
        numeric = {}
        for field in raw:
            numeric[field] = raw[field][raw['wedge']<0.5]
    wedge_angle = numpy.loadtxt('wedge_angle.txt')
    adiabatic_index =numpy.loadtxt('adiabatic_index.txt')
    mach = numpy.loadtxt('mach_number.txt')
    v_th = 0.3
    numeric['entropy'] = numeric['pressure']/numeric['density']**adiabatic_index
    shock_front = extract_contour(numeric['x_coordinate'],
                                  numeric['y_coordinate'],
                                  numeric['entropy'],
                                  1.01)
    x_far = shock_front.T[0][shock_front.T[0]>0]
    y_far = shock_front.T[1][shock_front.T[0]>0]
    lfit = numpy.polyfit(x_far, y_far,1)
    shock_angle = numpy.arctan(lfit[0])
    exact_shock_angle = oblique_shock.calc_shock_angle\
        (mach, adiabatic_index, wedge_angle)

    if True:
        print('mach number = '+str(mach))
        print('adiabatic index = '+str(adiabatic_index))
        print('wedge angle = '+str(wedge_angle*180/math.pi))
        print('shock angle = '+str(shock_angle*180/math.pi))
        print('calculated shock angle = '+\
              str(exact_shock_angle*180/math.pi))

    difrat = abs(shock_angle-exact_shock_angle)\
        /exact_shock_angle
    
    f = open('gradesheet.txt','w')
    f.write(str(difrat)+'\n')
    f.close()

    return difrat<0.35

if __name__=='__main__':
    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
