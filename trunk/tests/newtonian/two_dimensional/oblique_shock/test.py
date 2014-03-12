def main():

    import numpy
    import math
    import imp
    oblique_shock = imp.load_source(\
        'oblique_shock',\
            '../../analytic/oblique_shock.py')
    import h5py

    h5f = h5py.File('final.h5')
    x_list = h5f['x_coordinate']
    y_list = h5f['y_coordinate']
    v_list = h5f['y_velocity']
    wedge_angle = numpy.loadtxt('wedge_angle.txt')
    adiabatic_index =numpy.loadtxt('adiabatic_index.txt')
    mach = numpy.loadtxt('mach_number.txt')
    x_c = -0.5
    y_c = -0.5
    v_th = 0.3
    x_leftmost = x_list[0]
    y_downmost = y_list[0]
    shock_angle = 0
    for x,y,v in zip(x_list,y_list,v_list):
        if v>v_th:
            slope = (y-y_c)/(x-x_c)
            shock_angle = max([math.atan(slope),
                               shock_angle])

    exact_shock_angle = oblique_shock.calc_shock_angle\
        (mach, adiabatic_index, wedge_angle)

    if False:
        print 'mach number = '+str(mach)
        print 'adiabatic index = '+str(adiabatic_index)
        print 'wedge angle = '+str(wedge_angle*180/math.pi)
        print 'shock angle = '+str(shock_angle*180/math.pi)
        print 'calculated shock angle = '+\
            exact_shock_angle

    difrat = abs(shock_angle-exact_shock_angle)\
        /exact_shock_angle
    
    f = open('gradesheet.txt','w')
    f.write(str(difrat)+'\n')
    f.close()

    return difrat<0.03

if __name__=='__main__':
    print main()
