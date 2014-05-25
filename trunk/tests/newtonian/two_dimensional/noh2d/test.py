def l1_comparison(a1, a2):

    return sum([abs(x-y) for x,y in zip(a1,a2)])/len(a1)

def shock_velocity(v0,g):

    return (g-1.0)/2.0

def calc_analytic_d(d0,v0,g,r,t):

    vs = shock_velocity(v0,g)

    if r>vs*t:
        return d0*(1.0+v0*t/r)
    else:
        return d0*((g+1.0)/(g-1.0))**2

def calc_analytic_v(d0,v0,g,r,t):

    vs = shock_velocity(v0,g)

    if r>vs*t:
        return v0
    else:
        return 0

def calc_analytic_p(d0,v0,g,r,t):

    vs = shock_velocity(v0,g)

    if r>vs*t:
        return 0
    else:
        return 0.5*d0*v0**2*(g+1.0)**2/(g-1.0)

def main():

    import numpy
    import math
    import h5py

    f = h5py.File('final.h5')

    x_list = f['x_coordinate']
    y_list = f['y_coordinate']
    r_list = [math.sqrt(x**2+y**2) for x,y in zip(x_list,y_list)]
    d_list = f['density']
    p_list = f['pressure']
    vx_list = f['x_velocity']
    vy_list = f['y_velocity']
    vm_list = [math.sqrt(vx**2+vy**2) for vx,vy in zip(vx_list, vy_list)]
    t = numpy.array(f['time'])[0]
    g = numpy.loadtxt('adiabatic_index.txt')
    vc = numpy.loadtxt('collapse_velocity.txt')
    di = numpy.loadtxt('initial_density.txt')

    analytic_d = [calc_analytic_d(di,vc,g,r,t) for r in r_list]
    analytic_p = [calc_analytic_p(di,vc,g,r,t) for r in r_list]
    analytic_v = [calc_analytic_v(di,vc,g,r,t) for r in r_list]

    l1_density = l1_comparison(analytic_d, d_list)
    l1_pressure = l1_comparison(analytic_p, p_list)
    l1_velocity = l1_comparison(analytic_v, vm_list)

    f = open('gradesheet.txt','w')
    f.write(str(l1_density)+'\n')
    f.write(str(l1_pressure)+'\n')
    f.write(str(l1_velocity)+'\n')
    f.close()

    return l1_density<0.8 and \
        l1_pressure<0.23 and \
        l1_velocity<0.05

if __name__ == '__main__':

    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')

