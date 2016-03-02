def extract_contour(x_list, y_list, z_list, z_value):

    import matplotlib
    matplotlib.use('Agg')
    import pylab
    import numpy

    temp = pylab.tricontour(x_list, y_list, z_list, levels=[z_value])
    pylab.clf()
    return numpy.concatenate([itm.vertices for itm
                              in temp.collections[0].get_paths()])

def consolidate_single(fname):

    import h5py
    import numpy

    res = {}
    with h5py.File(fname,'r+') as f:
        structure = {'geometry':['x_coordinate','y_coordinate'],
                     'hydrodynamic':['density','pressure'],
                     'stickers':['wedge']}
        for field in structure:
            for subfield in structure[field]:
                res[subfield] = numpy.array(f[field][subfield])
    return res

def consolidate_multiple(fname_list):

    import numpy
    
    part_list = [consolidate_single(fname) for fname in fname_list]
    res = {}
    for field in part_list[0]:
        res[field] = numpy.concatenate([p[field] for p in part_list])
    return res

def main():

    import numpy
    import math
    import imp
    import os
    oblique_shock = imp.load_source(\
        'oblique_shock',\
            os.environ['RICH_ROOT']+'/analytic/oblique_shock.py')
    import h5py
    import glob

    if len(glob.glob('process_*_final.h5'))>0:
        raw = consolidate_multiple(glob.glob('process_*_final.h5'))
    else:
        raw = consolidate_single('final.h5')
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
                                  1.1)
    x_far = shock_front.T[0][shock_front.T[0]>0]
    y_far = shock_front.T[1][shock_front.T[0]>0]
    lfit = numpy.polyfit(x_far, y_far,1)
    shock_angle = numpy.arctan(lfit[0])
    exact_shock_angle = oblique_shock.calc_shock_angle\
        (mach, adiabatic_index, wedge_angle)

    difrat = abs(shock_angle-exact_shock_angle)\
        /exact_shock_angle
    
    f = open('gradesheet.txt','w')
    f.write(str(difrat)+'\n')
    f.close()

    return difrat<0.4

if __name__=='__main__':
    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
