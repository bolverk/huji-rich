#! /usr/bin/python

def goodness_of_fit(a1, a2):

    import numpy

    return sum(numpy.abs(numpy.array(a1)-numpy.array(a2)))/len(a1)

def consolidate_single(fname):

    import h5py
    import numpy

    with h5py.File(fname,'r+') as f:
        return {field:numpy.array(f[field]) for field in f}

def consolidate_multiple(f_list):

    assert(len(f_list)>1)

    import numpy

    parts_list = [consolidate_single(fname) for fname in f_list]
    return {key:numpy.concatenate([part[key] for part in parts_list]) for key in parts_list[0]}

def find_max_pos_par(x_list, y_list, window=2):

    import numpy

    mi = numpy.argmax(y_list)

    fit_data = numpy.polyfit(x_list[mi-window:mi+window+1],
                             y_list[mi-window:mi+window+1],
                             2)
    return -0.5*fit_data[1]/fit_data[0]

def calc_shock_front(numeric):

    import numpy

    radially_sorted = {field:numpy.array([x for (y,x) in sorted(zip(numeric['radius'],numeric[field]))])
                       for field in ['radius','density','pressure','velocity']}
    res = {'radius':find_max_pos_par(radially_sorted['radius'],
                                     radially_sorted['pressure'])}
    for field in ['density','pressure','velocity']:
        res[field] = numpy.interp(res['radius'],
                                  radially_sorted['radius'],
                                  radially_sorted[field])
    return res

class SedovTaylorProfiles:

    def __init__(self, upstream, shock_front, g, w, n ,nip=1000):

        import numpy
        import imp
        import os

        sedov_taylor = imp.load_source('sedov_taylor',
                                       os.environ['RICH_ROOT']+'/analytic/sedov_taylor.py')

        self.upstream = upstream
        ssv_list = numpy.linspace(1e-6+1./g,
                                  2./(g+1.),
                                  num=nip)
        self.shock_radius = shock_front['radius']
        self.tables = {'radius':numpy.array([shock_front['radius']*sedov_taylor.vtoz(v,w,g,n) for v in ssv_list]),
                       'density':numpy.array([shock_front['density']*sedov_taylor.vtod(v,w,g,n)/((g+1.)/(g-1.)) for v in ssv_list]),
                       'pressure':numpy.array([shock_front['pressure']*sedov_taylor.vtop(v,w,g,n)/(2./(g+1.)) for v in ssv_list]),
                       'velocity':numpy.array([shock_front['velocity']*sedov_taylor.vtoz(v,w,g,n)*v/(2./(g+1.)) for v in ssv_list])}

    def calc(self, field, r):

        import numpy

        if r>self.shock_radius:
            return self.upstream[field]
        else:
            return numpy.interp(r,
                                self.tables['radius'],
                                self.tables[field])

def main():

    import numpy
    import glob
    import pickle

    np = len(glob.glob('process_*_final.h5'))

    if np>0:
        numeric = consolidate_multiple(glob.glob('process_*_final.h5'))
    else:
        numeric = consolidate_single('final.h5')
    numeric['radius'] = numeric['grid']
    numeric['velocity'] = numeric['x_velocity']

    st_prof = SedovTaylorProfiles({'density':1,
                                   'pressure':1e-9,
                                   'velocity':0},
                                  calc_shock_front(numeric),
                                  5./3.,0,2)
    analytic = {field:numpy.array([st_prof.calc(field,r) for r in numeric['radius']]) for field in ['density','pressure','velocity']}

    l1_data = {field:goodness_of_fit(analytic[field],
                                     numeric[field])
               for field in analytic}

    f = open('gradesheet.txt','w')
    for field in ['density','pressure','velocity']:
        f.write(str(l1_data[field])+'\n')
    f.close()

    pickle.dump(numeric, open('numeric.pkl','wb'))
    pickle.dump(analytic, open('analytic.pkl','wb'))

    if False:
        import pylab
        for n,f in enumerate(['density','pressure','velocity']):
            pylab.subplot(3,1,n+1)
            pylab.plot(numeric['radius'],
                       numeric[f],'.')
            pylab.plot(numeric['radius'],
                       analytic[f],'.')
        pylab.show()

    return l1_data['density']<0.5 and l1_data['pressure']<14 and l1_data['velocity']<1.5

import sys
if __name__=='__main__':

    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
