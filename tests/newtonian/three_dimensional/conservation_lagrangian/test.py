#! /usr/bin/python

def goodness_of_fit(a1, a2):

    import math

    diff2 = [(x-y)**2 for x,y in zip(a1,a2)]
    return math.sqrt(sum(diff2)/(max(a1)-min(a1))/len(a1))

def load_file(fname):

    import h5py
    import numpy

    res = {}
    with h5py.File(fname, 'r') as f :
        for field in f['hydrodynamic']:
            res[field] = numpy.array(f['hydrodynamic'][field])
    return res

def main():

    import numpy
    import os

    initial = load_file('initial.h5')
    final = load_file('final.h5')

    for field in ['mass', 'energy', 'x_momentum', 'y_momentum', 'z_momentum']:
        if abs(sum(initial[field]) - sum(final[field]))>1e-5:
            return False

    return True

if __name__=='__main__':

    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
