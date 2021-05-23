#! /usr/bin/python

import logging

logging.basicConfig(level=logging.DEBUG)

def goodness_of_fit(a1, a2):

    import math

    diff2 = [(x-y)**2 for x,y in zip(a1,a2)]
    return math.sqrt(sum(diff2)/(max(a1)-min(a1))/len(a1))

def load_file(fname):

    import h5py
    import numpy

    res = {}
    with h5py.File(fname, 'r') as f :
        logging.debug([field for field in f])
        return {field:numpy.array(f[field])
                for field in f}

def main():

    import numpy
    import os

    initial = load_file('initial.h5')
    final = load_file('final.h5')

    for state in [initial, final]:
        state['mass'] = state['Volume']*state['Density']
        state['energy'] = (state['InternalEnergy']+
                           0.5*state['Vx']**2+
                           0.5*state['Vy']**2+
                           0.5*state['Vz']**2)*state['mass']
        state['x_momentum'] = state['mass']*state['Vx']
        state['y_momentum'] = state['mass']*state['Vy']
        state['z_momentum'] = state['mass']*state['Vz']

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
