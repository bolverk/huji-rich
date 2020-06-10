#! /usr/bin/python

def main():

    import os
    import numpy
    import h5py

    with h5py.File('initial.h5', 'r') as f:
        c_initial = numpy.array(f['colour'])
    with h5py.File('final.h5', 'r') as f:
        c_final = numpy.array(f['colour'])

    diff = numpy.sqrt(numpy.sum((c_final - c_initial)**2))/len(c_final)

    return diff<0.008

if __name__=='__main__':

    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')

