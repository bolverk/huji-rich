def main():

    import numpy
    import math
    import h5py
    import glob
	
    with h5py.File('final.h5','r+') as f:
        density_list = numpy.array(f['hydrodynamic']['density'])
        chi_2 = numpy.sqrt(sum((density_list-1)**2))/len(density_list)
    mass = numpy.loadtxt('mass.txt')
    return chi_2<4e-6 and abs(mass-4)<5e-8

 
if __name__ == '__main__':

    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')

