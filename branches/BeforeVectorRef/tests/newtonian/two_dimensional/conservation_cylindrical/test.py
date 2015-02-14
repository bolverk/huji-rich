#! /usr/bin/python

def all_equal(ar):
    """
    Checks that all terms in the array are equal
    Input:
    ar - Numerical array
    """

    for i in ar:
        if i!=ar[0]:
            return False
    return True

def main():

    import numpy

    mass, xmom, ymom, enr, trc = \
        numpy.loadtxt('res.txt',unpack=True);

    f = open('gradesheet.txt','w')
    f.write('mass '+str(all_equal(mass))+'\n')
    f.write('xmom '+str(all_equal(xmom))+'\n')
    f.write('ymom '+str(all_equal(ymom))+'\n')
    f.write('enr '+str(all_equal(enr))+'\n')
    f.write('trc '+str(all_equal(trc))+'\n')
    f.close()

    return all_equal(mass) and \
        all_equal(ymom) and \
        all_equal(enr) and \
        all_equal(trc)

if __name__=='__main__':
    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
