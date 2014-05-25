#! /usr/bin/python

def main():

    import numpy
    pa = 0.224614
    va = 0
    rawd = numpy.loadtxt('res.txt')
    pn = rawd[0]
    vn = rawd[1]

    f = open('report.txt','w')
    f.write('Analytic pressure = '+str(pa)+'\n')
    f.write('Numeric pressure = '+str(pn)+'\n')
    f.write('Analytic velocity = '+str(va)+'\n')
    f.write('Numeric velocity = '+str(vn)+'\n')
    f.close()

    return abs(pa-pn)/(pa+pn)<0.02 and abs(vn)<0.01

if __name__=='__main__':

    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')

