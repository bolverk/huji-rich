#! /usr/bin/python

def chi_squared(a1, a2):

    import math

    diff2 = [(x-y)**2 for x,y in zip(a1,a2)]
    return math.sqrt(sum(diff2)/len(a1))

def L1_error_norm(a1,a2):

    import math
    dfrnc = [abs(x-y) for x,y in zip(a1,a2)]
    return sum(dfrnc)/len(a1)

def error_norm(a1,a2,test_name='L1'):

    if('L1'==test_name):
        return L1_error_norm(a1,a2)
    elif('chi2'==test_name):
        return chi_squared(a1,a2)
    else:
        raise NameError('unknown test type '+test_name)

def main():

    import os
    import imp
    enrs = imp.load_source('enrs','../../analytic/enrs.py')

    import numpy

    left = enrs.Primitive(1,1,-1);
    right = enrs.Primitive(1,1,1);
    prof = enrs.RiemannProfile(left,right,5./3.)
    t = numpy.loadtxt('time.txt')
    offset = 0.5

    x = numpy.loadtxt('cell_centres.txt')
    d = numpy.loadtxt('densities.txt')
    p = numpy.loadtxt('pressures.txt')
    v = numpy.loadtxt('velocities.txt')
    da = [prof.CalcPrim((i-offset)/t).Density for i in x]
    pa = [prof.CalcPrim((i-offset)/t).Pressure for i in x]
    va = [prof.CalcPrim((i-offset)/t).Velocity for i in x]

    test_name = 'L1'
    gof1 = error_norm(d,da,test_name)
    gof2 = error_norm(p,pa,test_name)
    gof3 = error_norm(v,va,test_name)

    f = open('gradesheet.txt','w')
    f.write(str(gof1)+'\n')
    f.write(str(gof2)+'\n')
    f.write(str(gof3)+'\n')
    f.close()

    return gof1<0.06 and gof2<0.02 and gof3 < 0.04

import sys
if __name__=='__main__':
    print main()
