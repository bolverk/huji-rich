#! /usr/bin/python

def goodness_of_fit(a1, a2):

    import math

    diff2 = [(x-y)**2 for x,y in zip(a1,a2)]
    return math.sqrt(sum(diff2)/(max(a1)-min(a1))**2/len(a1))

def main():

    import os
    import imp
    enrs = imp.load_source('enrs',os.environ['RICH_ROOT']+'/analytic/enrs.py')

    import numpy
    import h5py

    left = enrs.Primitive(1,10,0);
    right = enrs.Primitive(1,1,0);
    prof = enrs.RiemannProfile(left,right,5./3.)
    h5f = h5py.File('final.h5')
    t = h5f['time']
    offset = 0.5

    x_list = h5f['grid']
    numeric = {}
    numeric['density'] = h5f['density']
    numeric['pressure'] = h5f['pressure']
    numeric['velocity'] = h5f['x_velocity']
    analytic = {}
    analytic['density'] = \
        [prof.CalcPrim((i-offset)/t).Density for i in x_list]
    analytic['pressure'] = \
        [prof.CalcPrim((i-offset)/t).Pressure for i in x_list]
    analytic['velocity'] = \
        [prof.CalcPrim((i-offset)/t).Velocity for i in x_list]

    gof = {}
    for varname in analytic:
        gof[varname] = \
            goodness_of_fit(numeric[varname],
                            analytic[varname])

    f = open('gradesheet.txt','w')
    for varname in ['density', 'pressure', 'velocity']:
        f.write(str(gof[varname])+'\n')
    f.close()

    return gof['density'] < 0.11 and \
        gof['pressure'] < 0.05 and \
        gof['velocity'] < 0.11

if __name__=='__main__':

    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')

