#! /usr/bin/python

def main():

    import numpy
    pa = 0.224614
    va = 0
    rawd = numpy.loadtxt('res.txt')
    pn = rawd[0]
    vn = rawd[1]
    return abs(pa-pn)/(pa+pn)<0.01 and abs(va-vn)/(va+vn)<0.01

import sys
if __name__=='__main__':
    print main()
