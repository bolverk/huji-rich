#! /usr/bin/python

def main():

    import numpy
    pa = 5.11232
    va = 1.53795
    rawd = numpy.loadtxt('res.txt')
    pn = rawd[0]
    vn = rawd[1]
    return abs(pa-pn)/(pa+pn)<0.01 and abs(va-vn)/(va+vn)<0.01

if __name__=='__main__':
    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')

