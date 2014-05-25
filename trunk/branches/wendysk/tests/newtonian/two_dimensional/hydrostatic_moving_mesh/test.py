#! /usr/bin/python

import numpy

def is_increasing(ar):
    """
    Checks that the array is in increasing order
    Input:
    ar - Numerical array
    """

    temp = numpy.diff(ar)
    for i in temp:
        if i<0:
            return False
    return True

def is_decreasing(ar):
    """
    Checks that the array is in decreasing order
    Input:
    ar - Numerical array
    """

    temp = numpy.diff(ar)
    for i in temp:
        if i>0:
            return False
    return True

def main():

    rawd = numpy.loadtxt('res.txt');
    return is_increasing(rawd[:,0]) and \
        is_decreasing(rawd[:,1])

import sys
if __name__=='__main__':
    print main()
