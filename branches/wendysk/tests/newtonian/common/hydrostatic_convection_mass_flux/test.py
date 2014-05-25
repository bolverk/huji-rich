#! /usr/bin/python

def main():

    import numpy
    
    rawd = numpy.loadtxt('res.txt')
    return rawd[0]==rawd[1]

import sys
if __name__=='__main__':
    print main()
