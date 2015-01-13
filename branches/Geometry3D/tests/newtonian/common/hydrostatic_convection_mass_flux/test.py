#! /usr/bin/python

def main():

    import numpy
    
    rawd = numpy.loadtxt('res.txt')
    return rawd[0]==rawd[1]

if __name__=='__main__':
    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')


