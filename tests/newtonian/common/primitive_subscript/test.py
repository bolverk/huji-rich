#! /usr/bin/python

def main():

    import numpy
    
    rawd = numpy.loadtxt('res.txt')
    return rawd[0]==1 and \
        rawd[1]==2 and \
        rawd[2]==3 and \
        rawd[3]==4 and \
        rawd[4]==5 and \
        rawd[5]==6

if __name__=='__main__':
    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')

