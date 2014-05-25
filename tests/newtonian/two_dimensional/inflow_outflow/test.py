def main():

    import numpy

    rawd = numpy.loadtxt('result.txt')
    return abs(rawd)<0.0002

if __name__=='__main__':

    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')


