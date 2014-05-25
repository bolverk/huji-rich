def main():

    import numpy

    results = numpy.loadtxt('result.txt')
    mass = results[0]
    momentum_x = results[1]
    momentum_y = results[2]
    min_density = results[3]
    max_density = results[4]

    return (abs(1-mass)<0.0001 and
            abs(1-momentum_x)<0.0001 and
            abs(momentum_y)<0.0001 and
            min_density>0.95 and
            max_density<1.02)

if __name__ == '__main__':

    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
            
