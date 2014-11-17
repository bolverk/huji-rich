def dist_sqr(arr):

    import numpy

    assert(len(arr)==2)

    return sum(numpy.array(arr)**2)

def main():

    import numpy

    threshold = 1e-6

    raw = numpy.loadtxt('output.txt')

    ref = [[2,3],
           [1,0],
           [-1,0],
           [1,0],
           [2,0],
           [0,1],
           [0,-1],
           [0,1],
           [0,2]]

    for n,a in zip(raw,ref):
        if dist_sqr(n-a)>threshold:
            return False

    return True

if __name__ == '__main__':

    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
    
