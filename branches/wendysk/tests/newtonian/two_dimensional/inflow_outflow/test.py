def main():

    import numpy

    rawd = numpy.loadtxt('result.txt')
    return abs(rawd)<0.0002

if __name__=='__main__':
    print main()
