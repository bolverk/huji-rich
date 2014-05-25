def l1_error_margin(a1,a2):

    abs_dif = [abs(x-y) for x,y in zip(a1,a2)]
    return sum(abs_dif)/len(abs_dif)

def main():

    graphic_flag = False

    if graphic_flag:
        import matplotlib
        matplotlib.use('Qt4Agg')
        import pylab
    
    import numpy
    import math
    import imp
    import h5py
    sedov_taylor = imp.load_source(
        'sedov_taylor',\
            '../../analytic/sedov_taylor.py')

    h5f = h5py.File('final.h5')
    x_list = h5f['x_coordinate']
    y_list = h5f['y_coordinate']
    p_list = h5f['pressure']

    x_c = 0
    y_c = 0.5
    r_list = [math.sqrt((y-y_c)**2+(x-x_c)**2)
              for x,y in zip(x_list,y_list)]
    p_front = numpy.max(p_list)
    r_front = r_list[numpy.argmax(p_list)]
    p_back = p_list[numpy.argmin(r_list)]

    r_relevant = [r for r in r_list if r<r_front]
    p_relevant = [p_list[i] for i in range(len(r_list)) 
                  if r_list[i]<r_front]

    # Analytic results
    g = 5./3.
    w = 0
    n = 3
    v_for_interp = numpy.linspace(1e-9+1./g,2./(g+1),num=1000)
    r_for_interp = [r_front*sedov_taylor.vtoz(v,w,g,n) 
                    for v in v_for_interp]
    p_for_interp = [p_back*sedov_taylor.vtop(v,w,g,n)/
                    sedov_taylor.vtop(1./g+1e-6,w,g,n)
                    for v in v_for_interp]

    l1 = l1_error_margin(p_relevant,
                         numpy.interp(r_relevant,
                                      r_for_interp,
                                      p_for_interp))

    f = open('gradesheet.txt','w')
    f.write(str(l1)+'\n')
    f.close()

    if graphic_flag:
        pylab.plot(r_relevant, p_relevant,'.')
        pylab.plot(r_for_interp,p_for_interp)
        pylab.show()

    return l1<20

if __name__=='__main__':
    print main()
