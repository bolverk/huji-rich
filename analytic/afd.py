import numpy
import cmath

def sgn(x):

    if x==0:
        return 0
    else:
        return x/abs(x)

def single_mode(x,t,k,omega,v,g0,dx,xi):

    return cmath.exp(1j*k*x)*g0*dx*xi*\
        (cmath.exp(1j*omega*t)-\
             (1-(1-cmath.exp(-1j*k*dx))*xi)**(t*v/xi/dx))/\
             (v*(-1+cmath.exp(1j*dx*xi*omega/v)+\
                      xi-xi*cmath.exp(-1j*k*dx)))

def primitives_to_conserved(hydro_data):
    """
    Converts primitive variables to variables that are conserved along streamlines 
    """

    g = hydro_data['adiabatic index']
    sound_speed = g*hydro_data['ambient']['pressure']/hydro_data['ambient']['density']
    res = {}
    res['positive riemann invariant'] = [dv + dp/(sound_speed*hydro_data['ambient']['density'])
                                         for dp, dv in zip(hydro_data['pert']['pressure'],
                                                           hydro_data['pert']['velocity'])]
    res['negative riemann invariant'] = [dv - dp/(sound_speed*hydro_data['ambient']['density'])
                                         for dp, dv in zip(hydro_data['pert']['pressure'],
                                                           hydro_data['pert']['velocity'])]
    res['entropy'] = [dp/hydro_data['ambient']['pressure'] - g*dd/hydro_data['ambient']['density']
                      for dd, dp in zip(hydro_data['pert']['density'],
                                        hydro_data['pert']['pressure'])]

    return res

def conserved_to_primitives(initial, conserved):

    import math

    g = initial['adiabatic index']
    res = {}
    res['velocity'] = [0.5*(jp+jm) for jp, jm in zip(conserved['positive riemann invariant'],
                                                     conserved['negative riemann invariant'])]
    sound_speed = math.sqrt(g*initial['ambient']['pressure']/initial['ambient']['density'])
    res['pressure'] = [0.5*sound_speed*initial['ambient']['density']*(jp-jm) 
                       for jp, jm in zip(conserved['positive riemann invariant'],
                                         conserved['negative riemann invariant'])]
    res['density'] = [(initial['ambient']['density']/g)*(dp/initial['ambient']['pressure']-ds)
                      for dp,ds in zip(res['pressure'],conserved['entropy'])]
    return res

def to_k_space(x_list, y_list):

    import numpy
    
    fy_list = numpy.fft.fft(y_list)
    dx = x_list[1] - x_list[0]
    k_list = 2*numpy.pi*numpy.fft.fftfreq(len(y_list),d=dx)
    return k_list, fy_list

def apply_filter(x_list, y_list, filter_func):

    k_list, fy_list = to_k_space(x_list, y_list)
    filter_list = [filter_func(k) for k in k_list]
    return numpy.fft.ifft([f*fy for f,fy in zip(filter_list, fy_list)])

def afd_advance_1(x_list, y_list, v, t, cfl=0.3):

    dx = x_list[1] - x_list[0]
    
    def filter_func(k):

        import cmath

        temp = 1 - cfl*(1.0-cmath.exp(-1j*k*dx*sgn(v)))
        return temp**(t*abs(v)/(cfl*dx))

    return [x.real for x in apply_filter(x_list, y_list, filter_func)]

def afd_advance_2(x_list, y_list, v, t, cfl=0.3):

    dx = x_list[1] - x_list[0]

    def filter_func(k):

        import cmath

        temp = 0.25*(4-cfl**2+cfl**2*cmath.cos(2*k*dx))-1j*cfl*cmath.sin(k*dx)
        return temp**(t*v/(cfl*dx))

    return [x.real for x in apply_filter(x_list, y_list, filter_func)]

def exact_advance(x_list, y_list, v, t):

    import cmath

    dx = x_list[1] - x_list[0]

    def filter_func(k):

        return cmath.exp(-1j*k*v*t)

    return [x.real for x in apply_filter(x_list, y_list, filter_func)]

def calc_propagation_speeds(initial):

    import math

    g = initial['adiabatic index']
    sound_speed = math.sqrt(g*initial['ambient']['pressure']/initial['ambient']['density'])
    return {'positive riemann invariant':initial['ambient']['velocity']+sound_speed,
            'negative riemann invariant':initial['ambient']['velocity']-sound_speed,
            'entropy':initial['ambient']['velocity']}

def time_advance(initial, time, scheme):

    initial_conserved = primitives_to_conserved(initial)
    propagation_speeds = calc_propagation_speeds(initial)
    final_conserved = {}
    for field in initial_conserved:
        final_conserved[field] = scheme(initial['grid'],
                                        initial_conserved[field],
                                        propagation_speeds[field],
                                        time)
    return conserved_to_primitives(initial, final_conserved)
    
def exact_time_advance(initial, time):

    return time_advance(initial, time, exact_advance)

def first_order_time_advance(initial, time, cfl=0.3):

    def my_scheme(x_list ,y_list, v, t):

        return afd_advance_1(x_list ,y_list, v, t, cfl)

    return time_advance(initial, time, my_scheme)

def second_order_time_advance(initial, time, cfl=0.3):

    def my_scheme(x_list ,y_list, v, t):

        return afd_advance_2(x_list ,y_list, v, t, cfl)

    return time_advance(initial, time, my_scheme)
