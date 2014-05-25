"""
Based on
http://www.grc.nasa.gov/WWW/k-12/airplane/oblique.html
"""

def arcsec(x):

    import math

    return math.acos(1/x)

def arccot(x):

    import math

    return math.atan(1/x)

def critical_shock_angle(m,g):

    """
    Calculates the shock angle at which the transition from weak to strong shock occurs
    Input:
    m - Mach number
    g - Adiabatic index
    """

    import math

    return arcsec(2*m*math.sqrt(-(g/(-4 + (1 - 3*g)*m**2 + math.sqrt((1 + g)*(16 + \
                                                                                  8*(-1 + g)*m**2 + (1 + g)*m**4))))))

def calc_deflection_angle(m,g,s):

    """
    Calculates the deflection angle
    Input:
    m - Mach number
    g - Adiabatic index
    s - Shock angle
    """

    import math
    if math.sin(s)==1./m:
        return 0
    cota = math.tan(s)*(0.5*(g+1)*m**2/(m**2*math.sin(s)**2-1)-1)
    return arccot(cota)

def critical_deflection_angle(m,g):

    """
    Calculates the critical deflection angle
    i.e. the angle at the transition from weak to
    strong shock wave
    Input:
    m - Mach number
    g - Adiabatic index
    """

    sc = critical_shock_angle(m,g)
    return calc_deflection_angle(m,g,sc)

def minimum_shock_angle(m):

    """
    Calculates the shock angle
    for which the deflection angle is zero
    Input:
    m - Mach number
    """

    import math

    return math.asin(1/float(m))

class ShockAngleEqn:

    def __init__(self,m,g,a):

        """
        Class constructor
        Input:
        m - Mach number
        g - Adiabatic index
        a - Deflection angle
        """

        self.m = m
        self.g = g
        self.a = a

    def eval(self, s):

        return self.a - \
            calc_deflection_angle(self.m,
                                  self.g,s)

def bisection(f,xl,xh,n):

    fl = f.eval(xl)
    fh = f.eval(xh)
    
    if fl*fh>0:
        raise NameError('Root is not bracketed')
    if fl==0:
        return xl
    if fh==0:
        return xh

    for i in range(n):
        xm = 0.5*(xh+xl);
        fm = f.eval(xm)
        if fm==0:
            return xm
        
        if fl*fm>0:
            fl = fm
            xl = xm
        elif fh*fm:
            fh = fm
            xh = xm
        else:
            raise NameError('Function is not monotonous')

    return xm

def calc_shock_angle(m,g,a):

    eqn = ShockAngleEqn(m,g,a)
    sh = critical_shock_angle(m,g)
    sl = minimum_shock_angle(m)
    return bisection(eqn,sh,sl,20)
