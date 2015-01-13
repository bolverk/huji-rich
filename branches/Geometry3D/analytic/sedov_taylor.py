def vtoz(v,w,g,n):
	"""
	Calculates the dimensionless radius
	Input:
	v - Dimensionless velocity
	w - Ambient density power law index
	g - Adiabatic index
	n - Geometry (1 - plane, 2 - cylinderical, 3 - spherical)
	"""
	
	return (4.**(1./(2. + n - w))*((-1. + g)/(1. + g))**((-1. + g)/(2. - n + g*(-2. + \
		w)))*(2. + n - (2. + (-1. + g)*n)*v - w)**((g**2.*(n**2. + (-2. + w)**2.) + \
		2.*(-2. + n)*(2.*n - w) + g*(-4. - 3.*n**2. - 2.*n*(-4. + w) + w**2.))/((2. + \
		(-1. + g)*n)*(2. - n + g*(-2. + w))*(2. + n - w))))/(((1. + g)*v)**(2./(2. + \
		n - w))*(-1. + g*v)**((-1. + g)/(2. - n + g*(-2. + w)))*(-((2. - 3.*n + w + \
		g*(-2. + n + w))/(1. + g)))**((g**2.*(n**2. + (-2. + w)**2.) + 2.*(-2. + \
		n)*(2.*n - w) + g*(-4. - 3.*n**2. - 2.*n*(-4. + w) + w**2.))/((2. + (-1. + \
		g)*n)*(2. - n + g*(-2. + w))*(2. + n - w))))
		
def vtod(v,w,g,n):
	"""
	Calculates the dimensionless density
	Input:
	v - Dimensionless velocity
	w - Ambient density power law index
	g - Adiabatic index
	n - Geometry (1 - plane, 2 - cylinderical, 3 - spherical)
	"""

	return (((-1 + g)/(1 + g))**(-1 + (n - g*w)/(2 - n + g*(-2 + w)) + (-2*n + \
		w + g*w)/((-2 + g)*n + w))*(1 + g)**((2*w)/(2 + n - w))*(1 - \
		v)**((2*n - (1 + g)*w)/((-2 + g)*n + w))*v**((2*w)/(2 + n - w))*(-1 + \
		g*v)**((n - g*w)/(-2 + n - g*(-2 + w)))*(2 + n - (2 + (-1 + g)*n)*v - \
		w)**((-(g**3*n*(n**2 + (-2 + w)**2)*w) + 2*(-2 + n)*(2*n**3 - \
		n**2*(-4 + w) - 6*n*w + 2*w**2) + g*(-3*n**4 + n**3*(2 - 6*w) - \
		2*w*(-4 + w**2) + 2*n*(-4 - 8*w + w**2) + n**2*(12 + 10*w + 3*w**2)) \
		+ g**2*(n**4 - 2*(-2 + w)**2*w + n**3*(2 + 3*w) + n**2*(4 - 14*w + \
		3*w**2) - n*(-8 + 4*w - 2*w**2 + w**3)))/((2 + (-1 + g)*n)*(2 + n - \
		w)*(-((-2 + g)*n**2) + n*(-4 - 2*g*(-3 + w) + g**2*(-2 + w) - w) + (2 \
		+ g*(-2 + w))*w))))/(4**(w/(2 + n - w))*(-((2 - 3*n + w + g*(-2 + n + \
		w))/(1 + g)))**((-(g**3*n*(n**2 + (-2 + w)**2)*w) + 2*(-2 + \
		n)*(2*n**3 - n**2*(-4 + w) - 6*n*w + 2*w**2) + g*(-3*n**4 + n**3*(2 - \
		6*w) - 2*w*(-4 + w**2) + 2*n*(-4 - 8*w + w**2) + n**2*(12 + 10*w + \
		3*w**2)) + g**2*(n**4 - 2*(-2 + w)**2*w + n**3*(2 + 3*w) + n**2*(4 - \
		14*w + 3*w**2) - n*(-8 + 4*w - 2*w**2 + w**3)))/((2 + (-1 + g)*n)*(2 \
		+ n - w)*(-((-2 + g)*n**2) + n*(-4 - 2*g*(-3 + w) + g**2*(-2 + w) - \
		w) + (2 + g*(-2 + w))*w))))
		
def vtop(v,w,g,n):
	"""
	Calculates the dimensionless pressure
	Input:
	v - Dimensionless velocity
	w - Ambient density power law index
	g - Adiabatic index
	n - Geometry (1 - plane, 2 - cylinderical, 3 - spherical)
	"""
	
	return ((1. + g)**((g*n*(2. - 3*n + w + g*(-2. + n + w)))/((2. + (-1. + \
		g)*n)*((-2. + g)*n + w)))*(1. - 2./(1. + g))**((-2.*n + w + g*w)/((-2. + \
		g)*n + w))*(1. - v)**((g*(n - w))/((-2. + g)*n + w))*v**((2.*n)/(2. + n - \
		w))*(-2. + 3*n - w - g*(-2. + n + w))**((n*(g**2.*(n**2. + (-2. + w)**2.) + \
		2.*(-2. + n)*(2.*n - w) + g*(-4 - 3*n**2. - 2.*n*(-4 + w) + w**2.)))/((2. + \
		(-1. + g)*n)*(2. + n - w)*((-2. + g)*n + w))))/(2.**((-2. + n + w)/(2. + n \
		- w))*(-1. + g)*(2. + n - (2. + (-1. + g)*n)*v - w)**((n*(g**2.*(n**2. + \
		(-2. + w)**2.) + 2.*(-2. + n)*(2.*n - w) + g*(-4 - 3*n**2. - 2.*n*(-4 + w) + \
		w**2.)))/((2. + (-1. + g)*n)*(2. + n - w)*((-2. + g)*n + w))))
