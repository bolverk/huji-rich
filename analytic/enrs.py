"""
Exact newtonian Riemann solver for an ideal gas
"""

def phugo(v2,d1,p1,g):
	"""
	Calculates the pressure along the hugoniot
	Input:
	v2 - Velocity jump
	d1 - Old density
	p1 - Old pressure
	g - Adiabatic index
	"""

	import math

	return (4*p1 + d1*(1 + g)*v2**2 + \
			math.sqrt(d1)*v2*\
			math.sqrt(16*g*p1 + \
					  d1*(1 + g)**2*v2**2))/4.

def vhugo(p2, d1, p1, g):

	"""
	Calculates the velocity along the hugoniot
	Input:
	p2 - New pressure
	d1 - Old density
	p1 - Old pressure
	g - Adiabatic index
	"""
	
	import math
	
	return (math.sqrt(2)*(-p1 + p2))/math.sqrt(d1*((-1 + g)*p1 + (1 + g)*p2))
	
def visen(p2, d1, p1, g):

	"""
	Calculates the velocity along the isentrope
	Input:
	p2 - New pressure
	d1 - Old density
	p1 - Old pressure
	g - Adiabatic index
	"""
	
	import math
	
	return (2*math.sqrt(g)*p1**(1/(2.*g))*(-p1**((-1 + g)/(2.*g)) + \
	p2**((-1 + g)/(2.*g))))/(math.sqrt(d1)*(-1 + g))
	
def vhydro(p2, d1, p1, g):

	"""
	Calculatest the velocity along the hydrodynamic trajectory
	Input:
	p2 - New pressure
	d1 - Old density
	p1 - Old pressure
	g - Adiabatic index
	"""
	
	if p2>p1:
		return vhugo(p2,d1,p1,g)
	else:
		return visen(p2,d1,p1,g)

class Primitive:
	"""
	Primitive variable (density, pressure and velocity)
	"""
	
	def __init__(self, d, p, v):
	
		"""
		Class constructor
		Input:
		d - Density
		p - Pressure
		v - Velocity
		"""
		
		self.Density = d
		self.Pressure = p
		self.Velocity = v
		
		return

	def clone(self):

		"""
		Creates a clone of the struct
		"""

		res = Primitive(self.Density,
				self.Pressure,
				self.Velocity)
		return res

def left_velocity(pnew, prim, g):
	
	"""
	Calculates the velocity on the left side of the interface
	Input:
	pnew - New pressure
	prim - Primitive variables
	g - Adiabatic index
	"""
	
	d1 = prim.Density;
	p1 = prim.Pressure;
	v1 = prim.Velocity;
	return v1 - vhydro(pnew, d1, p1, g)
		
def right_velocity(pnew, prim, g):
	
	"""
	Calculates the velocity on the right side of the interface
	Input:
	p2 - New pressure
	prim - Primitive variables
	g - Adiabatic index
	"""
	
	d1 = prim.Density
	p1 = prim.Pressure
	v1 = prim.Velocity
	return v1 + vhydro(pnew,d1,p1,g)
	
def eval_trans_eqn(p, left, right, g):

	"""
	Evaluates the transcendental equation
	Input:
	p - Pressure at the interface
	left - Primitive variables on the left side
	right - Primitive variables on the right side
	g - Adiabatic index
	"""
	
	return right_velocity(p,right,g) - left_velocity(p,left,g)
	
def find_upper_bracket(f, xl):

	"""
	Increases the range until the root is bracketed
	Input:
	f - Function
	xl - Lower boundary
	"""
	
	max_iter = 10
	iter = 0
	
	fl = f.Eval(xl)
	xr = xl
	fr = f.Eval(xr)
	while fl*fr>=0:
		xr = 2*xr
		fr = f.Eval(xr)
		if iter>max_iter:
			raise NameError('Maximum number of iteration exceeded in find_upper_bracket')
		else:
			iter = iter + 1
	return 2*xr
	
def bisection(f, xl, xr, tol):

	"""
	Solves an equation using the bisection method
	Input:
	f - Function
	lb - Left bound
	rb - Right bound
	tol - Tolerance
	"""
	
	max_iter = 100
	iter = 0
	
	fl = f.Eval(xl)
	if fl==0:
		return xl
	fr = f.Eval(xr)
	if fr==0:
		return xr
	
	if fl*fr>0:
		print('xl = '+str(xl))
		print('xr = '+str(xr))
		print('fl = '+str(fl))
		print('fr = '+str(fr))
		raise NameError("Root not bracketed")
	
	xm = 0.5*(xl+xr)
	while(abs(xl-xr)>tol):
		xm = 0.5*(xl+xr)
		fm = f.Eval(xm)
		if fm==0:
			return xm
		if fm*fl>0:
			xl = xm
		elif fm*fr>0:
			xr = xm
		else:
			raise NameError('Something bad happened. Probably a NaN')
			
		if iter>max_iter:
			raise NameError('Maximum number of iterations exceeded in bisection')
		else:
			iter = iter + 1
	
	return xm
	
class TransEqn:

	"""
	Transcendental equation
	"""
	
	def __init__(self, left, right, g):
	
		"""
		Class constructor
		Input:
		left - Primitive variables on the left side
		right - Primitive variables on the right side
		g - Adiabatic index
		"""
		
		self._left = left
		self._right = right
		self._g = g
		
		return
	
	def Eval(self, p):
	
		"""
		Evaluates the transcendental equation
		Input:
		p - Pressure at the interface
		"""
		
		return eval_trans_eqn(p,self._left, self._right, self._g)
		
def two_rarefaction(left, right, g):

	"""
	Pressure at the interface of a riemann problem in the case of two rarefactions
	Input:
	left - Primitive variables on the left side
	right - Primitive variables on the right side
	g - Adiabatic index
	"""
	
	import math
	
	dl = left.Density
	pl = left.Pressure
	vl = left.Velocity
	dr = right.Density
	pr = right.Pressure
	vr = right.Velocity
	ps = (4**(g/(1 - g))*(2*(math.sqrt(dr*g*pl) + math.sqrt(dl*g*pr)) + \
		(-1 + g)*math.sqrt(pl)*math.sqrt(pr)*(vl - vr))**((2*g)/(-1 + g)))/\
       (g**(g/(-1 + g))*(dr**(1/(2.*g))*math.sqrt(pl) + dl**(1/(2.*g))*math.sqrt(pr))**((2*g)/(-1 + g)))
	
	return ps

def two_strong_shocks(left, right, g):

	"""
	Pressure at the interface of the Riemann problem in the case of two strong shocks
	Input:
	left - Primitive variables on the left side
	right - Primitive variables on the right side
	g - Adiabatic index
	"""
	
	import math
	
	dl = left.Density
	vl = left.Velocity
	dr = right.Density
	vr = right.Velocity
	return (dl*dr*(1 + g)*(vl - vr)**2)/(2.*(math.sqrt(dl) + math.sqrt(dr))**2)
	
def acoustic_rp(left, right, g):

	"""
	Pressure at the interface of a riemann problem in the case of two acoustic waves
	Input:
	left - Primitive variable on the left side
	right - Primitive variable on the right side
	g - Adiabatic index
	"""
	
	import math
	
	dl = left.Density
	pl = left.Pressure
	vl = left.Velocity
	dr = right.Density
	pr = right.Pressure
	vr = right.Velocity
	return (math.sqrt(dl*g*pl)*pr + math.sqrt(dr*g*pr)*(pl + math.sqrt(dl*g*pl)*(vl - vr)))\
			/(math.sqrt(dl*g*pl) + math.sqrt(dr*g*pr))
	
def riemann_solve(left, right, g):

	"""
	Calculates the pressure and velocity at the interface
	Input:
	left - Primitive variables on the left side
	right - Primitive variable on the right side
	g - Adiabatic index
	"""
	
	tol = 1e-6
	
	vlpr = left_velocity(right.Pressure, left, g)
	vrpl = right_velocity(left.Pressure, right, g)
	if min(left.Velocity,vlpr)>max(right.Velocity,vrpl):
		# Two shocks
		eqn = TransEqn(left, right, g)
		p1 = max(left.Pressure, right.Pressure)
		plvr = phugo(left.Velocity-right.Velocity,
			     left.Density,left.Pressure,g)
		prvl = phugo(left.Velocity-right.Velocity,
			     right.Density,right.Pressure,g)
		p2 = max(plvr,prvl)
		ps = bisection(eqn,p1,p2,tol)
	elif max(left.Velocity,vlpr)<min(right.Velocity,vrpl):
		# Two rarefactions
		ps = two_rarefaction(left, right, g)
	else:
		# Shock rarefaction
		eqn = TransEqn(left, right, g)
		ps = bisection(eqn, min(left.Pressure, right.Pressure), max(left.Pressure, right.Pressure), tol)
	vs = 0.5*(left_velocity(ps,left,g)+right_velocity(ps,right,g))
	return ps, vs
	
# Self similar spatial profiles

# shocks

def shock_speed(p2, d1, p1, g):

	"""
	Calculates the shock speed in the upstream reference frame
	Input:
	p2 - Downstream pressure
	d1 - Upstream density
	p1 - Upstream pressure
	g - Adiabatic index
	"""
	
	import math
	
	return math.sqrt(-p1 + g*p1 + p2 + g*p2)/(math.sqrt(2)*math.sqrt(d1))
	
def dhugo(p2,d1,p1,g):

	"""
	Calculates the density on the hugoniot curve
	Input:
	p2 - Downstream pressure
	d1 - Upstream density
	p1 - Upstream pressure
	g - Adiabatic index
	"""
	
	return (d1*((-1 + g)*p1 + (1 + g)*p2))/((1 + g)*p1 + (-1 + g)*p2)

class ShockProf:
	"""
	Spatial profile of a shock wave
	"""

	def __init__(self, p2, d1, p1, g):
		"""
		Class constructor
		Input:
		p2 - Downstream pressure
		d1 - Upstream density
		p1 - Upstream pressure
		g - Adiabatic index
		"""
		
		self._vs = shock_speed(p2,d1,p1,g)
		self._upstream = Primitive(d1, p1, 0)
		self._downstream = Primitive(\
			dhugo(p2,d1,p1,g),\
				p2, \
				vhugo(p2,d1,p1,g))

		return
		
	def CalcPrim(self, v):
	
		"""
		Retrives the primitive variables
		Input:
		v - Dimensionless coordinate (v=x/t)
		"""
		
		if(v>self._vs):
			return Primitive(self._upstream.Density,\
							self._upstream.Pressure,\
							self._upstream.Velocity)
		else:
			return Primitive(self._downstream.Density,\
							self._downstream.Pressure,\
							self._downstream.Velocity)
			
# Isentropes

def disen(p2,d1,p1,g):

	"""
	Density on the isentrope
	Input:
	p2 - Downstream pressure
	d1 - Upstream density
	p1 - Upstream pressure
	g - Adiabatic index
	"""
	
	return d1*(p2/p1)**(1./g)
	
def visen(p2,d1,p1,g):

	"""
	Calculates the velocity on the hugoniot
	Input:
	p2 - Downstream pressure
	d1 - Upstream density
	p1 - Upstream pressure
	g - Adiabatic index
	"""
	
	import math
	
	return (2*math.sqrt(g)*p1**(1/(2.*g))*(-p1**((-1 + g)/(2.*g)) +\
			p2**((-1 + g)/(2.*g))))/(math.sqrt(d1)*(-1 + g))
	
def sound_speed(p, d, g):

	"""
	Calculates the speed of sound
	Input:
	p - Pressure
	d - Density
	g - Adiabatic index
	"""
	
	import math
	
	return math.sqrt(g*p/d)
	
class IsenProf:
	"""
	Spatial profile of a rarefaction wave
	"""
	
	def __init__(self,p2,d1,p1,g):
		"""
		Class constructor
	
		"""
		
		self._p1 = p1
		self._d1 = d1
		self._g = g
		self._upstream = Primitive(d1,p1,0)
		self._csmax = sound_speed(p1,d1,g)
		d2 = disen(p2,d1,p1,g)
		v2 = visen(p2,d1,p1,g)
		self._csmin = sound_speed(p2,d2,g)+v2
		self._downstream = Primitive(d2,p2,v2)
	
		return
	
	def CalcPrim(self, v):
		
		"""
		Calculates the primitive variables
		Input:
		v - Dimensionless coordinate (x/t)
		"""
		
		if v > self._csmax:
			return Primitive(self._upstream.Density,\
							self._upstream.Pressure,\
							self._upstream.Velocity)
		elif v < self._csmin:
			return Primitive(self._downstream.Density,\
							self._downstream.Pressure,\
							self._downstream.Velocity)
		else:

			u = (v-self._csmax)*2/(1+self._g)
			c = self._csmax+(self._g-1)*u/2;
			p = self._p1*(c/self._csmax)**(2*self._g/(self._g-1))
			d = disen(p,self._d1,self._p1,self._g)
			return Primitive(d,p,u)
			
class HydroProf:
	"""
	Spatial profile for both shock and rarefaction
	"""
	
	def __init__(self, p2, d1, p1, g):
		
		"""
		Class constructor
		Input:
		p2 - Downstream pressure
		d1 - Upstream density
		p1 - Upstream pressure
		g - Adiabatic index
		"""
		
		if p2>p1:
			self._prof = ShockProf(p2,d1,p1,g)
		else:
			self._prof = IsenProf(p2,d1,p1,g)
		
	def CalcPrim(self, v):
		
		"""
		Calculates the primitive variables
		"""
		
		return self._prof.CalcPrim(v)
		
class RiemannProfile:
	
	"""
	Complete Riemann problem profile
	"""
	
	def __init__(self, left, right, g):
		
		"""
		Class constructor
		Input:
		left - Primitive variables on the left side
		right - Primitive variables on the right side
		g - Adiabatic index
		"""
		
		self._ps, self._vs = riemann_solve(left, right, g)
		self._left_prof = HydroProf(self._ps,left.Density,\
						    left.Pressure, g)
		self._right_prof = HydroProf(self._ps,right.Density,\
						     right.Pressure,g)
		self._left = left
		self._right = right
		
		return
	
	def CalcPrim(self, v):
		
		"""
		Calculates the primitive variables
		Input:
		v - Self similar coordinate (x/t)
		"""
		
		if v>self._vs:
			res = self._right_prof.CalcPrim(v - self._right.Velocity)
			res.Velocity = res.Velocity + self._right.Velocity
		else:
			res = self._left_prof.CalcPrim(self._left.Velocity - v)			
			res.Velocity = self._left.Velocity-res.Velocity
		return res
