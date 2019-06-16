#include "LagrangianHLLC.hpp"
#include "hydrodynamic_variables.hpp"
#include "../../misc/universal_error.hpp"
#include "../../misc/utils.hpp"
#include <cmath>

using namespace std;

namespace 
{
	std::pair<double, double> HLLpu(Primitive const& left,Primitive const& right,EquationOfState const& eos)
	{
		double sl = std::min(left.Velocity.x - left.SoundSpeed, right.Velocity.x - right.SoundSpeed);
		double sr = std::max(left.Velocity.x + left.SoundSpeed, right.Velocity.x + right.SoundSpeed);
		Conserved Ul(left.Density, left.Density*left.Velocity, left.Density*(left.Energy + ScalarProd(left.Velocity, left.Velocity*0.5))),
			Ur(right.Density, right.Density*right.Velocity, right.Density*(right.Energy + ScalarProd(right.Velocity, right.Velocity*0.5)));
		Conserved Fl(left.Density*left.Velocity.x, Vector2D(left.Density*left.Velocity.x*left.Velocity.x + left.Pressure, left.Density*left.Velocity.x*left.Velocity.y), (Ul.Energy + left.Pressure)*left.Velocity.x),
			Fr(right.Density*right.Velocity.x, Vector2D(right.Density*right.Velocity.x*right.Velocity.x + right.Pressure, right.Density*right.Velocity.x*right.Velocity.y), (Ur.Energy + right.Pressure)*right.Velocity.x);
		Conserved Ull = (sr*Ur - sl * Ul + Fl - Fr) / (sr - sl);
		return std::pair<double, double> (eos.de2p(Ull.Mass, std::max(0.0,(Ull.Energy - ScalarProd(Ull.Momentum,Ull.Momentum)*0.5 / Ull.Mass))/Ull.Mass), Ull.Momentum.x / Ull.Mass);
	}

	class WaveSpeeds
	{
	public:

		WaveSpeeds(double left_i,
			double center_i,
			double right_i,double ps_i) :
			left(left_i),
			center(center_i),
			right(right_i),ps(ps_i) {}

		WaveSpeeds(WaveSpeeds const& other):left(other.left),center(other.center),
			right(other.right),ps(other.ps){}
		
		WaveSpeeds& operator=(WaveSpeeds const& ws)
		{
			left = ws.left;
			center = ws.center;
			right = ws.right;
			ps = ws.ps;
			return *this;
		}

		double left;
		double center;
		double right;
		double ps;
	};
}

namespace 
{
	WaveSpeeds estimate_wave_speeds(Primitive const& left, Primitive const& right,double pstar)
	{
		const double dl = left.Density;
		const double pl = left.Pressure;
		const double vl = left.Velocity.x;
		const double cl = left.SoundSpeed;
		const double dr = right.Density;
		const double pr = right.Pressure;
		const double vr = right.Velocity.x;
		const double cr = right.SoundSpeed;
		const double sl = vl - cl * (pstar > pl ? std::sqrt(0.8*(pstar / pl - 1) + 1) : 1);
		const double sr = vr + cr * (pstar > pr ? std::sqrt(0.8*(pstar / pr - 1) + 1) : 1);
		const double denom = 1.0 /(dl*(sl - vl) - dr * (sr - vr));
		const double ss = (pr - pl + dl * vl*(sl - vl) - dr * vr*(sr - vr)) *denom;
		const double ps =  std::max(0.0,dl * (sl - vl)*(pr - dr * (vr - vl)*(sr - vr)) *denom - pl * dr*(sr - vr) *denom);
		return WaveSpeeds(sl, ss, sr,ps);
	}
}

namespace 
{
	UniversalError invalid_wave_speeds(Primitive const& left,
		Primitive const& right,
		double velocity,
		double left_wave_speed,
		double center_wave_speed,
		double right_wave_speed)
	{
		UniversalError res("Invalid wave speeds in hllc solver");
		res.AddEntry("left density", left.Density);
		res.AddEntry("left pressure", left.Pressure);
		res.AddEntry("left x velocity", left.Velocity.x);
		res.AddEntry("left y velocity", left.Velocity.y);
		res.AddEntry("left sound speed", left.SoundSpeed);
		res.AddEntry("left energy", left.Energy);
		res.AddEntry("right density", right.Density);
		res.AddEntry("right pressure", right.Pressure);
		res.AddEntry("right x velocity", right.Velocity.x);
		res.AddEntry("right y velocity", right.Velocity.y);
		res.AddEntry("right sound speed", right.SoundSpeed);
		res.AddEntry("right energy", right.Energy);
		res.AddEntry("interface velocity", velocity);
		res.AddEntry("left wave speed", left_wave_speed);
		res.AddEntry("center wave speed", center_wave_speed);
		res.AddEntry("right wave speed", right_wave_speed);
		return res;
	}
}

namespace 
{
	Conserved starred_state
		(Primitive const& state, double sk, double ss)
	{
		const double dk = state.Density;
		const double pk = state.Pressure;
		const double uk = state.Velocity.x;
		const double vk = state.Velocity.y;
		const double ds = dk*(sk - uk) / (sk - ss);
		const double ek = TotalEnergyDensity(state);
		Conserved res;
		res.Mass = ds;
		res.Momentum.x = ds*ss;
		res.Momentum.y = ds*vk;
		res.Energy = ek*ds / dk +
			ds*(ss - uk)*(ss + pk / dk / (sk - uk));
		return res;
	}
}

LagrangianHLLC::LagrangianHLLC(EquationOfState const& eos,bool massflux,bool iter) :eos_(eos), massflux_(massflux),iter_(iter),energy(0)
{}

Conserved LagrangianHLLC::operator()
(Primitive const& left,
	Primitive const& right,
	double velocity) const
{
	if (is_nan(right.Velocity.x))
		throw UniversalError("Hllc::Solved entered with nan");

	const Vector2D normaldir(1, 0);
	Primitive local_left = left;
	Primitive local_right = right;

	local_left.Velocity -= velocity*normaldir;
	local_right.Velocity -= velocity*normaldir;

	Conserved f_gr;
	std::pair<double, double> p_u_star = HLLpu(local_left, local_right, eos_);
	double old_ps = p_u_star.first;

	//old_ps = 0;

	WaveSpeeds ws = estimate_wave_speeds(local_left, local_right, old_ps);
	if (iter_)
	{
		size_t counter = 0;
		while (ws.ps > 1.01 * old_ps || old_ps > 1.01 * ws.ps)
		{
			old_ps = ws.ps;
			ws = estimate_wave_speeds(local_left, local_right, ws.ps);
			++counter;
			if (counter > 54)
			{
				std::cout << "Too many iterations in HLLC" << std::endl;
				std::cout << "Normal " << normaldir.x << "," << normaldir.y << " velocity = " << velocity << std::endl;
				std::cout << " Left density = " << left.Density << " pressure = " << left.Pressure << " internal_energy = " << left.Energy << " vx = " << left.Velocity.x <<
					" vy = " << left.Velocity.y << std::endl;
				std::cout << " Right density = " << right.Density << " pressure = " << right.Pressure << " internal_energy = " << right.Energy << " vx = " << right.Velocity.x <<
					" vy = " << right.Velocity.y << std::endl;
				std::cout << "Old Pstar = " << old_ps << " new Pstar = " << ws.ps << std::endl;
				throw UniversalError("LagrangianHllc::No convergence");
			}

		}
	}
	if (!massflux_)
	{
		local_left.Velocity -= ws.center*normaldir;
		local_right.Velocity -= ws.center*normaldir;
		velocity += ws.center;
		energy = ws.center;
		ws = estimate_wave_speeds(local_left, local_right, ws.ps);
	}


	const Conserved ul = Primitive2Conserved(local_left);
	const Conserved ur = Primitive2Conserved(local_right);

	const Vector2D xdir(1, 0);
	const Conserved fl = Primitive2Flux(local_left, xdir);
	const Conserved fr = Primitive2Flux(local_right, xdir);

	const Conserved usl = starred_state(local_left, ws.left, ws.center);
	const Conserved usr = starred_state(local_right, ws.right, ws.center);

	if (ws.left > 0)
	{
		f_gr = fl;
		if(massflux_)
			energy = ScalarProd(local_left.Velocity, xdir)*local_left.Energy*local_left.Density;
	}
	else if (ws.left <= 0 && ws.center >= 0)
	{
		f_gr = fl + ws.left*(usl - ul);
		if (massflux_)
			energy = ScalarProd(local_left.Velocity, xdir)*local_left.Energy*local_left.Density +
				ws.left*local_left.Energy*((ws.left-local_left.Velocity.x)/(ws.left-ws.center)-1)*local_left.Density;
	}
	else if (ws.center < 0 && ws.right >= 0)
	{
		f_gr = fr + ws.right*(usr - ur);
		if (massflux_)
			energy = (ScalarProd(local_right.Velocity, xdir)*local_right.Energy + ws.right*local_right.Energy*
				((ws.right - local_right.Velocity.x) / (ws.right - ws.center) - 1))*local_right.Density;
	}
	else if (ws.right<0)
	{
		f_gr = fr;
		if (massflux_)
			energy = ScalarProd(local_right.Velocity, xdir)*local_right.Energy*local_right.Density;
	}
	else
		throw invalid_wave_speeds(local_left,
			local_right,
			velocity,
			ws.left,
			ws.center,
			ws.right);

	f_gr.Energy += ScalarProd(f_gr.Momentum, velocity*normaldir) +
		0.5*f_gr.Mass*velocity*velocity;
	f_gr.Momentum += velocity*f_gr.Mass*normaldir;
	return f_gr;
}
