#include "LagrangianHLLC3D.hpp"
#include <iostream>

using namespace std;

namespace
{
	std::pair<double, double> HLLpu(ComputationalCell3D const& left, ComputationalCell3D const& right, EquationOfState const& eos)
	{
		double cl = 0, cr = 0;
#ifdef RICH_DEBUG
		try
		{
#endif
			cl = eos.dp2c(left.density, left.pressure);
#ifdef RICH_DEBUG
		}
		catch (UniversalError &eo)
		{
			eo.AddEntry("Error in HLLpu left sound speed", 0);
			throw eo;
		}
#endif
#ifdef RICH_DEBUG
		try
		{
#endif
			cr = eos.dp2c(right.density, right.pressure);
#ifdef RICH_DEBUG
		}
		catch (UniversalError &eo)
		{
			eo.AddEntry("Error in HLLpu right sound speed", 0);
			throw eo;
		}
#endif
		double sl = std::min(left.velocity.x - cl, right.velocity.x - cr);
		double sr = std::max(left.velocity.x + cl, right.velocity.x + cr);
		Conserved3D Ul(left.density, left.density*left.velocity, left.density*(left.internal_energy + ScalarProd(left.velocity, left.velocity)*0.5), left.internal_energy),
			Ur(right.density, right.density*right.velocity, right.density*(right.internal_energy + ScalarProd(right.velocity, right.velocity)*0.5), right.internal_energy);
		Conserved3D Fl(left.density*left.velocity.x, Vector3D(left.density*left.velocity.x*left.velocity.x + left.pressure, left.density*left.velocity.x*left.velocity.y,
			left.density*left.velocity.x*left.velocity.z), (Ul.energy + left.pressure)*left.velocity.x, 0),
			Fr(right.density*right.velocity.x, Vector3D(right.density*right.velocity.x*right.velocity.x + right.pressure, right.density*right.velocity.x*right.velocity.y,
				right.density*right.velocity.x*right.velocity.z), (Ur.energy + right.pressure)*right.velocity.x, 0);
		Conserved3D Ull = (sr*Ur - sl * Ul + Fl - Fr) / (sr - sl);
		double ustar = Ull.momentum.x / Ull.mass;
		return std::pair<double, double>(eos.de2p(Ull.mass, std::max(Ull.mass*ustar*ustar*1e-16,Ull.energy - ScalarProd(Ull.momentum, Ull.momentum)*0.5 / Ull.mass)), ustar);
	}

	class WaveSpeeds
	{
	public:

		WaveSpeeds(double left_i,
			double center_i,
			double right_i, double ps_i) :
			left(left_i),
			center(center_i),
			right(right_i), ps(ps_i) {}

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

	WaveSpeeds estimate_wave_speeds(ComputationalCell3D const& left, ComputationalCell3D const& right,
		EquationOfState const &eos, TracerStickerNames const& tsn, double pstar)
	{
		double cl = 0, cr = 0;
		const double dl = left.density;
		const double pl = left.pressure;
		const double vl = left.velocity.x;
#ifdef RICH_DEBUG
		try
		{
#endif
			cl = eos.dp2c(dl, pl, left.tracers, tsn.tracer_names);
#ifdef RICH_DEBUG
		}
		catch (UniversalError &eo)
		{
			eo.AddEntry("Error in LagHLLC3D left sound speed", 0);
			throw eo;
		}
#endif
		const double dr = right.density;
		const double pr = right.pressure;
		const double vr = right.velocity.x;
#ifdef RICH_DEBUG
		try
		{
#endif
			cr = eos.dp2c(dr, pr, right.tracers, tsn.tracer_names);
#ifdef RICH_DEBUG
		}
		catch (UniversalError &eo)
		{
			eo.AddEntry("Error in LagHLLC3D right sound speed", 0);
			throw eo;
		}
#endif
		const double sl = vl - cl * (pstar > pl ? std::sqrt(0.8*(pstar / pl - 1) + 1) : 1);
		const double sr = vr + cr * (pstar > pr ? std::sqrt(0.8*(pstar / pr - 1) + 1) : 1);
		const double denom = 1.0 / (dl*(sl - vl) - dr * (sr - vr));
		const double ss = (pr - pl + dl * vl*(sl - vl) - dr * vr*(sr - vr)) *denom;
		const double ps = std::max(0.0,dl * (sl - vl)*(pr - dr * (vr - vl)*(sr - vr)) *denom - pl * dr*(sr - vr) *denom);
		return WaveSpeeds(sl, ss, sr, ps);
	}

	UniversalError invalid_wave_speeds(ComputationalCell3D const& left, ComputationalCell3D const& right,
		double velocity, double left_wave_speed, double center_wave_speed, double right_wave_speed)
	{
		UniversalError res("Invalid wave speeds in hllc solver");
		res.AddEntry("left density", left.density);
		res.AddEntry("left pressure", left.pressure);
		res.AddEntry("left x velocity", left.velocity.x);
		res.AddEntry("left y velocity", left.velocity.y);
		res.AddEntry("right density", right.density);
		res.AddEntry("right pressure", right.pressure);
		res.AddEntry("right x velocity", right.velocity.x);
		res.AddEntry("right y velocity", right.velocity.y);
		res.AddEntry("interface velocity", velocity);
		res.AddEntry("left wave speed", left_wave_speed);
		res.AddEntry("center wave speed", center_wave_speed);
		res.AddEntry("right wave speed", right_wave_speed);
		return res;
	}

	double TotalEnergyDensity3D(ComputationalCell3D const& p)
	{
		return p.density*(0.5*ScalarProd(p.velocity, p.velocity) + p.internal_energy);
	}

	Conserved3D starred_state(ComputationalCell3D const& state, double sk, double ss)
	{
		const double dk = state.density;
		const double pk = state.pressure;
		const double uk = state.velocity.x;
		const double vk = state.velocity.y;
		const double wk = state.velocity.z;
		const double ds = dk * (sk - uk) / (sk - ss);
		const double ek = TotalEnergyDensity3D(state);
		Conserved3D res;
		res.mass = ds;
		res.momentum.x = ds * ss;
		res.momentum.y = ds * vk;
		res.momentum.z = ds * wk;
		res.energy = ek * ds / dk + ds * (ss - uk)*(ss + pk / dk / (sk - uk));
		res.internal_energy = 0;
		return res;
	}

	Conserved3D PrimitiveToFlux(ComputationalCell3D const& p)
	{
		return Conserved3D(p.density*p.velocity.x,
			p.pressure*Vector3D(1, 0, 0) + p.density*p.velocity.x*p.velocity,
			(TotalEnergyDensity3D(p) + p.pressure)*	p.velocity.x, 0);
	}

}

LagrangianHLLC3D::LagrangianHLLC3D(bool mass_flux) :massflux_(mass_flux), ws(0) {}

Conserved3D LagrangianHLLC3D::operator()(ComputationalCell3D const& left, ComputationalCell3D const& right, double velocity,
	EquationOfState const& eos, TracerStickerNames const& tsn, Vector3D const& normaldir) const
{

	ComputationalCell3D local_left = left;
	ComputationalCell3D local_right = right;

	local_left.velocity.x -= velocity * normaldir.x;
	local_left.velocity.y -= velocity * normaldir.y;
	local_left.velocity.z -= velocity * normaldir.z;
	local_right.velocity.x -= velocity * normaldir.x;
	local_right.velocity.y -= velocity * normaldir.y;
	local_right.velocity.z -= velocity * normaldir.z;

	double par_left = fastabs(local_left.velocity - ScalarProd(local_left.velocity, normaldir)*normaldir);
	double par_right = fastabs(local_right.velocity - ScalarProd(local_right.velocity, normaldir)*normaldir);
	local_left.velocity.x = ScalarProd(local_left.velocity, normaldir);
	local_left.velocity.y = par_left;
	local_left.velocity.z = 0;
	local_right.velocity.x = ScalarProd(local_right.velocity, normaldir);
	local_right.velocity.y = par_right;
	local_right.velocity.z = 0;
	Conserved3D f_gr;
	std::pair<double, double> p_u_star = HLLpu(local_left, local_right, eos);
	p_u_star.first = std::max(p_u_star.first, 0.0);
	WaveSpeeds ws2 = estimate_wave_speeds(local_left, local_right, eos, tsn, p_u_star.first);
	double old_ps = ws2.ps;
	ws2 = estimate_wave_speeds(local_left, local_right, eos, tsn, ws2.ps);
	size_t counter = 0;
	while (ws2.ps > 1.1 * old_ps || old_ps > 1.1 * ws2.ps)
	{
		old_ps = ws2.ps;
		ws2 = estimate_wave_speeds(local_left, local_right, eos, tsn, ws2.ps);
		++counter;
		if (counter > 54)
		{
			std::cout << "Too many iterations in HLLC" << std::endl;
			std::cout << "Normal " << normaldir.x << "," << normaldir.y << "," << normaldir.z << " velocity = " << velocity << std::endl;
			std::cout << " Left density = " << left.density << " pressure = " << left.pressure << " internal_energy = " << left.internal_energy << " vx = " << left.velocity.x <<
				" vy = " << left.velocity.y << " vz = " << left.velocity.z << std::endl;
			std::cout << " Right density = " << right.density << " pressure = " << right.pressure << " internal_energy = " << right.internal_energy << " vx = " << right.velocity.x <<
				" vy = " << right.velocity.y << " vz = " << right.velocity.z << std::endl;
		}
		assert(counter < 55);
	}

	if (!massflux_)
	{
		ws = ws2.center;
		f_gr.energy = ws2.ps*ws;
		f_gr.momentum.Set(ws2.ps, 0, 0);
		f_gr.mass = 0;
	}
	else
	{
		Conserved3D ul, ur;
		PrimitiveToConserved(local_left, 1, ul);
		PrimitiveToConserved(local_right, 1, ur);

		const Conserved3D fl = PrimitiveToFlux(local_left);
		const Conserved3D fr = PrimitiveToFlux(local_right);

		const Conserved3D usl = starred_state(local_left, ws2.left, ws2.center);
		const Conserved3D usr = starred_state(local_right, ws2.right, ws2.center);


		if (ws2.left > 0)
			f_gr = fl;
		else if (ws2.left <= 0 && ws2.center >= 0)
			f_gr = fl + ws2.left*(usl - ul);
		else if (ws2.center < 0 && ws2.right >= 0)
			f_gr = fr + ws2.right*(usr - ur);
		else if (ws2.right < 0)
			f_gr = fr;
		else
			throw invalid_wave_speeds(local_left, local_right, velocity, ws2.left, ws2.center, ws2.right);
	}
	f_gr.energy += f_gr.momentum.x * velocity + 0.5*f_gr.mass*velocity*velocity;
	f_gr.momentum = (f_gr.momentum.x + f_gr.mass*velocity)*normaldir;
	if (f_gr.mass > 0)
		f_gr.momentum += (left.velocity - normaldir * ScalarProd(left.velocity, normaldir))*f_gr.mass;
	else
		f_gr.momentum += (right.velocity - normaldir * ScalarProd(right.velocity, normaldir))*f_gr.mass;
	f_gr.internal_energy = 0;
	return f_gr;
}
