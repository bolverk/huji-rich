#include "Hllc3D.hpp"
#include "../../misc/universal_error.hpp"
#include "../../misc/utils.hpp"

using namespace std;

namespace
{
	class WaveSpeeds
	{
	public:

		WaveSpeeds(double left_i, double center_i, double right_i) :left(left_i), center(center_i), right(right_i) {}
		double left;
		double center;
		double right;
	};

	WaveSpeeds estimate_wave_speeds(ComputationalCell3D const& left, ComputationalCell3D const& right,
		EquationOfState const &eos, TracerStickerNames const& tsn, double gamma)
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
			eo.AddEntry("Error in HLLC3D left sound speed", 0);
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
			eo.AddEntry("Error in HLLC3D right sound speed", 0);
			throw eo;
		}
#endif
		double sl = std::min(vr - cr, vl - cl);
		double sr = std::max(vr + cr, vl + cl);
		double pstar = 0;
		double denom = 1.0 / (dl*(sl - vl) - dr * (sr - vr));
		double ps = std::max(0.0, dl * (sl - vl)*(pr - dr * (vr - vl)*(sr - vr)) *denom - pl * dr*(sr - vr) *denom);
		size_t counter = 0;
		while ((ps > 1.1 * pstar || pstar > 1.1 * ps) && gamma > 0.001)
		{
			pstar = ps;
			double al = (pstar > pl ? fastsqrt(1 + gamma * (pstar / pl - 1)) : 1);
			double ar = (pstar > pr ? fastsqrt(1 + gamma * (pstar / pr - 1)) : 1);
			sl = vl - cl * al;
			sr = vr + cr *ar;
			if (sl > sr)
			{
				sl = std::min(vr - cr*ar, vl - cl*al);
				sr = std::max(vr + cr*ar, vl + cl*al);
			}
			denom = 1.0 / (dl*(sl - vl) - dr * (sr - vr));
			ps = std::max(0.0, dl * (sl - vl)*(pr - dr * (vr - vl)*(sr - vr)) *denom - pl * dr*(sr - vr) *denom);
			++counter;
			if (counter > 54)
			{
				std::cout << "Too many iterations in HLLC" << std::endl;
				//std::cout << "Normal " << normaldir.x << "," << normaldir.y << "," << normaldir.z << " velocity = " << velocity << std::endl;
				std::cout << " Left density = " << left.density << " pressure = " << left.pressure << " internal_energy = " << left.internal_energy << " vx = " << left.velocity.x <<
					" vy = " << left.velocity.y << " vz = " << left.velocity.z << std::endl;
				std::cout << " Right density = " << right.density << " pressure = " << right.pressure << " internal_energy = " << right.internal_energy << " vx = " << right.velocity.x <<
					" vy = " << right.velocity.y << " vz = " << right.velocity.z << std::endl;
			}
			assert(counter < 55);
		}
		double ss = (pr - pl + dl * vl*(sl - vl) - dr * vr*(sr - vr)) *denom;
		return WaveSpeeds(sl, ss, sr);
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

	void BoostBack(Conserved3D &f_gr, double velocity, Vector3D const& normaldir, ComputationalCell3D const& left,
		ComputationalCell3D const& right)
	{
		f_gr.energy += f_gr.momentum.x * velocity + 0.5*f_gr.mass*velocity*velocity;
		f_gr.momentum = (f_gr.momentum.x + f_gr.mass*velocity)*normaldir;
		if (f_gr.mass > 0)
			f_gr.momentum += (left.velocity - normaldir * ScalarProd(left.velocity, normaldir))*f_gr.mass;
		else
			f_gr.momentum += (right.velocity - normaldir * ScalarProd(right.velocity, normaldir))*f_gr.mass;
		f_gr.internal_energy = 0;
	}

	void RotateBack(Conserved3D& f_gr, Vector3D const& normaldir, ComputationalCell3D const& left,
		ComputationalCell3D const& right)
	{
		f_gr.momentum = f_gr.momentum.x * normaldir;
		if (f_gr.mass > 0)
			f_gr.momentum += (left.velocity - normaldir * ScalarProd(left.velocity, normaldir)) * f_gr.mass;
		else
			f_gr.momentum += (right.velocity - normaldir * ScalarProd(right.velocity, normaldir)) * f_gr.mass;
		f_gr.internal_energy = 0;
	}

	Conserved3D HLLstarred_state(ComputationalCell3D const& left, ComputationalCell3D const& right,
		double sl, double sr)
	{
		Conserved3D res;
		Conserved3D Fl = PrimitiveToFlux(left), Fr = PrimitiveToFlux(right);
		double denom = 1.0 / (sr - sl);
		res.mass = (sr * right.density - sl * left.density + Fl.mass - Fr.mass) * denom;
		res.momentum = (sr * right.density * right.velocity - sl * left.density * left.velocity
			+ Fl.momentum - Fr.momentum) * denom;
		res.energy = (sr * TotalEnergyDensity3D(right) - sl * TotalEnergyDensity3D(left)
			+ Fl.energy - Fr.energy) * denom;
		return res;
	}
	
}

Hllc3D::Hllc3D(double gamma) :gamma_((gamma + 1) / (2 * gamma))
{}

Conserved3D Hllc3D::operator()(ComputationalCell3D const& left, ComputationalCell3D const& right, double velocity,
	EquationOfState const& eos, TracerStickerNames const& tsn, Vector3D const& normaldir) const
{
	double face_v = 0;	
	double minv = std::min(fastabs(left.velocity), fastabs(right.velocity));
	bool fast_flow = minv < std::abs(velocity) * 1e-3;
	if (fast_flow)
	{
		face_v = velocity;
		velocity = 0;
	}

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

	Conserved3D ul, ur;
	PrimitiveToConserved(local_left, 1, ul);
	PrimitiveToConserved(local_right, 1, ur);

	WaveSpeeds ws = estimate_wave_speeds(local_left, local_right, eos, tsn, gamma_);

	Conserved3D f_gr;
	// check if bad wavespeed
	bool HLL = false;
	if ((ws.center<ws.left || ws.center>ws.right))
	{
		ws.center = 0.5*(ws.left + ws.right);
		HLL = true;
	}

	if (ws.left > face_v)
	{
		const Conserved3D fl = PrimitiveToFlux(local_left);
		f_gr = fl;
		if (fast_flow)
		{
			Conserved3D ustate;
			PrimitiveToConserved(local_left, 1, ustate);
			f_gr -= ustate * face_v;			
		}
	}
	else if (ws.left <= face_v && ws.center >= face_v)
	{
		const Conserved3D fl = PrimitiveToFlux(local_left);
		const Conserved3D usl = starred_state(local_left, ws.left, ws.center);
		f_gr = fl + ws.left*(usl - ul);
		if (fast_flow)
			f_gr -= usl * face_v;
	}
	else if (ws.center < face_v && ws.right >= face_v)
	{
		const Conserved3D fr = PrimitiveToFlux(local_right);
		const Conserved3D usr = starred_state(local_right, ws.right, ws.center);
		f_gr = fr + ws.right*(usr - ur);
		if (fast_flow)
			f_gr -= usr * face_v;
	}
	else if (ws.right < face_v)
	{
		const Conserved3D fr = PrimitiveToFlux(local_right);
		f_gr = fr;
		if (fast_flow)
		{
			Conserved3D ustate;
			PrimitiveToConserved(local_right, 1, ustate);
			f_gr -= ustate * face_v;
		}
	}
	else
		throw invalid_wave_speeds(local_left, local_right, velocity, ws.left, ws.center, ws.right);

	// check if bad wavespeed
	if (HLL && ws.left < face_v && ws.right>face_v)
	{
		const Conserved3D fr = PrimitiveToFlux(local_right);
		const Conserved3D fl = PrimitiveToFlux(local_left);
		f_gr = (ws.right*fl - ws.left*fr + ws.left*ws.right*(ur - ul))*(1.0 / (ws.right - ws.left)); // HLL flux
		if (fast_flow)
			f_gr -= face_v * HLLstarred_state(local_left, local_right, ws.left, ws.right);
	}
	if (fast_flow)
		RotateBack(f_gr, normaldir, left, right);
	else
		BoostBack(f_gr, velocity, normaldir, left, right);
	return f_gr;
}
