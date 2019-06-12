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
		EquationOfState const &eos, TracerStickerNames const& tsn)
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
		const double sl = std::min(vl - cl, vr - cr);
		const double sr = std::max(vl + cl, vr + cr);
		const double ss = (pr - pl + dl * vl*(sl - vl) - dr * vr*(sr - vr)) /
			(dl*(sl - vl) - dr * (sr - vr));
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

}

Conserved3D Hllc3D::operator()(ComputationalCell3D const& left, ComputationalCell3D const& right, double velocity,
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

	Conserved3D ul, ur;
	PrimitiveToConserved(local_left, 1, ul);
	PrimitiveToConserved(local_right, 1, ur);

	const Conserved3D fl = PrimitiveToFlux(local_left);
	const Conserved3D fr = PrimitiveToFlux(local_right);

	WaveSpeeds ws = estimate_wave_speeds(local_left, local_right, eos, tsn);

	const Conserved3D usl = starred_state(local_left, ws.left, ws.center);
	const Conserved3D usr = starred_state(local_right, ws.right, ws.center);

	Conserved3D f_gr;
	// check if bad wavespeed
	bool HLL = false;
	if ((ws.center<ws.left || ws.center>ws.right))
	{
		ws.center = 0.5*(ws.left + ws.right);
		HLL = true;
	}

	if (ws.left > 0)
		f_gr = fl;
	else if (ws.left <= 0 && ws.center >= 0)
		f_gr = fl + ws.left*(usl - ul);
	else if (ws.center < 0 && ws.right >= 0)
		f_gr = fr + ws.right*(usr - ur);
	else if (ws.right < 0)
		f_gr = fr;
	else
		throw invalid_wave_speeds(local_left, local_right, velocity, ws.left, ws.center, ws.right);

	// check if bad wavespeed
	if (HLL && ws.left < 0 && ws.right>0)
		f_gr = ws.right*fl - ws.left*fr + ws.left*ws.right*(ur - ul)*(1.0 / (ws.right - ws.left)); // HLL flux

	f_gr.energy += f_gr.momentum.x * velocity + 0.5*f_gr.mass*velocity*velocity;
	f_gr.momentum = (f_gr.momentum.x + f_gr.mass*velocity)*normaldir;
	if (f_gr.mass > 0)
		f_gr.momentum += (left.velocity - normaldir * ScalarProd(left.velocity, normaldir))*f_gr.mass;
	else
		f_gr.momentum += (right.velocity - normaldir * ScalarProd(right.velocity, normaldir))*f_gr.mass;

	f_gr.internal_energy = 0;
	return f_gr;
}
