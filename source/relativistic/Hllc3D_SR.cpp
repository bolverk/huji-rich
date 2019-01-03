#include "Hllc3D_SR.hpp"

Hllc3D_SR::Hllc3D_SR() : local_left_(ComputationalCell3D()), local_right_(ComputationalCell3D())
{}


Hllc3D_SR::~Hllc3D_SR()
{}

using namespace std;

namespace
{
	class WaveSpeeds
	{
	public:

		WaveSpeeds(double left_i,
			double center_i,
			double right_i) :
			left(left_i),
			center(center_i),
			right(right_i) {}

		const double left;
		const double center;
		const double right;
	};
}

namespace
{
	WaveSpeeds estimate_wave_speeds(ComputationalCell3D const& left, ComputationalCell3D const& right,
		EquationOfState const& eos, TracerStickerNames const& tsn, double &pstar)
	{
		double cl = 0, cr = 0;
#ifdef RICH_DEBUG
		try
		{
#endif
			cl = eos.dp2c(left.density, left.pressure, left.tracers, tsn.tracer_names);
#ifdef RICH_DEBUG
		}
		catch (UniversalError &eo)
		{
			eo.AddEntry("Error in HLLC3D_SR left sound speed", 0);
			throw eo;
		}
#endif
		const double sig_left = cl * cl*(1 - ScalarProd(left.velocity, left.velocity)) / (1 - cl * cl);
		double disct = std::sqrt(sig_left*(1 - left.velocity.x*left.velocity.x + sig_left));
		const double lamda_minus_left = (left.velocity.x - disct) / (1 + sig_left);
		const double lamda_plus_left = (left.velocity.x + disct) / (1 + sig_left);
#ifdef RICH_DEBUG
		try
		{
#endif
			cr = eos.dp2c(right.density, right.pressure, right.tracers, tsn.tracer_names);
#ifdef RICH_DEBUG
		}
		catch (UniversalError &eo)
		{
			eo.AddEntry("Error in HLLC3D_SR right sound speed", 0);
			throw eo;
		}
#endif
		const double sig_right = cr * cr*(1 - ScalarProd(right.velocity, right.velocity)) / (1 - cr * cr);
		disct = std::sqrt(sig_right*(1 - right.velocity.x*right.velocity.x + sig_right));
		const double lamda_minus_right = (right.velocity.x - disct) / (1 + sig_right);
		const double lamda_plus_right = (right.velocity.x + disct) / (1 + sig_right);
		const double sl = std::min(lamda_minus_right, lamda_minus_left);
		const double sr = std::max(lamda_plus_left, lamda_plus_right);
		// calculate hll flux components
		const double El = (left.internal_energy + 1)*left.density / (1 - ScalarProd(left.velocity, left.velocity)) - left.pressure;
		const double Er = (right.internal_energy + 1)*right.density / (1 - ScalarProd(right.velocity, right.velocity)) - right.pressure;
		const double mxl = (El + left.pressure)*left.velocity.x;
		const double mxr = (Er + right.pressure)*right.velocity.x;;
		const double E_hll = (sr*Er - sl * El + mxl - mxr) / (sr - sl);
		const double mx_hll = (sr*mxr - sl * mxl + (mxl*left.velocity.x + left.pressure) - (mxr*right.velocity.x + right.pressure)) / (sr - sl);
		const double F_E = (sr * mxl - sl * mxr + sr * sl*(Er - El)) / (sr - sl);
		const double F_mx = (sr * (mxl*left.velocity.x + left.pressure) - sl * (mxr*right.velocity.x + right.pressure) + sr * sl*(mxr - mxl)) / (sr - sl);
		const double b = -(E_hll + F_mx);
		const double ss = (std::abs(F_E) < std::max(1e-7*E_hll, 1e-30)) ? -mx_hll / b : (-b - std::sqrt(b*b - 4 * mx_hll*F_E)) / (2 * F_E);
		pstar = -F_E * ss + F_mx;
		assert(sr > sl);
		assert(sr > ss);
		assert(ss > sl);
		return WaveSpeeds(sl, ss, sr);
	}
}

namespace
{
	UniversalError invalid_wave_speeds(ComputationalCell3D const& left,
		ComputationalCell3D const& right,
		double velocity,
		double left_wave_speed,
		double center_wave_speed,
		double right_wave_speed)
	{
		UniversalError res("Invalid wave speeds in hllc solver");
		res.AddEntry("left density", left.density);
		res.AddEntry("left pressure", left.pressure);
		res.AddEntry("left x velocity", left.velocity.x);
		res.AddEntry("left y velocity", left.velocity.y);
		res.AddEntry("left z velocity", left.velocity.z);
		res.AddEntry("left Energy", (left.internal_energy + 1));
		res.AddEntry("right density", right.density);
		res.AddEntry("right pressure", right.pressure);
		res.AddEntry("right x velocity", right.velocity.x);
		res.AddEntry("right y velocity", right.velocity.y);
		res.AddEntry("right z velocity", right.velocity.z);
		res.AddEntry("right Energy", (right.internal_energy + 1));
		res.AddEntry("interface velocity", velocity);
		res.AddEntry("left wave speed", left_wave_speed);
		res.AddEntry("center wave speed", center_wave_speed);
		res.AddEntry("right wave speed", right_wave_speed);
		return res;
	}
}

namespace
{
	Conserved3D SR_Primitive2Flux(ComputationalCell3D const& p, double edge_vel)
	{
		const double g_2 = 1.0 / (1 - ScalarProd(p.velocity, p.velocity));
		double gamma = std::sqrt(g_2);
		return Conserved3D(p.density*gamma*(p.velocity.x - edge_vel), Vector3D(p.pressure, 0, 0) + p.density*(p.internal_energy + 1)*g_2*p.velocity*(p.velocity.x - edge_vel),
			p.density*p.internal_energy*g_2*p.velocity.x + p.density*p.velocity.x*(g_2 - gamma) - edge_vel * (p.density*p.internal_energy*g_2 - p.pressure + p.density*(g_2 - gamma)), 0);
	}

	Conserved3D starred_flux(ComputationalCell3D const& state, double lambda_star, double lambda, double edge_vel, double pstar)
	{
		const double g_2 = 1.0 / (1 - ScalarProd(state.velocity, state.velocity));
		const double dlambda_1 = 1.0 / (lambda - lambda_star);
		const double d = state.density*(lambda - state.velocity.x)*dlambda_1*std::sqrt(g_2);
		const double mx = (state.density*(state.internal_energy + 1)*state.velocity.x*g_2*(lambda - state.velocity.x) + pstar - state.pressure)*dlambda_1;
		//const double my = (state.density*(state.internal_energy + 1)*state.velocity.y*g_2*(lambda - state.velocity.x))*dlambda_1;
		//const double mz = (state.density*(state.internal_energy + 1)*state.velocity.z*g_2*(lambda - state.velocity.x))*dlambda_1;
		const double E = ((state.density*(state.internal_energy + 1)*g_2 - state.pressure)*(lambda - state.velocity.x) + pstar * lambda_star - state.pressure*state.velocity.x)*dlambda_1;
		Conserved3D res(d*(lambda_star - edge_vel), Vector3D(mx*(lambda_star - edge_vel) + pstar, 0, 0), (lambda_star - edge_vel)*(E - d) + pstar * lambda_star, 0);
		return res;
	}
}


Conserved3D Hllc3D_SR::operator()(ComputationalCell3D const & left, ComputationalCell3D const & right, double velocity, EquationOfState const & eos, TracerStickerNames const & tsn, Vector3D const & normaldir) const
{
	ReplaceComputationalCell(local_left_, left);
	double par_left = fastabs(local_left_.velocity - ScalarProd(local_left_.velocity, normaldir)*normaldir);
	ReplaceComputationalCell(local_right_, right);
	double par_right = fastabs(local_right_.velocity - ScalarProd(local_right_.velocity, normaldir)*normaldir);
	local_left_.velocity.x = ScalarProd(local_left_.velocity, normaldir);
	local_left_.velocity.y = par_left;
	local_left_.velocity.z = 0;
	local_right_.velocity.x = ScalarProd(local_right_.velocity, normaldir);
	local_right_.velocity.y = par_right;
	local_right_.velocity.z = 0;
	double pstar = 0;
	WaveSpeeds ws = estimate_wave_speeds(local_left_, local_right_, eos, tsn, pstar);

	Conserved3D f_gr;
	if (ws.left > velocity)
		f_gr = SR_Primitive2Flux(local_left_, velocity);
	else if (ws.left <= velocity && ws.center >= velocity)
		f_gr = starred_flux(local_left_, ws.center, ws.left, velocity, pstar);
	else if (ws.center < velocity &&ws.right >= velocity)
		f_gr = starred_flux(local_right_, ws.center, ws.right, velocity, pstar);
	else if (ws.right < velocity)
		f_gr = SR_Primitive2Flux(local_right_, velocity);
	else
		throw invalid_wave_speeds(local_left_, local_right_, velocity, ws.left, ws.center, ws.right);

	f_gr.momentum = f_gr.momentum.x * normaldir;
	if (ws.center >= velocity)
		f_gr.momentum += (left.velocity - normaldir * ScalarProd(left.velocity, normaldir))*f_gr.mass*(left.internal_energy + 1)*std::sqrt(1.0 / (1 - ScalarProd(left.velocity, left.velocity)));
	else
		f_gr.momentum += (right.velocity - normaldir * ScalarProd(right.velocity, normaldir))*f_gr.mass*(right.internal_energy + 1)*std::sqrt(1.0 / (1 - ScalarProd(right.velocity, right.velocity)));
#ifdef RICH_DEBUG
	bool good = true;
	if (!std::isfinite(f_gr.energy))
		good = false;
	if (!std::isfinite(f_gr.internal_energy))
		good = false;
	if (!std::isfinite(f_gr.momentum.x))
		good = false;
	if (!std::isfinite(f_gr.momentum.y))
		good = false;
	if (!std::isfinite(f_gr.momentum.z))
		good = false;
	if (!std::isfinite(f_gr.mass))
		good = false;
	if (!good)
	{
		UniversalError eo("Bad flux");
		eo.AddEntry("Energy flux", f_gr.energy);
		eo.AddEntry("Internal Energy flux", f_gr.internal_energy);
		eo.AddEntry("Mass flux", f_gr.mass);
		eo.AddEntry("Momentum x flux", f_gr.momentum.x);
		eo.AddEntry("Momentum y flux", f_gr.momentum.y);
		eo.AddEntry("Momentum z flux", f_gr.momentum.z);
		eo.AddEntry("Left cell density", left.density);
		eo.AddEntry("Left cell pressure", left.pressure);
		eo.AddEntry("Left cell Vx", left.velocity.x);
		eo.AddEntry("Left cell Vy", left.velocity.y);
		eo.AddEntry("Left cell Vz", left.velocity.z);
		eo.AddEntry("Left cell internal energy", left.internal_energy);
		eo.AddEntry("Left cell ID", left.ID);
		eo.AddEntry("Right cell density", right.density);
		eo.AddEntry("Right cell pressure", right.pressure);
		eo.AddEntry("Right cell Vx", right.velocity.x);
		eo.AddEntry("Right cell Vy", right.velocity.y);
		eo.AddEntry("Right cell Vz", right.velocity.z);
		eo.AddEntry("Right cell internal energy", right.internal_energy);
		eo.AddEntry("Right cell ID", right.ID);
		eo.AddEntry("Face velocity", velocity);
		throw eo;
	}
#endif

	return f_gr;
}
