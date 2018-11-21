#include "Hllc_SR.hpp"
#include <cmath>
#include "../newtonian/common/hydrodynamic_variables.hpp"
#include "../misc/universal_error.hpp"
#include "../misc/utils.hpp"

Hllc_SR::Hllc_SR()
{}



using namespace std;

namespace {
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
	WaveSpeeds estimate_wave_speeds(Primitive const& left, Primitive const& right,double &pstar)
	{
		const double sig_left = left.SoundSpeed*left.SoundSpeed*(1 - ScalarProd(left.Velocity, left.Velocity)) / (1 - left.SoundSpeed*left.SoundSpeed);
		double disct = std::sqrt(sig_left*(1 - left.Velocity.x*left.Velocity.x + sig_left));
		const double lamda_minus_left = (left.Velocity.x - disct) / (1 + sig_left);
		const double lamda_plus_left = (left.Velocity.x + disct) / (1 + sig_left);
		const double sig_right = right.SoundSpeed*right.SoundSpeed*(1 - ScalarProd(right.Velocity, right.Velocity)) / (1 - right.SoundSpeed*right.SoundSpeed);
		disct = std::sqrt(sig_right*(1 - right.Velocity.x*right.Velocity.x + sig_right));
		const double lamda_minus_right = (right.Velocity.x - disct) / (1 + sig_right);
		const double lamda_plus_right = (right.Velocity.x + disct) / (1 + sig_right);
		const double sl = std::min(lamda_minus_right, lamda_minus_left);
		const double sr = std::max(lamda_plus_left, lamda_plus_right);
		// calculate hll flux components
		const double El = (left.Energy+1)*left.Density / (1 - ScalarProd(left.Velocity, left.Velocity)) - left.Pressure;
		const double Er = (right.Energy+1)*right.Density / (1 - ScalarProd(right.Velocity, right.Velocity)) - right.Pressure;
		const double mxl = (El + left.Pressure)*left.Velocity.x;
		const double mxr = (Er + right.Pressure)*right.Velocity.x;;
		const double E_hll = (sr*Er - sl * El + mxl - mxr) / (sr - sl);
		const double mx_hll = (sr*mxr - sl * mxl + (mxl*left.Velocity.x + left.Pressure) - (mxr*right.Velocity.x + right.Pressure)) / (sr - sl);
		const double F_E = (sr * mxl - sl * mxr + sr * sl*(Er - El)) / (sr - sl);
		const double F_mx = (sr * (mxl*left.Velocity.x + left.Pressure) - sl * (mxr*right.Velocity.x + right.Pressure) + sr * sl*(mxr - mxl)) / (sr - sl);
		const double b = -(E_hll + F_mx);
		const double ss = (std::abs(F_E) < std::max(1e-14*E_hll, 1e-14)) ? -mx_hll / b : (-b - std::sqrt(b*b - 4 * mx_hll*F_E)) / (2 * F_E);
		//const double ss= (1e6*std::abs(F_E) < std::abs(E_hll+F_mx))	? mx_hll/(E_hll+F_mx) :	(((E_hll + F_mx) - std::sqrt((E_hll + F_mx)*(E_hll + F_mx) - 4 * mx_hll*F_E)) / (2 * F_E));
		pstar = -F_E * ss + F_mx;
		assert(sr > sl);
		assert(sr > ss);
		assert(ss > sl);
		return WaveSpeeds(sl, ss, sr);
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
		res.AddEntry("left Energy", (left.Energy+1));
		res.AddEntry("right density", right.Density);
		res.AddEntry("right pressure", right.Pressure);
		res.AddEntry("right x velocity", right.Velocity.x);
		res.AddEntry("right y velocity", right.Velocity.y);
		res.AddEntry("right sound speed", right.SoundSpeed);
		res.AddEntry("right Energy", (right.Energy+1));
		res.AddEntry("interface velocity", velocity);
		res.AddEntry("left wave speed", left_wave_speed);
		res.AddEntry("center wave speed", center_wave_speed);
		res.AddEntry("right wave speed", right_wave_speed);
		return res;
	}
}

namespace 
{
	Conserved SR_Primitive2Flux(Primitive const& p,double edge_vel)
	{
		const double g_2 = 1.0 / (1 - ScalarProd(p.Velocity, p.Velocity));
		double gamma =std::sqrt(g_2);
		return Conserved(p.Density*gamma*(p.Velocity.x-edge_vel),Vector2D(p.Density*(p.Energy+1)*g_2*p.Velocity.x*(p.Velocity.x-edge_vel) + p.Pressure, p.Density*(p.Energy+1)*g_2*p.Velocity.y*(p.Velocity.x-edge_vel)), 
			p.Density*(p.Energy+1)*g_2*p.Velocity.x-edge_vel*(p.Density*(p.Energy+1)*g_2-p.Pressure));
	}	

	Conserved starred_flux(Primitive const& state, double lambda_star, double lambda,double edge_vel,double pstar)
	{
		const double g_2 = 1.0/(1-ScalarProd(state.Velocity,state.Velocity));
		const double dlambda_1 = 1.0 / (lambda - lambda_star);
		const double d = state.Density*(lambda-state.Velocity.x)*dlambda_1*std::sqrt(g_2);
		const double mx = (state.Density*(state.Energy+1)*state.Velocity.x*g_2*(lambda - state.Velocity.x) + pstar - state.Pressure)*dlambda_1;
		const double my = (state.Density*(state.Energy+1)*state.Velocity.y*g_2*(lambda - state.Velocity.x))*dlambda_1;
		const double E = ((state.Density*(state.Energy+1)*g_2 - state.Pressure)*(lambda - state.Velocity.x) + pstar * lambda_star - state.Pressure*state.Velocity.x)*dlambda_1;
		Conserved res(d*(lambda_star-edge_vel), Vector2D(mx*(lambda_star-edge_vel) + pstar, my*(lambda_star-edge_vel)), mx - E*edge_vel);
		return res;
	}

/*	Conserved Hll_flux(Primitive const& left, Primitive const& right, double ll, double lr)
	{
		Conserved Fl = SR_Primitive2Flux(left,0);
		Conserved Fr = SR_Primitive2Flux(right,0);
		double g_2l = 1.0 / (1 - ScalarProd(left.Velocity, left.Velocity));
		Conserved Ul(left.Density*std::sqrt(g_2l), g_2l*left.Density*(left.Energy + 1)*left.Velocity, g_2l*left.Density*(left.Energy + 1)-left.Pressure);
		double g_2r = 1.0 / (1 - ScalarProd(right.Velocity, right.Velocity));
		Conserved Ur(right.Density*std::sqrt(g_2r), g_2r*right.Density*(right.Energy + 1)*right.Velocity, g_2r*right.Density*(right.Energy + 1)-right.Pressure);
		return (lr*Fl - ll * Fr + ll * lr*(Ur - Ul)) / (lr - ll);
	}*/
}

Conserved Hllc_SR::operator()(Primitive const& left,Primitive const& right,	double velocity) const
{
	if (is_nan(right.Velocity.x))
		throw UniversalError("Hllc_SR::Solved entered with nan");
	double pstar = 0;
	WaveSpeeds ws = estimate_wave_speeds(left,right,pstar);
	
	Conserved f_gr;
	if (ws.left > velocity)
		f_gr = SR_Primitive2Flux(left, velocity);
	else if (ws.left <= velocity && ws.center >= velocity)
		f_gr = starred_flux(left, ws.center, ws.left, velocity,pstar);
		//f_gr = Hll_flux(left, right, ws.left, ws.right);
	else if (ws.center < velocity &&ws.right >= velocity)
		f_gr = starred_flux(right, ws.center, ws.right, velocity,pstar);
		//f_gr = Hll_flux(left, right, ws.left, ws.right);
	else if (ws.right < velocity)
		f_gr = SR_Primitive2Flux(right, velocity);
	else
		throw invalid_wave_speeds(left,
			right,
			velocity,
			ws.left,
			ws.center,
			ws.right);
	// Remove rest mass Energy
	f_gr.Energy -= f_gr.Mass;
	return f_gr;
}
