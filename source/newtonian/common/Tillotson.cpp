#include "Tillotson.hpp"
#include <assert.h>
#include <boost/math/tools/roots.hpp>
#include <iostream>
#include "../../misc/universal_error.hpp"

Tillotson::Tillotson(double a, double b, double A, double B, double rho0, double E0, double EIV, double ECV,double
	alpha,double beta):
	a_(a),b_(b),A_(A),B_(B),rho0_(rho0),E0_(E0),EIV_(EIV),ECV_(ECV),alpha_(alpha),beta_(beta),temp_d_(0),temp_p_(0){}

double Tillotson::dp2EI(double d, double p) const
{
	double eta = d / rho0_;
	double mu = eta - 1;
	double c = E0_*eta*eta;
	double A = A_*mu;
	double B = B_*mu*mu;
	double sqr = sqrt(4 * a_*c*d*(p - A - B) + std::pow(A + B + (a_ + b_)*c*d - p, 2));
	double first_part = p - A - B - a_*c*d - b_*c*d;
	double E = (first_part + sqr) / 	(2 * a_*d);
	E = std::max(E, 0.0);
	return E;
}

double Tillotson::dp2EIV(double d, double p) const
{
	double eta = d / rho0_;
	double mu = eta - 1;
	double c = E0_*eta*eta;
	double A = A_*mu;
	double exp_alpha = std::exp(-alpha_*std::pow(rho0_ / d - 1, 2));
	double exp_beta = A*std::exp(-beta_*(rho0_ / d - 1));
	double b = b_*exp_alpha;
	double AB = exp_alpha*exp_beta;
	double E = (p - AB - a_*c*d - b*c*d + sqrt(4 * a_*c*d*(p - AB) + std::pow(AB + (a_ + b)*c*d - p, 2))) /
		(2 * a_*d);
	assert(E > 0);
	return E;
}

double Tillotson::de2pI(double d, double e)const
{
	double eta = d / rho0_;
	double mu = eta - 1;
	double c = E0_*eta*eta;
	double A = A_*mu;
	double B = B_*mu*mu;
	double res = (a_ + b_ / (e / c + 1))*d*e + A + B;
	return std::max(res,(a_+b_)*d*e*1e-8);
}

double Tillotson::de2pII(double d, double e)const
{
	double P2 = de2pI(d, e);
	double P3 = de2pIV(d, e);
	return ((e - EIV_)*P3 + (ECV_ - e)*P2) / (ECV_ - EIV_);
}

double Tillotson::de2pIV(double d, double e)const
{
	double eta = d / rho0_;
	double mu = eta - 1;
	double c = E0_*eta*eta;
	double A = A_*mu;
	double exp_alpha = std::exp(-alpha_*std::pow(rho0_ / d - 1, 2));
	double exp_beta = A*std::exp(-beta_*(rho0_ / d - 1));
	return std::max(a_*d*e + exp_alpha*(b_*d*e / (e / c + 1) + exp_beta), (a_+b_)*d*e*1e-8);
}

double Tillotson::dep2cI(double d, double e, double p) const
{
	double eta = d / rho0_;
	double w0 = e / (E0_*eta*eta) + 1;
	double gamma = a_ + b_ / w0;
	return ((gamma + 1)*p + (A_ + B_*(eta*eta - 1))) / d + b_*(w0 - 1)*(2 * e - p / d) / (w0*w0);
}

double Tillotson::dep2cIV(double d, double e, double p) const
{
	double eta = d / rho0_;
	double mu = eta - 1;
	double w0 = e / (E0_*eta*eta) + 1;
	double z = 1.0 / eta - 1;
	double gamma = a_ + b_*std::exp(-beta_*z*z);
	return (gamma + 1)*p / d + A_*std::exp(-(alpha_*z + beta_*z*z))*(1 + mu*(alpha_ + 2 * beta_*z - eta) / (eta*eta)) / rho0_
		+ b_*d*e*std::exp(-beta_*z*z)*(2 * beta_*z*w0 / rho0_ + (p / d - 2 * e) / (E0_*d)) / (w0*w0*eta*eta);
}

/*struct dp2eII
{ 
	dp2eII(Tillotson const& eos) : eos_(eos)
	{}

	std::pair<double, double> operator()(double e)
	{
		double eta = eos_.temp_d_ / eos_.rho0_;
		double mu = eta - 1;
		double c = eos_.E0_*eta*eta;
		double A = eos_.A_*mu;
		double B = eos_.B_*mu*mu;
		double exp_alpha = std::exp(-eos_.alpha_*std::pow(eos_.rho0_ / eos_.temp_d_ - 1, 2));
		double exp_beta = std::exp(-eos_.beta_*(eos_.rho0_ / eos_.temp_d_ - 1));
		double res = eos_.temp_p_ - eos_.de2pII(eos_.temp_d_, e);
		double slope = (B*std::pow(c + e, 2) - A*(-1 + exp_alpha*exp_beta)*std::pow(c + e, 2) +
			eos_.temp_d_*(-eos_.a_*std::pow(c + e, 2)*(eos_.ECV_ - eos_.EIV_) + eos_.b_* c*(-(-1 + exp_alpha)*e*e +
				c*(-eos_.ECV_ + 2 * e - 2 * exp_alpha*e + exp_alpha*eos_.EIV_)))) / (std::pow(c + e, 2)* (eos_.ECV_ - eos_.EIV_));
		return std::make_pair(res, slope);
	}
private:
	Tillotson const& eos_;                               
};
*/
struct dp2eII
{
	dp2eII(Tillotson const& eos) : eos_(eos)
	{}

	double operator()(double e)
	{
		double res = 1 - eos_.de2pII(eos_.temp_d_, e)/ eos_.temp_p_;
		return res;
	}
private:
	Tillotson const& eos_;
};


double Tillotson::dp2e(double d, double p, tvector const & /*tracers*/, vector<string> const & /*tracernames*/) const
{
	double eta = d / rho0_;
	double mu = eta - 1;
	double c = E0_*eta*eta;
	double A = A_*mu;
	double B = B_*mu*mu;
	if (d >= rho0_)
	{
		return dp2EI(d, p);
	}
	else
	{
		double PIV = std::max((a_ + b_ / (EIV_ / c + 1))*d*EIV_ + A + B, (a_+b_)*d*EIV_*1e-8);
		if (p <= PIV)
		{
			return dp2EI(d, p);
		}
		double exp_alpha = std::exp(-alpha_*std::pow(rho0_ / d - 1, 2));
		double exp_beta = A*std::exp(-beta_*(rho0_ / d - 1));
		double PCV = std::max(a_*d*ECV_ + exp_alpha*(b_*d*ECV_ / (ECV_ / c + 1) + exp_beta), (a_+b_)*d*ECV_*1e-8);
		if (p >= PCV)
		{
			return dp2EIV(d, p);
		}
		else
		{
			temp_d_ = d;
			temp_p_ = p;
			boost::uintmax_t it = 50;
			std::pair<double, double> res;
			if (d * 1000 > rho0_)
			{
				res=boost::math::tools::bisect(dp2eII(*this), EIV_, ECV_,
					boost::math::tools::eps_tolerance<double>(30), it);
			}
			else
			{
				it = 100;
				res = boost::math::tools::bisect(dp2eII(*this), EIV_, ECV_,
					boost::math::tools::eps_tolerance<double>(40), it);
			}
			double result = 0.5*(res.first + res.second);
			double newp = de2p(d, result);
			if (newp > PIV * 2)
			{
				if (std::abs(p - newp) > 0.001*p)
				{
					UniversalError eo("No dp2e convergence");
					eo.AddEntry("Density", d);
					eo.AddEntry("Pressure", p);
					eo.AddEntry("New Pressure", newp);
					eo.AddEntry("EIV", EIV_);
					eo.AddEntry("ECV", ECV_);
					eo.AddEntry("First energy", res.first);
					eo.AddEntry("Second energy", res.second);
					eo.AddEntry("First pressure", de2p(d,res.first));
					eo.AddEntry("Second pressure", de2p(d, res.second));
					throw eo;
				}
			}
			assert(result > 0);
			return result;
		}
	}
}

double Tillotson::de2p(double d, double e, tvector const& /*tracers*/, vector<string> const& /*tracernames*/) const
{
	if (d >= rho0_)
	{
		return de2pI(d, e);
	}
	else
	{
		if (e <= EIV_)
			return de2pI(d, e);
		if (e >= ECV_)
			return de2pIV(d, e);
		return de2pII(d, e);
	}
}

double Tillotson::de2c(double d, double e, tvector const & tracers, vector<string> const & tracernames) const
{
	double p = de2p(d, e, tracers, tracernames);
	return dp2c(d, p, tracers, tracernames);
}

double Tillotson::dp2c(double d, double p, tvector const & /*tracers*/, vector<string> const & /*tracernames*/) const
{
	double e = dp2e(d, p);
	if (d >= rho0_ || e<=EIV_)
	{
		double res = dep2cI(d, e, p);
		assert(res > 0);
		return sqrt(res);
	}
	else
	{
		if (e >= ECV_)
		{
			double res = dep2cIV(d, e, p);
			assert(res > 0);
			return sqrt(res);
		}
		else
		{
			double res = (dep2cIV(d, e, p)*(ECV_-e)+ dep2cI(d, e, p)*(e-EIV_))/(ECV_-EIV_);
			assert(res > 0);
			return sqrt(res);
		}
	}
}

double Tillotson::dp2s(double /*d*/, double /*p*/, tvector const & /*tracers*/, vector<string> const & /*tracernames*/) const
{
	std::cout << "dp2s not implemented" << std::endl;
	assert(false);
	return 0;
}

double Tillotson::sd2p(double /*s*/, double /*d*/, tvector const & /*tracers*/, vector<string> const & /*tracernames*/) const
{
	std::cout << "sd2p not implemented" << std::endl;
	assert(false);
	return 0;
}

Tillotson::~Tillotson(void)
{}
