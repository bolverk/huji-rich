#include "Tillotson.hpp"
#include <assert.h>
#include <boost/math/tools/roots.hpp>
#include <iostream>
#include "../../misc/universal_error.hpp"

Tillotson::Tillotson(double a, double b, double A, double B, double rho0, double E0, double EIV, double ECV, double
	alpha, double beta, bool negative_pressure) :
	a_(a), b_(b), A_(A), B_(B), rho0_(rho0), E0_(E0), EIV_(EIV), ECV_(ECV), alpha_(alpha), beta_(beta),
	negative_pressure_(negative_pressure), temp_d_(0), temp_p_(0) {}

double Tillotson::dp2EI(double d, double p) const
{
	double eta = d / rho0_;
	double mu = eta - 1;
	double c = E0_*eta*eta;
	double A = A_*mu;
	double B = B_*mu*mu;
	double sqr = sqrt(4 * a_*c*d*(p - A - B) + std::pow(A + B + (a_ + b_)*c*d - p, 2));
	double first_part = p - A - B - a_*c*d - b_*c*d;
	double E = (first_part + sqr) / (2 * a_*d);
	if (E < 0)
	{
		std::cout << "d " << d << " p " << p << std::endl;
	}
	assert(E > 0);
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
	if (E < 0)
	{
		std::cout << "d " << d << " p " << p << std::endl;
	}
	assert(E > 0);
	return E;
}

double Tillotson::de2pI(double d, double e)const
{
	double eta = d / rho0_;
	double c = E0_*eta*eta;
	double AB = (A_ - 2 * B_)*eta + (B_ - A_) + B_*eta*eta;
	double res = 0;
	if (negative_pressure_)
		res = (a_ + b_ / (e / c + 1))*d*e + AB;
	else
		res = std::max((a_ + b_ / (e / c + 1))*d*e + AB, (a_ + b_)*d*e*1e-7);
	return res;
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
	if (alpha_ > 35 * eta*eta)
		return a_*d*e;
	double mu = eta - 1;
	double c = E0_*eta*eta;
	double A = A_*mu;
	double exp_alpha = std::exp(-alpha_*std::pow(rho0_ / d - 1, 2));
	double exp_beta = A*std::exp(-beta_*(rho0_ / d - 1));
	if (negative_pressure_)
		return a_*d*e + exp_alpha*(b_*d*e / (e / c + 1) + exp_beta);
	else
		return std::max(a_*d*e + exp_alpha*(b_*d*e / (e / c + 1) + exp_beta), (a_ + b_)*d*e*1e-7);
}

double Tillotson::dep2cI(double d, double e, double /*p*/) const
{
	double eta = d / rho0_;
	double w0 = e / (E0_*eta*eta) + 1;
	double mu = eta - 1;
	double res = 2 * B_*mu / rho0_ + A_ / rho0_ + 2 * b_*e*e / (E0_*eta*eta*w0*w0) + (b_ + a_*w0*w0)*(A_*mu + B_*mu*mu + d*e*(a_ + b_ / w0))
		/ (d*w0*w0) + e*(a_+b_/w0);
	res = std::max(res, 1e-10*E0_);
	return res;
}

double Tillotson::dep2cIV(double d, double e, double p) const
{
	double eta = d / rho0_;
	double w0 = e / (E0_*eta*eta) + 1;
	double z = 1.0 / eta - 1.0;
	double afactor = std::exp(-alpha_*z*z);
	double res0 = p*(a_ + b_*afactor / w0 + 1) / d;
	double res1 = A_*std::exp(-beta_*z)*afactor*(1 + (eta - 1)*(beta_ + 2 * alpha_*z - eta) / (eta*eta)) / rho0_;
	res1 += b_*d*e*afactor*(2 * alpha_*z*w0 / rho0_ + (p / d - 2 * e) / (E0_*d)) / (w0*w0*eta*eta);
	double res = std::max(res0 + res1, 1e-10*E0_);
	return res;
}

struct dp2eII
{
	dp2eII(Tillotson const& eos) : eos_(eos)
	{}

	double operator()(double e)
	{
		double res = 1 - eos_.de2pII(eos_.temp_d_, e) / eos_.temp_p_;
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
	if (d >= rho0_)
	{
		return dp2EI(d, p);
	}
	else
	{
		double PIV = de2pI(d, EIV_);
		if (p <= PIV)
		{
			return dp2EI(d, p);
		}
		double exp_alpha = std::exp(-alpha_*std::pow(rho0_ / d - 1, 2));
		double exp_beta = A*std::exp(-beta_*(rho0_ / d - 1));
		double PCV = 0;
		if (negative_pressure_)
			PCV = a_*d*ECV_ + exp_alpha*(b_*d*ECV_ / (ECV_ / c + 1) + exp_beta);
		else
			PCV = std::max(a_*d*ECV_ + exp_alpha*(b_*d*ECV_ / (ECV_ / c + 1) + exp_beta), (a_ + b_)*d*ECV_*1e-7);
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
			try
			{
				res = boost::math::tools::toms748_solve(dp2eII(*this), EIV_, ECV_, boost::math::tools::eps_tolerance<double>(30), it);
			}
			catch (boost::exception const& eo)
			{
				std::cout << " EIV_ " << EIV_ << " ECV_ " << ECV_ << " density " << d << " pressure " << p << " PIV " << PIV << " PCV " << PCV << std::endl;
			}
			double result = 0.5*(res.first + res.second);
			double newp = de2p(d, result);
			if (newp > PIV * 2)
			{
				if (std::abs(p - newp) > 0.001*std::abs(p))
				{
					UniversalError eo("No dp2e convergence");
					eo.AddEntry("Density", d);
					eo.AddEntry("Pressure", p);
					eo.AddEntry("New Pressure", newp);
					eo.AddEntry("EIV", EIV_);
					eo.AddEntry("ECV", ECV_);
					eo.AddEntry("First energy", res.first);
					eo.AddEntry("Second energy", res.second);
					eo.AddEntry("First pressure", de2p(d, res.first));
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
	if (d >= rho0_ || e <= EIV_)
	{
		double res = dep2cI(d, e, p);
		assert(res > 0);
		return std::sqrt(res);
	}
	else
	{
		if (e >= ECV_)
		{
			double res = dep2cIV(d, e, p);
			assert(res > 0);
			return std::sqrt(res);
		}
		else
		{
			double c1 = dep2cI(d, e, de2pI(d,e));
			double c4 = dep2cIV(d, e, de2pIV(d,e));
			double res = std::sqrt((c4*(ECV_ - e) + c1*(e - EIV_)) / (ECV_ - EIV_));
			assert(res > 0);
			return res;
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
