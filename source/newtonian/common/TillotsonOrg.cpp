#include "TillotsonOrg.hpp"
#include <assert.h>
#include <boost/math/tools/roots.hpp>
#include <iostream>
#include "../../misc/universal_error.hpp"
#include "../../misc/utils.hpp"

TillotsonOrg::TillotsonOrg(double a, double b, double A, double B, double rho0, double E0, double EIV, double ECV, double
	alpha, double beta):
	a_(a), b_(b), A_(A), B_(B), rho0_(rho0), E0_(E0), EIV_(EIV), ECV_(ECV), alpha_(alpha), beta_(beta),
	temp_d_(0), temp_p_(0) {}

double TillotsonOrg::dp2EI(double d, double p) const
{
	double eta = d / rho0_;
	double c = E0_ * eta*eta;
	double AB = (A_ - 2 * B_)*eta + (B_ - A_) + B_ * eta*eta;
	double A_B = A_ - B_;
	double sqr = std::sqrt(d*c * 2 * (a_ - b_)*(p - (A_B + eta * B_)*(eta - 1)) + d * d*c*c*(a_ + b_)*(a_ + b_) + (p - (eta - 1)*(A_B + eta * B_))
		*(p - (eta - 1)*(A_B + eta * B_)));
	double first_part = p - AB - a_ * c*d - b_ * c*d;
	double E = (first_part + sqr) / (2 * a_*d);
	assert(E > 0);
	double newp = de2pI(d, E);
	if (newp > 1.00001*p || newp < 0.99999*p)
		E = p * 1e7 / (a_*d);
	
	return E;
}

double TillotsonOrg::dp2EIV(double d, double p) const
{
	double eta = d / rho0_;
	double mu = eta - 1;
	double c = E0_ * eta*eta;
	double A = A_ * mu;
	double eta2 = alpha_ * (rho0_ / d - 1) * (rho0_ / d - 1);
	double exp_alpha = (eta2 > 100) ? 0 : std::exp(-eta2);
	double eta3 = beta_ * (rho0_ / d - 1);
	double exp_beta = (eta3 > 100) ? 0 : A * std::exp(-eta3);
	double b = b_ * exp_alpha;
	double AB = exp_alpha * exp_beta;
	double E = (p - AB - a_ * c*d - b * c*d + sqrt(4 * a_*c*d*(p - AB) + std::pow(AB + (a_ + b)*c*d - p, 2))) /
		(2 * a_*d);
	assert(E > 0);
	double newp = de2pIV(d, E);
	if (newp > 1.00001*p || newp < 0.99999*p)
		E = p * 1e5 / (a_*d);
	return E;
}

double TillotsonOrg::de2pI(double d, double e)const
{
	double eta = d / rho0_;
	double c = E0_ * eta*eta;
	double AB = (A_ - 2 * B_)*eta + (B_ - A_) + B_ * eta*eta;
	double res = 0;
	res = std::max((a_ + b_ / (e / c + 1))*d*e + AB, a_*d*e*1e-7);
	return res;
}

double TillotsonOrg::de2pII(double d, double e)const
{
	double P2 = de2pI(d, e);
	double P3 = de2pIV(d, e);
	return ((e - EIV_)*P3 + (ECV_ - e)*P2) / (ECV_ - EIV_);
}

double TillotsonOrg::de2pIV(double d, double e)const
{
	double eta = d / rho0_;
	if (alpha_ > 100 * eta*eta)
		return a_ * d*e;
	double mu = eta - 1;
	double c = E0_ * eta*eta;
	double A = A_ * mu;
	double eta2 = alpha_ * (rho0_ / d - 1) * (rho0_ / d - 1);
	double exp_alpha = (eta2 > 100) ? 0 : std::exp(-eta2);
	double eta3 = beta_ * (rho0_ / d - 1);
	double exp_beta = (eta3 > 100) ? 0 : A * std::exp(-eta3);
	return std::max(a_*d*e + exp_alpha * (b_*d*e / (e / c + 1) + exp_beta), a_ *d*e*1e-5);
}

double TillotsonOrg::dep2cI(double d, double e, double p) const
{
	double eta = d / rho0_;
	double c = E0_ * eta*eta;
	double w0 = e / c + 1;
	double gamma = a_ + b_ / w0;
	double AB = (A_ - 2 * B_)*eta + (B_ - A_) + B_ * eta*eta;
	p = (a_ + b_ / (e / c + 1))*d*e + AB;
	double res = (gamma + 1)*p / d + (A_ + B_ * (eta*eta - 1)) / d + b_ * (w0 - 1)*(2 * e - p / d) / (w0*w0);
	res = std::max(res, 1e-10*E0_);
	return res;
}

double TillotsonOrg::dep2cIV(double d, double e, double p) const
{
	double eta = d / rho0_;
	double w0 = e / (E0_*eta*eta) + 1;
	double z = 1.0 / eta - 1.0;
	double afactor = (alpha_*z*z > 100) ? 0 : std::exp(-alpha_ * z*z);
	double res0 = p * (a_ + b_ * afactor / w0 + 1) / d;
	double bfactor = (beta_*z > 100) ? 0 : std::exp(-beta_ * z);
	double res1 = A_ * bfactor*afactor*(1 + (eta - 1)*(beta_ + 2 * alpha_*z - eta) / (eta*eta)) / rho0_;
	res1 += b_ * d*e*afactor*(2 * alpha_*z*w0 / rho0_ + (p / d - 2 * e) / (E0_*d)) / (w0*w0*eta*eta);
	double res = std::max(res0 + res1, 1e-10*E0_);
	return res;
}

	struct dp2eIIOrg
	{
		dp2eIIOrg(TillotsonOrg const& eos) : eos_(eos)
		{}

		double operator()(double e)
		{
			double res = 1 - eos_.de2pII(eos_.temp_d_, e) / eos_.temp_p_;
			return res;
		}
	private:
		TillotsonOrg const& eos_;
	};



double TillotsonOrg::dp2e(double d, double p, tvector const & /*tracers*/, vector<string> const & /*tracernames*/) const
{
	double eta = d / rho0_;
	double mu = eta - 1;
	double c = E0_ * eta*eta;
	double A = A_ * mu;
	if (d >= rho0_)
	{
		return dp2EI(d, p);
	}
	else
	{
		double e1 = dp2EI(d, p);
		if (e1<EIV_)
			return e1;
		double e4 = dp2EIV(d, p);
		if (e4 > ECV_)
			return e4;
		double PIV = de2pI(d, EIV_);
		double eta2 = alpha_ * (rho0_ / d - 1) * (rho0_ / d - 1);
		double exp_alpha = (eta2 > 100) ? 0 : std::exp(-eta2);
		double eta3 = beta_ * (rho0_ / d - 1);
		double exp_beta = (eta3 > 100) ? 0 : A * std::exp(-eta3);
		double PCV = 0;
		PCV = std::max(a_*d*ECV_ + exp_alpha * (b_*d*ECV_ / (ECV_ / c + 1) + exp_beta), (a_ + b_)*d*ECV_*1e-5);
		temp_d_ = d;
		temp_p_ = p;
		boost::uintmax_t it = 50;
		std::pair<double, double> res;
		try
		{
			res = boost::math::tools::toms748_solve(dp2eIIOrg(*this), EIV_, ECV_, boost::math::tools::eps_tolerance<double>(30), it);
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

double TillotsonOrg::de2p(double d, double e, tvector const& /*tracers*/, vector<string> const& /*tracernames*/) const
{
	if (d >= rho0_)
	{
		return de2pI(d, e);
	}
	else
	{
		if (e <= EIV_ )
			return de2pI(d, e);
		if (e >= ECV_)
			return de2pIV(d, e);
		return de2pII(d, e);
	}
}

double TillotsonOrg::de2c(double d, double e, tvector const & tracers, vector<string> const & tracernames) const
{
	double p = de2p(d, e, tracers, tracernames);
	return dp2c(d, p, tracers, tracernames);
}

double TillotsonOrg::dp2c(double d, double p, tvector const & tracers, vector<string> const & tracernames) const
{
	double e = dp2e(d, p, tracers, tracernames);
	if (d >= rho0_  || e<EIV_)
	{
		double res = dep2cI(d, e, p);
		assert(res > 0);
		return std::sqrt(res);
	}
	else
	{
		if (e >= ECV_ )
		{
			double res = dep2cIV(d, e, p);
			assert(res > 0);
			return std::sqrt(res);
		}
		else
		{
			double c1 = dep2cI(d, e, p);
			double c4 = dep2cIV(d, e, p);
			double res = std::sqrt((c1*(ECV_ - e) + c4 * (e - EIV_)) / (ECV_ - EIV_));
			assert(res > 0);
			return res;
		}
	}
}

double TillotsonOrg::dp2s(double /*d*/, double /*p*/, tvector const & /*tracers*/, vector<string> const & /*tracernames*/) const
{
	std::cout << "dp2s not implemented" << std::endl;
	assert(false);
	return 0;
}

double TillotsonOrg::sd2p(double /*s*/, double /*d*/, tvector const & /*tracers*/, vector<string> const & /*tracernames*/) const
{
	std::cout << "sd2p not implemented" << std::endl;
	assert(false);
	return 0;
}

TillotsonOrg::~TillotsonOrg(void)
{}

