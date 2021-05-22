/*! \brief Tillotson EOS
\file Tillotson.hpp
\author Elad Steinberg
*/

#ifndef TILLOTSON_HPP
#define TILLOTSON_HPP 1
#include "equation_of_state.hpp"

//! \brief Helper class
struct dp2eII;

//! \brief Tillotson equation of state
class Tillotson : public EquationOfState
{
  //! \brief Auxiliary struct
	friend  struct dp2eII;
private:
	double a_, b_, A_, B_, rho0_, E0_, EIV_, ECV_,alpha_,beta_;
	const bool negative_pressure_;
	mutable double temp_d_, temp_p_;
	const size_t e_index_;

	double dp2EI(double d, double p)const;
	double dp2EIV(double d, double p)const;
	double de2pI(double d, double e)const;
	double de2pII(double d, double e)const;
	double de2pIV(double d, double e)const;
	double dep2cI(double d, double e, double p)const;
	double dep2cIV(double d, double e, double p)const;
public:
  /*! \brief Class constructor
    \param a a parameter
    \param b b parameter
    \param A A parameter
    \param B B parameter
    \param rho0 Refernce density
    \param E0 Refernce energy
    \param EIV EIV parameter
    \param ECV ECV parameter
    \param alpha alpha Parameter
    \param beta beta Parameter
    \param negative_pressure Flag for handling negative pressure
    \param e_index TBD
   */
	Tillotson(double a, double b, double A, double B, double rho0, double E0, double EIV, double ECV,double alpha,
		double beta,bool negative_pressure=false, size_t e_index = 0);

	 double dp2e(double d, double p, tvector const& tracers = tvector(), 
		 vector<string> const& tracernames = vector<string>())const override;

	double de2p(double d, double e,	tvector const& tracers = tvector(), 
		vector<string> const& tracernames = vector<string>()) const override;

	double de2c(double d, double e,	tvector const& tracers = tvector(), 
		vector<string> const& tracernames = vector<string>()) const override;
	
	double dp2c(double d, double p,tvector const& tracers = tvector(), 
		vector<string> const& tracernames = vector<string>()) const override;

	double dp2s(double d, double p,	tvector const& tracers = tvector(), 
		vector<string> const& tracernames = vector<string>()) const override;

	double sd2p(double s, double d,tvector const& tracers = tvector(), 
		vector<string> const& tracernames = vector<string>()) const override;

	~Tillotson(void) override;
};

#endif //TILLOTSON_HPP

