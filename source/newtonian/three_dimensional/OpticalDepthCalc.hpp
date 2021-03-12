#ifndef OPTDEPTHCALC_HPP
#define OPTDEPTHCALC_HPP 1

#include "../../3D/GeometryCommon/Tessellation3D.hpp"
#include "computational_cell.hpp"
#include <string>

//! \brief Optical depth calculator
class OpticalDepthCalc
{
private:
	const double opening_;
#ifdef RICH_MPI
	Tessellation3D const* tproc_;
#endif // RICH_MPI
	std::string const d_name_;
  OpticalDepthCalc(const OpticalDepthCalc&/*other*/) :opening_(0)
#ifdef RICH_MPI
						     , tproc_(0)
#endif // RICH_MPI
						     , d_name_("") {}
	OpticalDepthCalc& operator=(const OpticalDepthCalc /*other*/) 
  { 
#ifdef RICH_MPI
    this->tproc_=0;
#endif // RICH_MPI
    return *this; 
  }
public:
	OpticalDepthCalc(double opening = 0.25
			 #ifdef RICH_MPI
			 , Tessellation3D const* tproc = 0
#endif // RICH_MPI
			 , const std::string& debug_name = "");

	// returns dz * Sigma, so for time need to multiply by kappa and divide by c
	void operator()(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells,
		std::vector<double>& res) const;
	~OpticalDepthCalc() {}
};

#endif //OPTDEPTHCALC_HPP

