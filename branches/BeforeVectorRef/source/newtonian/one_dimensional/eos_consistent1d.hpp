/*! \file eos_consistent1d.hpp
  \brief A class that guarantees that a spatial reconstruction would satisfy the equation of state
  \author Almog Yalinewich
 */

#ifndef EOSCONSISTENT_HPP
#define EOSCONSISTENT_HPP 1

#include "spatial_reconstruction1d.hpp"
#include "../common/equation_of_state.hpp"

using std::vector;

//! \brief Interpolations for one dimensional simulation
namespace interpolations1d{

  //! \brief A class that guarantees that a spatial reconstruction would satisfy the equation of state
  class EOSConsistent: public SpatialReconstruction1D
  {
  public:

    /*! \brief Class constructor
      \param naive Base interpolation
      \param eos Equation of state
    */
    EOSConsistent(SpatialReconstruction1D const& naive,
		  EquationOfState const& eos);

    Primitive InterpState(vector<double> const& vp,
			  vector<Primitive> const& hv,
			  double interface_speed,
			  size_t i, int dir, double dt) const;

  private:
    SpatialReconstruction1D const& naive_;
    EquationOfState const& eos_;
  };
}

#endif // EOSCONSISTENT_HPP
