#ifndef EOSCONSISTENT_HPP
#define EOSCONSISTENT_HPP 1

#include "spatial_reconstruction1d.hpp"
#include "../common/equation_of_state.hpp"

namespace interpolations1d{
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
			  int i, int dir, double dt) const;

  private:
    SpatialReconstruction1D const& naive_;
    EquationOfState const& eos_;
  };
}

#endif // EOSCONSISTENT_HPP
