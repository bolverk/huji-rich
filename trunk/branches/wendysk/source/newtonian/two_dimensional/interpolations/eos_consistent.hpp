#ifndef EOS_CONSISTENT_HPP
#define EOS_CONSISTENT_HPP 1

#include "../spatial_reconstruction.hpp"
#include "../../common/equation_of_state.hpp"
#include "../HydroBoundaryConditions.hpp"

/*! \brief A wrapper for spatial reconstruction that guarantees that the interpolated values will satisfy the equation of state
  \details The class uses the interpolated values of the density and pressure and calculates the energy and sound speed from the equation of state
 */
class EOSConsistent: public SpatialReconstruction
{
public:
	/*!
	\brief Class constructor
	\param sr The original spatial reconstruction
	\param eos The equation of state
	*/
  EOSConsistent(SpatialReconstruction& sr,
		EquationOfState const& eos);

  void Prepare(Tessellation const* tess,vector<Primitive> const& cells,
	       double dt,vector<bool> const& mask,double time);

  Primitive Interpolate(Tessellation const* tess,
			vector<Primitive> const& cells,
			double dt,Edge const& edge,int side,InterpolationType interptype
			) const;

  bool WasSlopeLimited(int index)const;

private:

  SpatialReconstruction& sr_;
  EquationOfState const& eos_;
};

#endif // EOS_CONSISTENT_HPP
