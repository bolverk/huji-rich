/*! \file RiemannSolver3D.hpp
\brief Base class for Riemann solver in 3D
\author Elad Steinberg
*/

#ifndef RIEMANN_SOLVER3D_HPP
#define RIEMANN_SOLVER3D_HPP 1

#include "computational_cell.hpp"
#include "conserved_3d.hpp"

//! \brief Base class for Riemann solver
class RiemannSolver3D
{
public:

  RiemannSolver3D(void);

	/*! \brief Solve Riemann porblme
	\param left Primitive variables on the left side
	\param right Primitive variables on the right side
	\param velocity Velocity of the face
	\param eos The equation of state
	\param tsn The names of the tracers and stickers
	\param normaldir Normal to the face
	\return Corrected fluxes
	*/
	virtual Conserved3D operator()
	(ComputationalCell3D const& left,
	 ComputationalCell3D const& right,
	 double velocity,
	 EquationOfState const& eos,
	 Vector3D const& normaldir) const = 0;

  /*
  RiemannSolver3D(void);
  RiemannSolver3D(const RiemannSolver3D&) = default;
  RiemannSolver3D(RiemannSolver3D&&) = default;
  RiemannSolver3D& operator=(const RiemannSolver3D&) = default;
  RiemannSolver3D& operator=(RiemannSolver3D&&) = default;
  */

  virtual ~RiemannSolver3D(void);
protected:

  /*! \brief Copy constructor
    \param rs Riemann solver
   */
  RiemannSolver3D(const RiemannSolver3D& rs);
};

#endif // RIEMANN_SOVLER3D_HPP
