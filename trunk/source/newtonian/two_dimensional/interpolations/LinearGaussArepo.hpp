/*! \file LinearGaussArepo.hpp
  \brief Interpolation described in the Arepo paper
  \author Elad Steinberg
*/

#ifndef LINEAR_GAUSS_AREPO
#define LINEAR_GAUSS_AREPO 1

#include <string>
#include <fstream>
#include <algorithm>
#include "../spatial_reconstruction.hpp"
#include "../OuterBoundary.hpp"
#include "../HydroBoundaryConditions.hpp"
#include "../../../misc/universal_error.hpp"
#include "../../common/equation_of_state.hpp"
#include "../../../misc/utils.hpp"
#include "../source_terms/ConservativeForce.hpp"
#include <boost/array.hpp>
#include <boost/container/static_vector.hpp>

//! Gradients of the primitive variables without energy and sound speed
class ReducedPrimitiveGradient2D
{
public:

  //! \brief Null constructor (sets everything to zero)
  ReducedPrimitiveGradient2D(void);

  /*! \brief Class constructor
    \param d Density gradient
    \param p Pressure gradient
    \param vx X velocity gradient
    \param vy Y velocity gradient
    \param trace Tracers gradient
   */
  ReducedPrimitiveGradient2D
  (Vector2D const& d,
   Vector2D const& p,
   Vector2D const& vx,
   Vector2D const& vy,
   vector<Vector2D> const& trace);

  //! \brief Density gradient
  Vector2D density;

  //! \brief Pressure gradient
  Vector2D pressure;

  //! \brief x velocity gradient
  Vector2D xvelocity;

  //! \brief y velocity gradient
  Vector2D yvelocity;

  //! \brief Tracers gradient
  vector<Vector2D> tracers;

  /*! \brief Term by term addition
    \param source Gradient to add
    \return sum
   */
  ReducedPrimitiveGradient2D& operator+=
  (ReducedPrimitiveGradient2D const& source);

  /*! \brief Scalar multiplication
    \param s Scalar
    \return Reference to gradient multiplied by a scalar
   */
  ReducedPrimitiveGradient2D& operator*=
  (double s);
};

/*! \brief Spatial reconstruction based on Gauss' theorem
  \details Details are in the Arepo paper
  \author Almog Yalinewich
*/class LinearGaussArepo: public SpatialReconstruction
  {
  public:

    /*! \brief Class constructor
      \param eos Equation of state
      \param obc The outer boundary conditions
      \param hbc The hydro boundary conditions
      \param acc Acceleration
      \param slf Slope limiter flag
      \param delta_v The GradV*L/Cs ratio needed for slope limiter
      \param theta The theta from tess in slope limiter.
      \param delta_P The pressure ratio for shock detection
      \param rigidflag Flag whether not to use reflected cell in rigid boundary for slope construction and limit
    */
    LinearGaussArepo
    (EquationOfState const& eos,
     OuterBoundary const& obc,
     HydroBoundaryConditions const* hbc,Acceleration& acc,
     bool slf=true,double delta_v=0.2,double theta=0.5,
     double delta_P=0.7,bool rigidflag=false);

    void Prepare(Tessellation const* tessellation,vector<Primitive> const& cells,
		 vector<vector<double> > const& tracers,double dt,double time);

    Primitive Interpolate(Tessellation const* tess,
			  vector<Primitive> const& cells,double dt,Edge const& edge,
			  int side,InterpolationType interptype,Vector2D const& vface) const;

    vector<double> interpolateTracers
    (Tessellation const* tess,vector<Primitive> const& cells,
     vector<vector<double> > const& tracers,double dt,Edge const& edge,
     int side,InterpolationType interptype,Vector2D const& vface) const;

    ~LinearGaussArepo(void);

  private:
    EquationOfState const& eos_;
    vector<ReducedPrimitiveGradient2D> rslopes_;
    OuterBoundary const* obc_;
    HydroBoundaryConditions const* hbc_;
    Acceleration& acc_;
    bool slf_;
    double shockratio_,diffusecoeff_,pressure_ratio_;
    bool _rigidflag;
    double time_;

    LinearGaussArepo(const LinearGaussArepo& origin);
    LinearGaussArepo& operator=(const LinearGaussArepo& origin);
};

#endif // LINEAR_GAUSS_AREPO
