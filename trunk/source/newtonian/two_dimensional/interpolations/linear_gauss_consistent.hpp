/*! \file linear_gauss_consistent.hpp
  \brief Linear interpolation that guarantees compliance with the equation of state
  \author Almog Yalinewich
*/

#ifndef LINEAR_GAUSS_CONSISTENT
#define LINEAR_GAUSS_CONSISTENT 1

#include <string>
#include <fstream>
#include <algorithm>
#include "../spatial_reconstruction.hpp"
#include "../OuterBoundary.hpp"
#include "../HydroBoundaryConditions.hpp"
#include "../../../misc/universal_error.hpp"
#include "../../common/equation_of_state.hpp"
#include "../../../misc/utils.hpp"
#include <boost/array.hpp>
#include <boost/container/static_vector.hpp>

//! \brief Gradient of the hydrodynamic variables, without energy and sound speed
class ReducedPrimitiveGradient2D
{
public:

  //! \brief Null constructor (sets everything to zero)
  ReducedPrimitiveGradient2D(void);

  /*! \brief Class constructor
    \param d Density gradient
    \param p Pressure gradient
    \param vx x velocity gradient
    \param vy y Velocity gradient
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

  //! \brief Gradient of the x velocity
  Vector2D xvelocity;

  //! \brief Gradient of the y velocity
  Vector2D yvelocity;

  //! \brief Gradient of the tracers
  vector<Vector2D> tracers;

  /*! \brief Term by term sum
    \param source Gradient to add
    \return sum
   */
  ReducedPrimitiveGradient2D& operator+=
  (ReducedPrimitiveGradient2D const& source);

  /*! \brief Scalar multiplication
    \param s Scalar
    \return Gradient multiplied by a scalar
   */
  ReducedPrimitiveGradient2D& operator*=
  (double s);
};

/*! \brief Spatial reconstruction based on Gauss' theorem
  \details Details are in the Arepo paper
  \author Almog Yalinewich
*/class LinearGaussConsistent: public SpatialReconstruction
  {
  public:

    /*! \brief Class constructor
      \param eos Equation of state
      \param obc The outer boundary conditions
      \param hbc The hydro boundary conditions
      \param slf Slope limiter flag
      \param soitf Second order in time flag
      \param delta_v The GradV*L/Cs ratio needed for slope limiter
      \param theta The theta from tess in slope limiter.
      \param delta_P The pressure ratio for shock detection
      \param rigidflag Flag whether not to use reflected cell in rigid boundary for slope construction and limit
    */
    LinearGaussConsistent
    (EquationOfState const& eos,
     OuterBoundary const& obc,
     HydroBoundaryConditions const* hbc,
     bool slf=true, bool soitf=false,double delta_v=0.2,double theta=0.5,
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

    ~LinearGaussConsistent(void);

  private:
    EquationOfState const& eos_;
    vector<ReducedPrimitiveGradient2D> rslopes_;
    OuterBoundary const* obc_;
    HydroBoundaryConditions const* hbc_;
    bool slf_;
    bool soitf_;
    double shockratio_,diffusecoeff_,pressure_ratio_;
    bool _rigidflag;

    LinearGaussConsistent(const LinearGaussConsistent& origin);
    LinearGaussConsistent& operator=(const LinearGaussConsistent& origin);
};

#endif // LINEAR_GAUSS_CONSISTENT
