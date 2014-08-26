/*! \file LinearGaussArepo.hpp
  \brief Interpolation described in the Arepo paper
  \author Elad Steinberg
*/

#ifndef LINEAR_GAUSS_AREPO
#define LINEAR_GAUSS_AREPO 1

#include <string>
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
#include "../ReducedPrimitiveGradient2D.hpp"

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
    LinearGaussArepo(EquationOfState const& eos,OuterBoundary const& obc,
     HydroBoundaryConditions const& hbc,Acceleration& acc,
     bool slf=true,double delta_v=0.2,double theta=0.5,
     double delta_P=0.7,bool rigidflag=false);

    void Prepare(Tessellation const& tessellation,vector<Primitive> const& cells,
		 vector<vector<double> > const& tracers,vector<bool> const& isrelevant,
		 double dt,double time);

    Primitive Interpolate(Tessellation const& tess,
			  vector<Primitive> const& cells,double dt,Edge const& edge,
			  int side,InterpolationType interptype,Vector2D const& vface) const;

    vector<double> interpolateTracers
    (Tessellation const& tess,vector<Primitive> const& cells,
     vector<vector<double> > const& tracers,double dt,Edge const& edge,
     int side,InterpolationType interptype,Vector2D const& vface) const;

	vector<ReducedPrimitiveGradient2D>& GetGradients(void);

    ~LinearGaussArepo(void);

  private:
    EquationOfState const& eos_;
    vector<ReducedPrimitiveGradient2D> rslopes_;
    OuterBoundary const& obc_;
    HydroBoundaryConditions const& hbc_;
    Acceleration& acc_;
    bool slf_;
    double shockratio_,diffusecoeff_,pressure_ratio_;
    bool _rigidflag;
    double time_;

    LinearGaussArepo(const LinearGaussArepo& origin);
    LinearGaussArepo& operator=(const LinearGaussArepo& origin);
};

#endif // LINEAR_GAUSS_AREPO