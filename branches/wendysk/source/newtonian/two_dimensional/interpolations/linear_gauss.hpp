#ifndef LINEAR_GAUSS
#define LINEAR_GAUSS 1

#include "../spatial_reconstruction.hpp"
#include "../OuterBoundary.hpp"
#include "../HydroBoundaryConditions.hpp"
#include "../../../misc/universal_error.hpp"
#include <string>
#include <fstream>

/*! \brief Spatial reconstruction based on Gauss' theorem
  \details Details are in the Arepo paper
  \author Almog Yalinewich
 */class LinearGauss: public SpatialReconstruction
{
public:

  /*! \brief Class constructor
    \param obc The outer boundary conditions
    \param hbc The hydro boundary conditions
    \param slf Slope limiter flag
    \param soitf Second order in time flag
    \param delta_v The GradV*L/Cs ratio needed for slope limiter
    \param theta The theta from tess in slope limiter.
    \param delta_P The pressure ratio for shock detection
    \param rigidflag Flag whether not to use reflected cell in rigid boundary for slope construction and limit
   */
  LinearGauss(OuterBoundary const& obc,HydroBoundaryConditions const* hbc,
	  bool slf=true, bool soitf=false,double delta_v=0.2,double theta=0.5,
	  double delta_P=0.7,bool rigidflag=false);

  void Prepare(Tessellation const* tessellation,
	       vector<Primitive> const& cells,
	       double dt,vector<bool> const& mask,double time);

  Primitive Interpolate(Tessellation const* tess,
	  vector<Primitive> const& cells,double dt,Edge const& edge,
			int side,InterpolationType interptype) const;
  /*!
  \brief Returns the cell's slope
  \param cell_index The cell's index
  \return The slope
  */
  PrimitiveGradient2D GetSlope(int cell_index)const;

  bool WasSlopeLimited(int index)const;

  ~LinearGauss(void);

private:
  vector<bool> slope_limited_;
  vector<PrimitiveGradient2D> slopes_;
  OuterBoundary const* obc_;
  HydroBoundaryConditions const* hbc_;
  bool slf_;
  bool soitf_;
  double shockratio_,diffusecoeff_,pressure_ratio_;
  bool _rigidflag;
  
  LinearGauss(const LinearGauss& origin);
  LinearGauss& operator=(const LinearGauss& origin);
};

#endif // LINEAR_GAUSS
