#ifndef SPATIAL_RECONSTRUCTION_HPP
#define SPATIAL_RECONSTRUCTION_HPP 1

#include "../common/hydrodynamic_variables.hpp"
#include "../../tessellation/tessellation.hpp"

enum InterpolationType{InBulk,Boundary};

//! \brief Two dimensional gradient of primitive variables
class PrimitiveGradient2D
{
public:

  /*! \brief Class constructor
    \param xi x component
    \param yi y component
   */
  PrimitiveGradient2D(Primitive const& xi,
		      Primitive const& yi);

  /*! \brief Null constructor
    \details Sets everything to 0
   */
  PrimitiveGradient2D(void);

  /*! \brief Increment operation
    \param pg Incremental values
   */
  PrimitiveGradient2D& operator+=(PrimitiveGradient2D const& pg);

  //! \brief x component
  Primitive x;

  //! \brief y component
  Primitive y;
};

Vector2D CalcCentroid(Edge const& edge);

/*! \brief Spatial reconstruction of the primitive functions
  \author Almog Yalinewich
 */ 
class SpatialReconstruction
{
public:

	/*!
	\brief Calculates the slopes for all of the cells
	\param tessellation The tessellation
	\param cells The primitive cells
	\param dt The time step
	\param mask Whether tocalculate cell or not
	\param time The simulation time
	*/
  virtual void Prepare(Tessellation const* tessellation,
		       vector<Primitive> const& cells,
		       double dt,vector<bool> const& mask,double time) = 0;

  /*! \brief Interpolates the hydrodynamic variables near the edge
    \param tessellation Point and edge positions
    \param cells Hydrodynamic variables
	\param dt Time step
    \param edge The face on which the Riemann problem will be solved
    \param side The side of the edge
	\param interptype Type of interpolation
	\return The interpolated primitive
   */
  virtual Primitive Interpolate(Tessellation const* tessellation,
				vector<Primitive> const& cells,double dt,Edge const& edge,
				int side,InterpolationType interptype) const = 0;

  /*!
	\brief Returns whether a cell was slope limited aggresively or not
	\param index The cell index
	\return True if the the slope was limited only aggresivly
  */
  virtual bool WasSlopeLimited(int index)const=0;

  //! \brief Virtual destructor
  virtual ~SpatialReconstruction(void);
};

#endif // SPATIAL_RECONSTRUCTION_HPP
