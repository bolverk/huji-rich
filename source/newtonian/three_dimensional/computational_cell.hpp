/*! \file computational_cell.hpp
  \brief A container for the hydrodynamic variables
  \author Almog Yalinewich
 */

#ifndef COMPUTATIONAL_CELL3D_HPP
#define COMPUTATIONAL_CELL3D_HPP 1

#include "../../3D/GeometryCommon/Vector3D.hpp"
#include "../two_dimensional/computational_cell_2d.hpp"
#ifdef RICH_MPI
#include "../../misc/serializable.hpp"
#endif // RICH_MPI

 //! \brief Container for the hydrodynamic variables
class ComputationalCell3D
{
public:

	//! \brief Density
	double density;

	//! \brief Pressure
	double pressure;

	//! \brief Velocity
	Vector3D velocity;

	//! \brief Tracers
	vector<double> tracers;

	//! \brief Stickers
	vector<bool> stickers;

	//! \brief Null constructor
	ComputationalCell3D(void);

	/*! \brief Class constructor
	  \param density_i Density
	  \param pressure_i Pressure
	  \param velocity_i Velocity
	  */
	ComputationalCell3D(double density_i,
		double pressure_i,
		const Vector3D& velocity_i);

	/*! \brief Class constructor
	  \param density_i Density
	  \param pressure_i Pressure
	  \param velocity_i Velocity
	  \param tracers_i Tracers
	  */
	ComputationalCell3D(double density_i,
		double pressure_i,
		const Vector3D& velocity_i,
		const vector<double>& tracers_i,
		const vector<bool>& stickers_i);
};

//! \brief Class for 3D spatial interpolations
class Slope3D
#ifdef RICH_MPI
	: public Serializable
#endif // RICH_MPI
{
public:
	//! \brief Slope in the x direction
	ComputationalCell3D xderivative;

	//! \brief Slope in the y direction
	ComputationalCell3D yderivative;
	/*!
	\brief Class constructor
	\param x The x derivative
	\param y The y derivative
	*/
	Slope3D(ComputationalCell3D const& x, ComputationalCell3D const& y);
	//! \brief Default constructor
	Slope3D(void);
#ifdef RICH_MPI
	size_t getChunkSize(void) const;

	vector<double> serialize(void) const;

	void unserialize(const vector<double>& data);
#endif//RICH_MPI
};


#endif // COMPUTATIONAL_CELL3D_HPP
