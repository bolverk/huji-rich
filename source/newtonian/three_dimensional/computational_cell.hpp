/*! \file computational_cell.hpp
  \brief A container for the hydrodynamic variables
  \author Almog Yalinewich
 */

#ifndef COMPUTATIONAL_CELL3D_HPP
#define COMPUTATIONAL_CELL3D_HPP 1

#include <array>
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

	//! \brief Internal energy
	double internal_energy;

	//! \brief Unique ID
	size_t ID;

	//! \brief Velocity
	Vector3D velocity;

	//! \brief Tracers
	std::array<double,MAX_TRACERS> tracers;

	//! \brief Stickers
	std::array<bool,MAX_STICKERS> stickers;

	//! \brief Null constructor
	ComputationalCell3D(void);

	/*! \brief Class constructor
	  \param density_i Density
	  \param pressure_i Pressure
	  \param velocity_i Velocity
	  \param internal_energy_i Internal energy per unit mass
	  \param ID_i The ID
	  */
	ComputationalCell3D(double density_i,
		double pressure_i,double internal_energy_i, size_t ID_i,
		const Vector3D& velocity_i);

	/*! \brief Class constructor
	  \param density_i Density
	  \param pressure_i Pressure
	  \param internal_energy_i Internal energy per unit mass
	  \param velocity_i Velocity
	  \param tracers_i Tracers
	  \param stickers_i Stickers
	  \param ID_i The ID
	  */
	ComputationalCell3D(double density_i,
		double pressure_i,
		double internal_energy_i,
		size_t ID_i,
		const Vector3D& velocity_i,
		const std::array<double, MAX_TRACERS>& tracers_i,
		const std::array<bool, MAX_STICKERS>& stickers_i);

  ComputationalCell3D(const ComputationalCell3D& other);

	/*! \brief Self increment operator
	\param other Addition
	\return Reference to self
	*/
	ComputationalCell3D& operator+=(ComputationalCell3D const& other);
	/*! \brief Self reduction operator
	\param other Reduction
	\return Reference to self
	*/
	ComputationalCell3D& operator-=(ComputationalCell3D const& other);

	/*! \brief Self multiplication operator
	\param s The scalar to multiply
	\return Reference to self
	*/
	ComputationalCell3D& operator*=(double s);

	/*! \brief Self decrement operator
	\param other difference
	\return Reference to self
	*/
	ComputationalCell3D& operator=(ComputationalCell3D const& other);

#ifdef RICH_MPI
	size_t getChunkSize(void) const;

	vector<double> serialize(void) const;

	void unserialize
		(const vector<double>& data);
#endif // RICH_MPI

};

/*! \brief Term by term addition
\param p1 Computational Cell
\param p2 Computational Cell
\return Computational Cell
*/
ComputationalCell3D operator+(ComputationalCell3D const& p1, ComputationalCell3D const& p2);

/*! \brief Term by term subtraction
\param p1 Computational Cell
\param p2 Computational Cell
\return Computational Cell
*/
ComputationalCell3D operator-(ComputationalCell3D const& p1, ComputationalCell3D const& p2);

/*! \brief Scalar division
\param p Computational Cell
\param s Scalar
\return Computational Cell
*/
ComputationalCell3D operator/(ComputationalCell3D const& p, double s);

/*! \brief Scalar multiplication on the right
\param p Computational Cell
\param s Scalar
\return Computational Cell
*/
ComputationalCell3D operator*(ComputationalCell3D const& p, double s);

/*! \brief Scalar multiplication on the left
\param s Scalar
\param p Computational Cell
\return Computational Cell
*/
ComputationalCell3D operator*(double s, ComputationalCell3D const& p);

void ComputationalCellAddMult(ComputationalCell3D &res, ComputationalCell3D const& other, double scalar);

void ReplaceComputationalCell(ComputationalCell3D &cell, ComputationalCell3D const& other);

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

	//! \brief Slope in the z direction
	ComputationalCell3D zderivative;
	/*!
	\brief Class constructor
	\param x The x derivative
	\param y The y derivative
	\param z The z derivative
	*/
	Slope3D(ComputationalCell3D const& x, ComputationalCell3D const& y, ComputationalCell3D const& z);
	//! \brief Default constructor
	Slope3D(void);
#ifdef RICH_MPI
	size_t getChunkSize(void) const;

	vector<double> serialize(void) const;

	void unserialize(const vector<double>& data);
#endif//RICH_MPI
};


#endif // COMPUTATIONAL_CELL3D_HPP
