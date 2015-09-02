/*! \file computational_cell.hpp
\author Almog Yalinewich
\brief Computational cell
*/

#ifndef COMPUTATIONAL_CELL_HPP
#define COMPUTATIONAL_CELL_HPP 1

#include <map>
#include <boost/container/flat_map.hpp>
#include <string>
#include "../../tessellation/geometry.hpp"

//! \brief Computational cell
class ComputationalCell
{
public:

	//! \brief Density
	double density;

	//! \brief Pressure
	double pressure;

	//! \brief Velocity
	Vector2D velocity;

	//! \brief Tracers (can transfer from one cell to another)
	// std::map<std::string,double> tracers;
	boost::container::flat_map<std::string, double> tracers;

	//! \brief Stickers (stick to the same cell)
	//std::map<std::string,bool> stickers;
	boost::container::flat_map<std::string, bool> stickers;

	/*!
	\brief Copy constructor
	\param other The cell to copy
	*/
	ComputationalCell(ComputationalCell const& other);
	/*!
	\brief Default constructor
	*/
	ComputationalCell(void);

	ComputationalCell& operator+=(ComputationalCell const& other);
};

/*! \brief Term by term addition
\param p1 Computational Cell
\param p2 Computational Cell
\return Computational Cell
*/
ComputationalCell operator+(ComputationalCell const& p1, ComputationalCell const& p2);

/*! \brief Term by term subtraction
\param p1 Computational Cell
\param p2 Computational Cell
\return Computational Cell
*/
ComputationalCell operator-(ComputationalCell const& p1, ComputationalCell const& p2);

/*! \brief Scalar division
\param p Computational Cell
\param s Scalar
\return Computational Cell
*/
ComputationalCell operator/(ComputationalCell const& p, double s);

/*! \brief Scalar multiplication on the right
\param p Computational Cell
\param s Scalar
\return Computational Cell
*/
ComputationalCell operator*(ComputationalCell const& p, double s);

/*! \brief Scalar multiplication on the left
\param s Scalar
\param p Computational Cell
\return Computational Cell
*/
ComputationalCell operator*(double s, ComputationalCell const& p);
#endif // COMPUTATIONAL_CELL_HPP
