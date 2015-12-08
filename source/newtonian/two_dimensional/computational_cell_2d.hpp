/*! \file computational_cell.hpp
\author Almog Yalinewich
\brief Computational cell
*/

#ifndef COMPUTATIONAL_CELL_HPP
#define COMPUTATIONAL_CELL_HPP 1

#include <boost/container/flat_map.hpp>
#include <string>
#include "../../tessellation/geometry.hpp"
#ifdef RICH_MPI
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include "flat_map_serial.hpp"
#endif // RICH_MPI

using std::string;

//! \brief Computational cell
class ComputationalCell
#ifdef RICH_MPI
  : public Serializable
#endif // RICH_MPI
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

  /*! \brief Self increment operator
    \param other Addition
    \return Reference to self
   */
	ComputationalCell& operator+=(ComputationalCell const& other);
	/*! \brief Self reduction operator
	\param other Reduction
	\return Reference to self
	*/
	ComputationalCell& operator-=(ComputationalCell const& other);

	/*! \brief Self multiplication operator
	\param s The scalar to multiply
	\return Reference to self
	*/
	ComputationalCell& operator*=(double s);

  /*! \brief Self decrement operator
    \param other difference
    \return Reference to self
   */
  	ComputationalCell& operator=(ComputationalCell const& other);

#ifdef RICH_MPI
  /*! \brief Serializer
    \param ar Archive
    \param int version
   */
	template<class Archive>
	void serialize
		(Archive& ar,
		 const unsigned int /*version*/)
	{
		ar & density;
		ar & pressure;
		ar & velocity;
		ar & tracers;
		ar & stickers;
	}

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

void ComputationalCellAddMult(ComputationalCell &res, ComputationalCell const& other, double scalar);

void ReplaceComputationalCell(ComputationalCell &cell, ComputationalCell const& other);

//! \brief Class for spatial interpolations
class Slope
#ifdef RICH_MPI
	: public Serializable
#endif // RICH_MPI
{
public:
	ComputationalCell xderivative;
	ComputationalCell yderivative;

	/*!
	\brief Class constructor
	\param x The x derivative 
	\param y The y derivative
	*/
	Slope(ComputationalCell const& x, ComputationalCell const& y);
	//! \brief Default constructor
	Slope(void);

	size_t getChunkSize(void) const;

	vector<double> serialize(void) const;

	void unserialize(const vector<double>& data);
};
#endif // COMPUTATIONAL_CELL_HPP
