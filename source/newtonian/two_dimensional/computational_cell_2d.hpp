/*! \file computational_cell.hpp
  \author Almog Yalinewich
  \brief Computational cell
*/

#ifndef COMPUTATIONAL_CELL_HPP
#define COMPUTATIONAL_CELL_HPP 1

/*#include <boost/container/small_vector.hpp>
  typedef boost::container::small_vector<double,0> tvector;
  typedef boost::container::small_vector<bool,0> svector;*/
#include <vector>
#include <array>
#include <string>
#define MAX_TRACERS 5
#define MAX_STICKERS 5

using std::vector;
using std::string;

typedef std::array<double, MAX_TRACERS> tvector;
typedef std::array<bool, MAX_STICKERS> svector;
#include <string>
#include "../../tessellation/geometry.hpp"
#ifdef RICH_MPI
#include "../../misc/serializable.hpp"
#endif // RICH_MPI

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

  static vector<string> tracerNames;
  static vector<string> stickerNames;

  //! \brief Tracers (can transfer from one cell to another)
  tvector tracers;

  //! \brief Stickers (stick to the same cell)
  svector stickers;

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
  size_t getChunkSize(void) const override;

  vector<double> serialize(void) const override;

  void unserialize
  (const vector<double>& data) override;
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
  //! \brief Slope in the x direction
  ComputationalCell xderivative;

  //! \brief Slope in the y direction
  ComputationalCell yderivative;

  /*!
    \brief Class constructor
    \param x The x derivative 
    \param y The y derivative
  */
  Slope(ComputationalCell const& x, ComputationalCell const& y);
  //! \brief Default constructor
  Slope(void);
#ifdef RICH_MPI
  size_t getChunkSize(void) const override;

  vector<double> serialize(void) const override;

  void unserialize(const vector<double>& data) override;
#endif//RICH_MPI
};

#endif // COMPUTATIONAL_CELL_HPP
