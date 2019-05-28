/*! \file extensive.hpp
  \author Almog Yalinewich
  \brief Extensive variables
 */

#ifndef EXTENSIVE_HPP
#define EXTENSIVE_HPP 1

#include "../../tessellation/geometry.hpp"
#include "computational_cell_2d.hpp"
#ifdef RICH_MPI
#include "../../misc/serializable.hpp"
#endif // RICH_MPI

using std::string;

//! \brief Extensive variables
class Extensive
#ifdef RICH_MPI
  : public Serializable
#endif // RICH_MPI
{
public:
  //! \brief mass
  double mass;

  //! \brief energy
  double energy;

  //! \brief momentum
  Vector2D momentum;

  //! \brief tracers
  tvector tracers;

  Extensive(const Extensive& other);

  /*! \brief Assignment operator
    \param origin Original extensives variables
    \return Copy
   */
  Extensive& operator=(const Extensive& origin);

  /*! \brief Self difference operator
    \param diff Difference
    \return Reference to self
   */
  Extensive& operator-=(const Extensive& diff);

  /*! \brief Self addition operator
    \param diff Addition
    \return Reference to self
   */
  Extensive& operator+=(const Extensive& diff);

  /*! \brief Self multiplication by scalar
  \param scalar The scalar to multiply
  \return Reference to self
  */
  Extensive& operator*=(const double scalar);

  //! \brief Null constructor
  Extensive(void);

  /*!
  \brief constructor for extensive with a tracer list. All tracers start with zero.
  \param Tracers The tracers 
  */
  explicit Extensive(tvector const& Tracers);

#ifdef RICH_MPI
  size_t getChunkSize(void) const;

  vector<double> serialize(void) const;

  void unserialize
  (const vector<double>& data);

#endif

};

/*! \brief Scalar product
  \param s Scalar
  \param e Source extensive variables
  \return Extensive multiplied by s
 */
Extensive operator*(const double s,
		    const Extensive& e);

/*! \brief Addition
  \param e1 First argument
  \param e2 Second argument
  \return Sum of two extensives
 */
Extensive operator+(const Extensive& e1,
		    const Extensive& e2);

/*! \brief Difference
  \param e1 First argument
  \param e2 Second argument
  \return Difference of two extensives
 */
Extensive operator-(const Extensive& e1,
		    const Extensive& e2);
/*!
\brief Replaces the data in the extensive. The tracers should already be allocated
\param toreplace The extensive to change
\param other The data to copy
*/
void ReplaceExtensive(Extensive &toreplace, Extensive const& other);

#endif // EXTENSIVE_HPP
