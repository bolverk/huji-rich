/*! \file extensive.hpp
  \author Almog Yalinewich
  \brief Extensive variables
 */

#ifndef EXTENSIVE_HPP
#define EXTENSIVE_HPP 1

#include "boost/container/flat_map.hpp"
#include "../../tessellation/geometry.hpp"

//! \brief Extensive variables
class Extensive
{
public:
  
  //! \brief mass
  double mass;

  //! \brief energy
  double energy;

  //! \brief momentum
  Vector2D momentum;

  //! \brief tracers
  boost::container::flat_map<std::string,double> tracers;

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

  //! \brief Null constructor
  Extensive(void);

  /*!
  \brief constructor for extensive with a tracer list. All tracers start with zero.
  \param Tracers The tracers 
  */
  Extensive(boost::container::flat_map<std::string, double> const& Tracers);
  
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

#endif // EXTENSIVE_HPP
