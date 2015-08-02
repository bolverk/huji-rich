/*! \file extensive.hpp
  \author Almog Yalinewich
  \brief Extensive variables
 */

#ifndef EXTENSIVE_HPP
#define EXTENSIVE_HPP 1

#include <map>
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
  std::map<std::string,double> tracers;

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
};

/*! \brief Scalar product
  \param s Scalar
  \param e Source extensive variables
  \return Extensive multiplied by s
 */
Extensive operator*(const double s,
		    const Extensive& e);

Extensive operator+(const Extensive& e1,
		    const Extensive& e2);

Extensive operator-(const Extensive& e1,
		    const Extensive& e2);

#endif // EXTENSIVE_HPP
