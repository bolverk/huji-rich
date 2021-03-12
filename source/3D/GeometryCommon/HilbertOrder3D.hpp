/*! \file HilbertOrder3D.hpp
\brief Hilbert 3D-space filling curve
\author Itai Linial
*/

// Need t do, add option to define number of points the same for early break, add protection against two points the same

#ifndef HILBERTORDER3D_HPP
#define HILBERTORDER3D_HPP 1
#include <vector>
#include "Vector3D.hpp"
using std::vector;

/*!
\brief Returns the 3D-Hilbert curve ordering
\param cor The points
\return The Hilbert order indices
*/
vector<std::size_t> HilbertOrder3D(vector<Vector3D> const& cor);

/*! \brief Get the reorderd list of indices
  \param cor List of points
  \param ll Low left corner
  \param ur Upper right corner
  \param Hmax Upper index of real indices
  \return List of indices
 */
vector<std::size_t> GetGlobalHibertIndeces(vector<Vector3D> const& cor,Vector3D const& ll, Vector3D const& ur,size_t &Hmax);

#endif // HILBERTORDER_HPP
