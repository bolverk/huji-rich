/*! \file random_pert.hpp
  \brief Adds a random perturbation to the position of all mesh generating points
  \author Almog Yalinewich
  \deprecated Replace with methods in mesh_generator.hpp file
 */

#ifndef RANDOM_PERT
#define RANDOM_PERT 1

#include <vector>
#include "../../tessellation/geometry.hpp"

using std::vector;

/*! \brief Adds a random perturbation to the position of mesh generating points
  \param original Unperturbed mesh generating points
  \param amp_x Amplitude in the x direction
  \param amp_y Amplitude in the y direction
  \return Perturbed positions of the mesh generating points
 */
vector<Vector2D> add_random_pert(vector<Vector2D> const& original,
				 double amp_x, double amp_y);

#endif // RANDOM_PERT
