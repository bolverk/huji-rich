#ifndef RANDOM_PERT
#define RANDOM_PERT 1

#include <vector>
#include "../../tessellation/geometry.hpp"

using namespace std;

vector<Vector2D> add_random_pert(vector<Vector2D> const& original,
				 double amp_x, double amp_y);

#endif // RANDOM_PERT
