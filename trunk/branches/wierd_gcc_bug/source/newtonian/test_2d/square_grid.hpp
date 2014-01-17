#ifndef SQUARE_GRID_HPP
#define SQUARE_GRID_HPP 1

#include <vector>
#include "../../tessellation/geometry.hpp"

using namespace std;

vector<Vector2D> square_grid(double side, int np);

vector<Vector2D> offset_grid(vector<Vector2D> const& original,
			     Vector2D const& offset);

#endif // SQUARE_GRID_HPP
