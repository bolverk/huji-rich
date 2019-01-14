#ifndef INTERSECTION3D_HPP
#define INTERSECTION3D_HPP 1

#include <array>

#include "../GeometryCommon/Tessellation3D.hpp"


extern "C" 
{
#include "r3d.h"
}

void GetPlanes(vector<r3d_plane> &res, Tessellation3D const& tess, size_t index);

bool GetPoly(Tessellation3D const & oldtess, size_t oldcell, r3d_poly &poly, point_vec &itemp, vector<size_t> &all_indeces,
	vector<vector<int> > &faceinds);

std::pair<bool, std::array<double,4> > PolyhedraIntersection(Tessellation3D const& newtess, size_t newcell,r3d_poly &poly,
	vector<r3d_plane> *planes = 0);

#endif //INTERSECTION3D_HPP