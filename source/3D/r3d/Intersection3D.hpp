#ifndef INTERSECTION3D_HPP
#define INTERSECTION3D_HPP 1

#include "../GeometryCommon/Tessellation3D.hpp"


extern "C" 
{
#include "r3d.h"
}

void GetPlanes(vector<r3d_plane> &res, Tessellation3D const& tess, size_t index);

std::pair<bool, double> PolyhedraIntersection(Tessellation3D const& oldtess, Tessellation3D const& newtess, size_t newcell,
	size_t oldcell, r3d_poly &poly, vector<Vector3D> &vtemp, vector<size_t> &itemp, vector<size_t> &all_indeces,
	vector<vector<int> > &faceinds,vector<r3d_plane> *planes = 0);

#endif //INTERSECTION3D_HPP