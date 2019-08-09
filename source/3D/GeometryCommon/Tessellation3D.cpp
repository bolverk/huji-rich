#include "Tessellation3D.hpp"

Tessellation3D::~Tessellation3D(void)
{}

point_vec_v VectorValues(std::vector<Vector3D> const&v, point_vec const &index)
{
	if (index.empty() || v.empty())
		return point_vec_v();

	point_vec_v result;
	size_t N = index.size();
	for (std::size_t i = 0; i < N; ++i)
		result.push_back(v[static_cast<std::size_t>(index[i])]);
	return result;
}
