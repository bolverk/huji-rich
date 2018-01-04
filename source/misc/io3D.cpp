#include "io3D.hpp"

void write_vec3d(std::vector<Vector3D> const & vec, std::string const & fname)
{
	std::ofstream file_handle(fname.c_str(), std::ios::out | std::ios::binary);
	assert(file_handle.is_open());
	size_t stemp = vec.size();
	binary_write_single_int(static_cast<int>(stemp), file_handle);
	for (std::size_t i = 0; i < stemp; ++i)
	{
		binary_write_single_double(vec[i].x, file_handle);
		binary_write_single_double(vec[i].y, file_handle);
		binary_write_single_double(vec[i].z, file_handle);
	}
	file_handle.close();
}

std::vector<Vector3D> read_vec3d(std::string fname)
{
	vector<Vector3D> res;
	std::ifstream fh(fname.c_str(), std::ios::binary);
	int npoints;
	fh.read(reinterpret_cast<char*>(&npoints), sizeof(int));
	for (int i = 0; i < npoints; ++i)
	{
		double x = 0, y = 0, z = 0;
		fh.read(reinterpret_cast<char*>(&x), sizeof(double));
		fh.read(reinterpret_cast<char*>(&y), sizeof(double));
		fh.read(reinterpret_cast<char*>(&z), sizeof(double));
		res.push_back(Vector3D(x, y, z));
	}
	fh.close();
	return res;
}
