#include "BinaryLogger.hpp"

using namespace std;

BinaryLogger::BinaryLogger(string const& filename):filename_(filename){}

BinaryLogger::~BinaryLogger(void){}

void BinaryLogger::Log(Tessellation3D const& tess)const
{
	ofstream file_handle(filename_.c_str(),ios::binary);
	if(!file_handle.good())
		throw UniversalError("Error opening BinaryLogger file!!");

	size_t N=tess.GetPointNo();
	binary_write_single_size_t(N,file_handle);

	// Write the points
	for(size_t i=0;i<N;++i)
	{
		Vector3D temp=tess.GetMeshPoint(i);
		binary_write_single_double(temp.x,file_handle);
		binary_write_single_double(temp.y,file_handle);
		binary_write_single_double(temp.z,file_handle);
	}
	
	// Write the faces
	size_t nfaces=tess.GetTotalFacesNumber();
	for(size_t i=0;i<nfaces;++i)
	{
		Face const& f=tess.GetFace(i);
		// write neighbors
		binary_write_single_size_t(f.neighbors.first,file_handle);
		binary_write_single_size_t(f.neighbors.second,file_handle);
	}

}