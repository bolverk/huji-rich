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
	binary_write_single_size_t(nfaces,file_handle);
	for(size_t i=0;i<nfaces;++i)
	{
		Face const& f=tess.GetFace(i);
		// write neighbors
		binary_write_single_size_t(f.neighbors.first,file_handle);
		binary_write_single_size_t(f.neighbors.second,file_handle);
		// write the vertices
		binary_write_single_size_t(f.vertices.size(),file_handle);
		for(size_t j=0;j<f.vertices.size();++j)
		{
			binary_write_single_double(f.vertices[j].x,file_handle);
			binary_write_single_double(f.vertices[j].y,file_handle);
			binary_write_single_double(f.vertices[j].z,file_handle);
		}
	}
	// write which faces each cell has
	for(size_t i=0;i<N;++i)
	{
		vector<size_t> temp=tess.GetCellFaces(i);
		binary_write_single_size_t(temp.size(),file_handle);
		for(size_t j=0;j<temp.size();++j)
			binary_write_single_size_t(temp[j],file_handle);
	}
	file_handle.close();
}
