#include <fstream>
#include "delaunay_logger.hpp"
#include "../misc/simple_io.hpp"

using namespace delaunay_loggers;

DelaunayLogger::DelaunayLogger() {}

void DelaunayLogger::output(vector<Vector2D> const& /*cor*/,
			    vector<facet> const& /*f*/) {}

DelaunayLogger::~DelaunayLogger(void) {}

BinaryLogger::BinaryLogger(string const& file_name):
  file_name_(file_name) {}

void BinaryLogger::output(vector<Vector2D> const& cor,
			  vector<facet> const& f)
{
  ofstream file_handle(file_name_.c_str(),std::ostream::binary);

  binary_write_single_int((int)cor.size(),file_handle);

  binary_write_single_int((int)f.size(),file_handle);

  for(int i=0;i<(int)cor.size();++i){
    binary_write_single_double(cor[i].x,file_handle);
    binary_write_single_double(cor[i].y,file_handle);
  }

  for(int i=0;i<(int)f.size();++i){
    for(int j=0;j<3;++j)
      binary_write_single_int(f[i].vertices[j],file_handle);
  }

  for(int i=0;i<(int)f.size();++i){
    for(int j=0;j<3;++j)
      binary_write_single_int(f[i].get_friend(j),file_handle);
  }

  file_handle.close();
}
