#include "index2filename.hpp"
#include "../../misc/int2str.hpp"
#ifdef RICH_MPI
#include "../../mpi/mpi_macro.hpp"
#endif

Index2FileName::~Index2FileName(void) {}

Rubric::Rubric(const string& prefix,
	       const string& postfix):
  prefix_(prefix), postfix_(postfix) {}

string Rubric::operator()(int index)
{
  #ifdef RICH_MPI
  return prefix_+int2str(index)+"_"+int2str(get_mpi_rank())+postfix_;
  #else
  return prefix_+int2str(index)+postfix_;
  #endif
}
