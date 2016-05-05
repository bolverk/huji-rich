#include "index2filename.hpp"
#include "../../misc/int2str.hpp"
#ifdef RICH_MPI
#include <mpi.h>
#endif

Index2FileName::~Index2FileName(void) {}

Rubric::Rubric(const string& prefix,
	       const string& postfix):
  prefix_(prefix), postfix_(postfix) {}

string Rubric::operator()(int index)
{
#ifndef RICH_MPI
  return prefix_+int2str(index)+postfix_;
#else
	int rank=0;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	return prefix_+int2str(index)+"_"+int2str(rank)+postfix_;
#endif
}
