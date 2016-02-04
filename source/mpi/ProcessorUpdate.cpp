#include "ProcessorUpdate.hpp"
#ifdef RICH_MPI
#include <mpi.h>
#endif

ProcessorUpdate::~ProcessorUpdate(void){}

#ifdef RICH_MPI
double ProcessorUpdate::GetLoadImbalance(Tessellation const& tlocal,int &total)const
{
	int ws;
	MPI_Comm_size(MPI_COMM_WORLD, &ws);
	vector<int> N(ws,0);
	int n=tlocal.GetPointNo();
	MPI_Allgather(&n, 1, MPI_INT, &N[0], 1, MPI_INT, MPI_COMM_WORLD);
	total=0;
	for(size_t i=0;i<N.size();++i)
		total+=N[i];
	return 1.0*(*std::max_element(N.begin(),N.end())*ws)/(1.0*total);
}
#endif
