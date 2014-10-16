#include "ProcessorUpdate.hpp"

ProcessorUpdate::~ProcessorUpdate(void){}

double ProcessorUpdate::GetLoadImbalance(Tessellation const& tlocal)const
{
	int ws=get_mpi_size();
	vector<int> N0(ws,0),N(ws,0);
	int n=tlocal.GetPointNo();
	N0[get_mpi_rank()]=n;
	MPI_Allgather(&n,1,MPI_INT,&N[0],1,MPI_INT,MPI_COMM_WORLD);
	int total=0;
	for(size_t i=0;i<N.size();++i)
		total+=N[i];
	return 1.0*(*std::max_element(N.begin(),N.end())*ws)/(1.0*total);
}
