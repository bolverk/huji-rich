#include "ProcessorUpdate.hpp"
#ifdef RICH_MPI
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#endif

ProcessorUpdate::~ProcessorUpdate(void){}

#ifdef RICH_MPI
double ProcessorUpdate::GetLoadImbalance(Tessellation const& tlocal)const
{
	const boost::mpi::communicator world;
	int ws=world.size();
	vector<int> N(ws,0);
	int n=tlocal.GetPointNo();
	boost::mpi::all_gather(world, n, N);
	int total=0;
	for(size_t i=0;i<N.size();++i)
		total+=N[i];
	return 1.0*(*std::max_element(N.begin(),N.end())*ws)/(1.0*total);
}
#endif
