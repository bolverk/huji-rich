#ifdef RICH_MPI
#include <mpi.h>
#endif
#include <boost/foreach.hpp>
#include "source/newtonian/test_2d/triangle_step.hpp"
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/newtonian/test_2d/main_loop_2d.hpp"

using namespace std;
using namespace simulation2d;

namespace {

  class WriteConserved: public DiagnosticFunction
  {
  public:

    WriteConserved(string const& fname):
      cons_(), fname_(fname) {}

    void operator()(hdsim const& sim)
    {
      cons_.push_back(total_conserved(sim));
    }

    ~WriteConserved(void)
    {
#ifdef RICH_MPI
      if(get_mpi_rank()==0){
#endif
	ofstream f(fname_.c_str());
	BOOST_FOREACH(Conserved c,cons_)
	  {
	    f << c.Mass << " "
	      << c.Momentum.x << " "
	      << c.Momentum.y << " "
	      << c.Energy << "\n";
	  }
	f.close();
#ifdef RICH_MPI
      }
#endif
    }

  private:
    vector<Conserved> cons_;
    const string fname_;
  };
}

int main(void)
{
#ifdef RICH_MPI
  MPI_Init(NULL, NULL);
#endif
  TriangleStep sim_data("lagrangian");
  hdsim& sim = sim_data.getSim();

  {
    SafeTimeTermination term_cond(0.05,1e6);
    WriteConserved diag("res.txt");
    main_loop(sim, 
	      term_cond, 
	      1,
	      &diag);
  }

#ifdef RICH_MPI
  MPI_Finalize();
#endif

  return 0;
}

