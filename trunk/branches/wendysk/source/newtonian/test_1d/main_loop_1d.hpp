#include <string>
#include "../one_dimensional/hdsim.hpp"

using namespace std;

namespace simulation1d{

  class TerminationCondition
  {
  public:

    virtual bool should_continue(hdsim1D const& sim) = 0;

    virtual ~TerminationCondition(void) = 0;

  private:
  };

  class SafeTimeTermination: public TerminationCondition
  {
  public:

    SafeTimeTermination(double termination_time,
			int max_cycles);

    bool should_continue(hdsim1D const&  sim);

  private:
    
    const double termination_time_;
    const int max_cycles_;
  };

  class DiagnosticsFunction
  {
  public:

    virtual void diagnose(hdsim1D const& sim) = 0;
    
    virtual ~DiagnosticsFunction(void) = 0;
  };

  class WriteTime
  {
  public:

    WriteTime(string const& fname);

    void diagnose(hdsim1D const& sim);

  private:
    string fname_;
  };

  void main_loop(hdsim1D& sim,
		 TerminationCondition& term_cond,
		 int time_order,
		 DiagnosticsFunction* diag = 0);

  void main_loop(hdsim1D& sim,
		 double final_time,
		 int max_iter=1e6,
		 int time_order=1,
		 string const& time_log="");
}
