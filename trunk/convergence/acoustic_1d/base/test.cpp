#include <iostream>
#include <fstream>
#include <cmath>
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/one_dimensional/hdsim.hpp"
#include "source/newtonian/one_dimensional/pcm1d.hpp"
#include "source/newtonian/one_dimensional/plm1d.hpp"
#include "source/newtonian/one_dimensional/eos_consistent.hpp"
#include "source/newtonian/one_dimensional/spatial_distribution1d.hpp"
#include "source/newtonian/one_dimensional/eulerian1d.hpp"
#include "source/newtonian/one_dimensional/lagrangian1d.hpp"
#include "source/newtonian/one_dimensional/periodic_1d.hpp"
#include "source/newtonian/one_dimensional/diagnostics_1d.hpp"
#include "source/misc/int2str.hpp"
#include "source/misc/utils.hpp"
#include "source/convergence/external_data.hpp"
#include "source/convergence/choose_between.hpp"
#include "source/misc/simple_io.hpp"

using namespace std;

class SineWave: public SpatialDistribution1D
{
public:

  /*! \brief Class constructor
    \param amp Amplitude
    \param k Wavenumber
    \param ph Phase
    \param offset Offset
   */
  SineWave(double amp, double k, double ph, double offset):
    amp_(amp), k_(k), ph_(ph), offset_(offset) {}

  double EvalAt(double x) const
  {
    return amp_*sin(k_*x+ph_)+offset_;
  }

private:

  const double amp_;
  const double k_;
  const double ph_;
  const double offset_;
};

/*! \brief Initial conditions for simple waves
 */
class SimpleWaveIC
{
public:

  /*! \brief Class constructor
    \param d0 Mean density
    \param p0 Mean pressure 
    \param eos Equation state
    \param dd Density perturbation amplitude
    \param wavelength Wavelength
  */
  SimpleWaveIC(double d0, double p0,
	       EquationOfState const& eos,
	       double dd,
	       double wavelength):
    c0_(eos.dp2c(d0,p0)),
    k_(2*M_PI/wavelength),
    density_(dd,k_,0,d0),
    pressure_(dd*pow(c0_,2),
	      k_,0,p0),
    xvelocity_(c0_*dd/d0,
	       k_,0,0),
    yvelocity_(0) {}

  /*! \brief Returns spatial profile
    \param pname Name of hydrodynamic variable
   */
  SpatialDistribution1D const& getProfile(string const& pname)
  {
    if("density"==pname)
      return density_;
    else if("pressure"==pname)
      return pressure_;
    else if("xvelocity"==pname)
      return xvelocity_;
    else if("yvelocity"==pname)
      return yvelocity_;
    else
      throw "Unknown profile name "+pname;
  }

private:
  const double c0_;
  const double k_;
  const SineWave density_;
  const SineWave pressure_;
  const SineWave xvelocity_;
  const Uniform yvelocity_;
};

class SimData
{
public:

  /*! \brief Class constructor
    \param exd Data from external file
   */
  SimData(ExternalData const& exd):
    time_int_ord_(exd.getTimeIntegrationOrder()),
    vertices_(linspace(0,1,exd.getPointNumber())),
    eos_(5./3.),
    plm_(plm_naive_,eos_),
    init_cond_(1,3./5.,eos_,1e-6,1),
    lagrangian_(true),
    sim_(vertices_,
	 choose_between<SpatialReconstruction1D>
	 (exd.getInterpolationFlag(),
	  plm_,pcm_),
	 init_cond_.getProfile("density"),
	 init_cond_.getProfile("pressure"),
	 init_cond_.getProfile("xvelocity"),
	 init_cond_.getProfile("yvelocity"),
	 eos_,
	 rs_,
	 choose_between<VertexMotion>
	 (exd.getVertexMotionFlag(),
	  lagrangian_,eulerian_),
	 bc_)
  {}

  /*! \brief Returns time integration order
   */
  int getTimeIntegrationOrder(void)
  {
    return time_int_ord_;
  }

  /*! \brief Returns a reference to the simulation
   */
  hdsim1D& getSim(void)
  {
    return sim_;
  }

private:
  int time_int_ord_;
  vector<double> vertices_;
  IdealGas eos_;
  PCM1D pcm_;
  PLM1D plm_naive_;
  EOSConsistent plm_;
  SimpleWaveIC init_cond_;
  Hllc rs_;
  Eulerian1D eulerian_;
  Lagrangian1D lagrangian_;
  Periodic1D bc_;
  hdsim1D sim_;
};

vector<string> get_output_hydro_vars(void)
{
  vector<string> res;
  res.push_back("center");
  res.push_back("density");
  res.push_back("pressure");
  res.push_back("xvelocity");
  return res;
}

void write_snapshot(hdsim1D const& sim, int cycle)
{
  vector<string> output_vars = get_output_hydro_vars();
  for(int i=0;i<(int)output_vars.size();++i)
    write_cells_property(sim,output_vars[i],
			 output_vars[i]+"_list."+
			 int2str(cycle)+".txt",10); 
}

void main_loop(hdsim1D& sim, int tio)
{
  const double tf = 1;
  const int max_cycle = (int)1e6;
  Timer timer(0,tf/10);

  write_snapshot(sim,0);

  while(sim.GetTime()<tf){
    
    if(tio==0)
      sim.TimeAdvance();
    else
      sim.TimeAdvanceRK(tio);

    if(sim.GetCycle()>max_cycle){
      cout << "Max number of cycles exceeded" << endl;
      throw;
    }

    if(timer.isTime(sim.GetTime()))
      write_snapshot(sim,timer.getCycle());
  }
}

void write_output(hdsim1D const& sim)
{
  const int prec = 14;
  write_number(sim.GetTime(),"time.txt");
  write_cells_property(sim,"center","cell_centres.txt",prec);
  write_cells_property(sim,"density","densities.txt",prec);
  write_cells_property(sim,"pressure","pressures.txt",prec);
  write_cells_property(sim,"xvelocity","velocities.txt",prec);
}

int main(void)
{
  try{
  SimData sim_data(ExternalData("data_file.txt"));
  hdsim1D& sim = sim_data.getSim();
			
  main_loop(sim,sim_data.getTimeIntegrationOrder());

  write_output(sim);
  }
  catch(string const& eo){
    cout << eo << endl;
    throw;
  }
  catch(char const* eo){
    cout << eo << endl;
    throw;
  }
  
  return 0;
}
