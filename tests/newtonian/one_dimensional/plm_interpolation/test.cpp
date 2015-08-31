#include <iostream>
#include <cmath>
#include "source/newtonian/one_dimensional/spatial_distribution1d.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/one_dimensional/hdsim.hpp"
#include "source/misc/utils.hpp"
#include "source/newtonian/one_dimensional/pcm1d.hpp"
#include "source/newtonian/one_dimensional/plm1d.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/test_1d/sine_wave.hpp"

using namespace std;

namespace {
  /*! \brief Spatial profile parabolic in the coordinate
   */
  class Parabola: public SpatialDistribution1D
  {
  public:

    Parabola(double a, double b, double c):
      a_(a), b_(b), c_(c) {}

    double operator()(double x) const
    {
      return a_*pow(x,2)+b_*x+c_;
    }

  private:

    double a_;
    double b_;
    double c_;
  };
}

namespace {
  Primitive CalcPrimitive(double density, double pressure,
			  Vector2D const& velocity, 
			  EquationOfState const& eos)
  {
    Primitive res;
    res.Density = density;
    res.Pressure = pressure;
    res.Velocity = velocity;
    res.Energy = eos.dp2e(res.Density, res.Pressure);
    res.SoundSpeed = eos.dp2c(res.Density, res.Pressure);
    return res;
  }

  vector<Primitive> InitialiseCells
  (vector<double> const& vertices,
   SpatialDistribution1D const& density,
   SpatialDistribution1D const& pressure,
   SpatialDistribution1D const& paravelocity,
   SpatialDistribution1D const& perpvelocity,
   EquationOfState const& eos)
  {
    vector<Primitive> res(vertices.size()-1);
    for(size_t i = 0; i<vertices.size() - 1; i++){
      double r = 0.5*(vertices[i] + vertices[i+1]);
      double d = density(r);
      double p = pressure(r);
      Vector2D v(paravelocity(r),
		 perpvelocity(r));
      res[i] = CalcPrimitive(d, p, v, eos);
    }
    return res;
  }
}

/*! \brief Contains all data needed for a test
 */
class TestData
{
public:
  
  TestData(int np):
    edges_(linspace(0,1,np)),
    density_(1),
    pressure_(1),
    xvelocity1_(1,1,0,0),
    xvelocity2_(0.01,1,0),
    xvelocity_(xvelocity1_),
    yvelocity_(0),
    eos_(5./3.),
    cells_(InitialiseCells(edges_,
			   density_,
			   pressure_,
			   xvelocity_,
			   yvelocity_,
			   eos_)) {}

  vector<Primitive> const& getCells(void) const
  {
    return cells_;
  }

  vector<double> const& getEdges(void) const
  {
    return edges_;
  }

  SpatialDistribution1D const& getDistribution
  (string const& vname) const
  {
    if("density"==vname)
      return density_;
    else if("pressure"==vname)
      return pressure_;
    else if("xvelocity"==vname)
      return xvelocity_;
    else if("yvelocity"==vname)
      return yvelocity_;
    else
      throw "Unknown variable name "+vname;
  }

private:
  vector<double> edges_;
  Uniform density_;
  Uniform pressure_;
  SineWave xvelocity1_;
  Parabola xvelocity2_;
  SpatialDistribution1D const& xvelocity_;
  Uniform yvelocity_;
  IdealGas eos_;
  vector<Primitive> cells_;
};

namespace {
  double test_interp(SpatialReconstruction1D const& sr,
		     TestData const& test_data)
  {
    double res = 0;
    int counter = 0;
    //  for(int i=0;i<(int)test_data.getEdges().size();++i){
    const int margin = 5;
    for(size_t i=margin;i<test_data.getEdges().size()-margin;++i){
      double v_a = test_data.getDistribution("xvelocity")(test_data.getEdges()[i]);
      if(i>0){
	double v_l = sr.InterpState
	  (test_data.getEdges(),
	   test_data.getCells(),
	   0,
	   i,0,0).Velocity.x;
	res += pow(v_a-v_l,2);
	++counter;
      }
      if(i<test_data.getEdges().size()-1){
	double v_r = sr.InterpState
	  (test_data.getEdges(),
	   test_data.getCells(),
	   0,
	   i,1,0).Velocity.x;
	res += pow(v_a-v_r,2);
	++counter;
      }
    }
    return sqrt(res/counter);
  }
}

namespace{
  void write_interp_test(SpatialReconstruction1D const& sr,
			 TestData const& test_data,
			 string const& fname)
  {
    ofstream f(fname.c_str());
  
    for(size_t i=0;i<test_data.getEdges().size();++i){
      if(i>0)
	f << test_data.getEdges()[i] << " "
	  << test_data.getDistribution("xvelocity")(test_data.getEdges()[i]) << " "
	  << sr.InterpState
	  (test_data.getEdges(), test_data.getCells(),0,i,0,0).Velocity.x << endl;
      if(i<test_data.getEdges().size()-1)
	f << test_data.getEdges()[i] << " "
	  << test_data.getDistribution("xvelocity")(test_data.getEdges()[i]) << " "
	  << sr.InterpState
	  (test_data.getEdges(), test_data.getCells(),0,i,1,0).Velocity.x << endl;
    }
    f.close();
  }
}

namespace {
  vector<int> arange(int v_init,
		     int v_stop,
		     int v_step)
  {
    vector<int> res((v_stop-v_init)/v_step,0);
    for(size_t i=0;i<res.size();++i){
      res[i] = v_init+v_step*static_cast<int>(i);
    }
    return res;
  }
}

namespace {
  vector<double> convergence_curve
  (vector<int> const& np_list,
   SpatialReconstruction1D const& sr)
  {
    vector<double> res(np_list.size(),0);
    for(size_t i=0;i<np_list.size();++i)
      res[i] = test_interp(sr,
			   TestData(np_list[i]));
    return res;
  }
}

int main(void)
{
  write_interp_test(PCM1D(),
		    TestData(20),
		    "pcm_20_profile.txt");
  write_interp_test(PLM1D(),
		    TestData(20),
		    "plm_20_profile.txt");

  const vector<int> np_list = arange(20,10000,100);
  const vector<double> cc_pcm = convergence_curve
    (np_list,PCM1D());
  const vector<double> cc_plm = convergence_curve
    (np_list,PLM1D());

  write_vector(np_list,"np_list.txt");
  write_vector(cc_pcm,"cc_pcm.txt");
  write_vector(cc_plm,"cc_plm.txt");
  
  return 0;
}
