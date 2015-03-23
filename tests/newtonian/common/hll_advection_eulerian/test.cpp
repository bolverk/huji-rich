#include <iostream>
#include <fstream>
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/common/hydrodynamics.hpp"

using namespace std;

namespace {
void write_primitive(Primitive const& p,
		     string const& fname)
{
  ofstream f(fname.c_str());
  for(int i=0;i<p.GetVarNo();++i){
    f << p[i] << endl;
  }
  f.close();
}

void write_flux(Conserved const& flux,
		string const& fname)
{
  ofstream f(fname.c_str());
  f << flux.Mass << endl;
  f << flux.Momentum.x << endl;
  f << flux.Momentum.y << endl;
  f << flux.Energy << endl;
  f.close();
}

void write_output(Primitive const& left,
		  Primitive const& right,
		  Conserved const& flux)
{
  write_primitive(left,"left.txt");
  write_primitive(right,"right.txt");
  write_flux(flux,"flux.txt");
}
}

int main(void)
{
  double left_density = 1;
  double right_density = 2;
  double pressure = 3;
  Vector2D velocity(4,0);
  IdealGas eos(5./3.);
  Primitive left = CalcPrimitive(left_density,pressure,velocity,eos);
  Primitive right = CalcPrimitive(right_density,pressure,velocity,eos);
  double edge_velocity = 0;
  Hllc rs;
  Conserved flux = rs(left,right,edge_velocity);

  write_output(left,right,flux);

  return 0;
}
