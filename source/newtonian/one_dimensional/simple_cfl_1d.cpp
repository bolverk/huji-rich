#include "simple_cfl_1d.hpp"

using std::max;

SimpleCFL1D::SimpleCFL1D(double cfl): cfl_(cfl)
{
  assert(cfl>0);
  assert(cfl<1);
}

double SimpleCFL1D::operator()
  (const SimulationState1D& ss,
   const EquationOfState& eos) const
{
  double idt = 0;
  for(size_t i=0;i<ss.getCells().size();++i){
    const ComputationalCell& cell = ss.getCells().at(i);
    const double c = eos.dp2c(cell.density, cell.pressure);
    const double dx = ss.getVertices().at(i+1) - ss.getVertices().at(i);
    const double idt_candidate = (c+abs(cell.velocity))/dx;
    idt = max(idt, idt_candidate);
  }
  assert(idt>0);
  return cfl_/idt;
}
