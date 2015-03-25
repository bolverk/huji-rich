#include <cmath>
#include "collela.hpp"

using namespace std;

Collela::Collela(double ref, double a, double l, double offset):
  ref_(ref), a_(a), l_(l), offset_(offset) {}

double Collela::operator()(double x) const
{
  const double xc = x - offset_;
  if(abs(l_)<abs(xc))
    return ref_;
  else
    return ref_*(1+a_*pow((pow(xc/l_,2)-1),4));
}
