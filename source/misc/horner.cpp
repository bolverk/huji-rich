#include "horner.hpp"

using std::size_t;

double horner(const vector<double> coef_list, double x)
{
  double res = 0;
  for(size_t i=0;i<coef_list.size();++i){
    res *= x;
    res += coef_list[i];
  }
  return res;
}
