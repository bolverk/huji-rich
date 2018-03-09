#include "utils.hpp"
#if defined(_MSC_VER)
/* Microsoft C/C++-compatible compiler */
#include <intrin.h>
#else
#include <x86intrin.h>
#endif
#include <cfloat>
#include <cmath>

bool is_nan(double x)
{
  int b1 = (x>=0);
  int b2 = (x<0);
  return (1 != b1+b2);
}

vector<double> linspace(double xl, double xh, int n)
{
  vector<double> res(n,0);
  for(size_t i=0;i<size_t(n);++i)
    res[i] = xl + (xh-xl)*static_cast<double>(i)/
      static_cast<double>(n-1);
  return res;
}

vector<double> arange(double x_min, double x_max, double dx)
{
  assert((x_max-x_min)/dx>0 && "dx has wrong sign");
  vector<double> res;
  for(double x=x_min;x<x_max;x+=dx)
    res.push_back(x);
  return res;
}

double min(vector<double> const& v)
{
  double res = v[0];
  for(size_t i=1;i<v.size();++i)
    res = std::min(res,v[size_t(i)]);
  return res;
}

double max(vector<double> const& v)
{
  double res = v[0];
  for(size_t i=1;i<v.size();++i)
    res = std::max(res,v[size_t(i)]);
  return res;
}

double fastsqrt(double x)
{
	if (x<FLT_MIN || x > FLT_MAX)
		return std::sqrt(x);
	double res = static_cast<double>(_mm_cvtss_f32(_mm_rsqrt_ss(_mm_set_ss(static_cast<float>(x)))));
	return x*res*(1.5 - 0.5*res*res*x);
}
