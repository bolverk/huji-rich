#include "utils.hpp"

bool is_nan(double x)
{
  int b1 = (x>0);
  int b2 = (x<0);
  int b3 = (x==0);
  return (1 != b1+b2+b3);
}

vector<double> linspace(double xl, double xh, int n)
{
  vector<double> res(n,0);
  for(int i=0;i<n;++i)
    res[i] = xl + (xh-xl)*(double)i/(double)(n-1);
  return res;
}

double min(vector<double> const& v)
{
  double res = v[0];
  for(int i=1;i<(int)v.size();++i)
    res = min(res,v[i]);
  return res;
}

double max(vector<double> const& v)
{
  double res = v[0];
  for(int i=1;i<(int)v.size();++i)
    res = max(res,v[i]);
  return res;
}
