#include "simple_io.hpp"
#include "universal_error.hpp"

void write_number(double num,
		  string const& fname,
		  int prec)
{
  ofstream f(fname.c_str());
  f.precision(prec);
  f << num << endl;
  f.close();
}

void write_vector(vector<double> const& v,
		  string const& fname,
		  int prec)
{
  ofstream f(fname.c_str());
  f.precision(prec);
  for(int i=0;i<(int)v.size();++i)
    f << v[i] << endl;
  f.close();
}

void write_vector(vector<int> const& v,
		  string const& fname)
{
  ofstream f(fname.c_str());
  for(int i=0;i<(int)v.size();++i)
    f << v[i] << endl;
  f.close();
}

double read_number(string const& fname)
{
  double buf = 0;
  ifstream f(fname.c_str());
  if(!f)
    throw UniversalError("Could not find file "+fname);
  f >> buf;
  f.close();
  return buf;
}

int read_int(string const& fname)
{
  int buf = 0;
  ifstream f(fname.c_str());
  if(!f)
    throw UniversalError("Could not find file "+fname);
  f >> buf;
  f.close();
  return buf;
}
