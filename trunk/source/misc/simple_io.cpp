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
  for(size_t i=0;i<v.size();++i)
    f << v[i] << "\n";
  f.close();
}

void write_vector(vector<int> const& v,
		  string const& fname)
{
  ofstream f(fname.c_str());
  for(size_t i=0;i<v.size();++i)
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

void binary_write_single_int(int n, ofstream& fh)
{
  fh.write((const char*)&n,sizeof(int));
}

void binary_write_single_double(double d, ofstream& fh)
{
  fh.write((const char*)&d,sizeof(double));
}