#include "simple_io.hpp"
#include "universal_error.hpp"
#include <cassert>

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

void write_vector(vector<size_t> const& v,
	string const& fname)
{
	ofstream f(fname.c_str());
	for (size_t i = 0; i<v.size(); ++i)
		f << v[i] << endl;
	f.close();
}
namespace {
  bool missing_file_data(string const& fname)
  {
    std::cout << "Could not find file " << fname << std::endl;
    return false;
  }
}

vector<double> read_vector(string const& fname)
{
	double buf = 0;
	vector<double> res;
	ifstream f(fname.c_str());
	assert(f || missing_file_data(fname));
	while (f >> buf)
		res.push_back(buf);
	f.close();
	return res;
}

double read_number(string const& fname)
{
  double buf = 0;
  ifstream f(fname.c_str());
  assert(f || missing_file_data(fname));
  f >> buf;
  f.close();
  return buf;
}

int read_int(string const& fname)
{
  int buf = 0;
  ifstream f(fname.c_str());
  assert(f || missing_file_data(fname));
  f >> buf;
  f.close();
  return buf;
}

void binary_write_single_int(int n, ofstream& fh)
{
  fh.write(reinterpret_cast<const char*>(&n),sizeof(int));
}

void binary_write_single_double(double d, ofstream& fh)
{
  fh.write(reinterpret_cast<const char*>(&d),sizeof(double));
}

void binary_write_single_size_t(size_t n,ofstream& fh)
{
  fh.write(reinterpret_cast<const char*>(&n),sizeof(size_t));
}
