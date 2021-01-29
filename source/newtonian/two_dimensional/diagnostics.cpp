#include "../../misc/simple_io.hpp"
#include "../../misc/lazy_list.hpp"
#include "diagnostics.hpp"

using std::cout;
using std::fstream;
using std::ios;

void DisplayError(UniversalError const& eo)
{
  cout.precision(16);
  cout << eo.GetErrorMessage() << endl;
  for(size_t i=0;i<eo.GetFields().size();++i)
    cout << eo.GetFields()[i] << " = "<< eo.GetValues()[i] << endl;
}

void write_error(const string& fname,
		 const UniversalError& eo)
{
  ofstream f(fname.c_str());
  f << eo.GetErrorMessage() << endl;
  for(size_t i=0;eo.GetFields().size();++i)
    f << eo.GetFields()[i] << " = " << eo.GetValues()[i] << endl;
  f.close();
}

vector<Vector2D> ReadVector2DFromFile(string filename)
{
  fstream myFile (filename.c_str(),ios::in | ios::binary);
  if(!myFile.good())
    throw UniversalError("Error opening Vector2D file!!");
  int N;
  myFile.read(reinterpret_cast<char*>(&N),sizeof (int));
  vector<Vector2D> res(static_cast<size_t>(N));
  for(int i=0;i<N;++i)
    {
      double x,y;
      myFile.read(reinterpret_cast<char*>(&x),sizeof(double));
      myFile.read(reinterpret_cast<char*>(&y),sizeof(double));
      res[static_cast<size_t>(i)]=Vector2D(x,y);
    }
  myFile.close();
  return res;
}

void WriteVector2DToFile(vector<Vector2D> const& vec,string filename)
{
  if(vec.empty())
    throw UniversalError("Attempted to write a vector of Vector2D to file with zero length");
  fstream myFile (filename.c_str(),ios::out | ios::binary);
  int n=static_cast<int>(vec.size());
  myFile.write (reinterpret_cast<char*>(&n),sizeof(int));
  for(int i=0;i<n;++i)
    {
      myFile.write (reinterpret_cast<const char*>(&vec[static_cast<size_t>(i)].x),
		    sizeof(double));
      myFile.write (reinterpret_cast<const char*>(&vec[static_cast<size_t>(i)].y),
		    sizeof(double));
    }
  myFile.close();
}
