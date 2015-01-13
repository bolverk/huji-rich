#include "kill_switch.hpp"
#include <fstream>

using std::ofstream;
using std::ifstream;
using std::endl;

KillSwitch::KillSwitch(const string& fname,
		       TerminationCondition& tc):
  fname_(fname), tc_(tc)
{
  ofstream f(fname.c_str());
  f << '0' << endl;
  f.close();
}

bool KillSwitch::operator()(const hdsim& sim)
{
  char buf;
  ifstream f(fname_.c_str());
  f >> buf;
  return (buf=='0')&&(tc_(sim));
}
