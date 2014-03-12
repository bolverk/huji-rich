#include <sstream>
#include "int2str.hpp"

using namespace std;

string int2str(int n)
{
	stringstream ss;
	ss << n;
	return ss.str();
}
