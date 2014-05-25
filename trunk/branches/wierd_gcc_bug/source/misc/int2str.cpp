#include <sstream>
#include "int2str.hpp"

string int2str(int n)
{
	stringstream ss;
	ss << n;
	return ss.str();
}
