#include "index2filename.hpp"
#include "../../misc/int2str.hpp"

Rubric::Rubric(const string& prefix,
	       const string& postfix):
  prefix_(prefix), postfix_(postfix) {}

string Rubric::operator()(int index)
{
  return prefix_+int2str(index)+postfix_;
}
