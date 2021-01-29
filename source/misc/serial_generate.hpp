#ifndef SERIAL_GENERATE_HPP
#define SERIAL_GENERATE_HPP

#include <vector>
#include <functional>

using std::vector;
using std::function;

template<class S, class T> vector<T> serial_generate
(const vector<S>& source,
 function<T(const S&)> func)
{
  vector<T> res(source.size());
  transform(source.begin(), source.end(),
	    res.begin(),
	    func);
  /*
  for(size_t i=0;i<source.size();++i)
    res.at(i) = func(source.at(i));
  */
  return res;
}

#endif // SERIAL_GENERATE_HPP
