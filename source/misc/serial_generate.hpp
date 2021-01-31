#ifndef SERIAL_GENERATE_HPP
#define SERIAL_GENERATE_HPP

#include <vector>
#include <functional>
#include <numeric>

using std::vector;
using std::function;

template<class S, class T> vector<T> serial_generate
(const vector<S>& source,
 function<T(const S&)> func)
{
  vector<T> res(source.size());
  transform(source.begin(),
	    source.end(),
	    res.begin(),
	    func);
  return res;
}

template<class S1, class S2, class T> vector<T> serial_generate
(const vector<S1>& s1,
 const vector<S2>& s2,
 function<T(const S1&, const S2&)> func)
{
  assert(s1.size()==s2.size());
  vector<T> res(s1.size());
  transform(s1.begin(),
	    s1.end(),
	    s2.begin(),
	    res.begin(),
	    func);
  return res;
}

template<class T> vector<T> diff(const vector<T>& source){
  vector<T> res(source.size()-1);
  transform(next(source.begin(),1),
	    source.end(),
	    source.begin(),
	    res.begin(),
	    [&](const T& t1, const T& t2)
	    {return t1-t2;});
  return res;
}

template<class S> vector<S> create_range(const S& start, const S& length)
{
  vector<S> res(static_cast<size_t>(length));
  iota(res.begin(), res.end(), start);
  return res;
}

template<class T> vector<T> run_along_index(size_t num,
					    function<T(size_t i)> func)
{
  vector<T> res(num);
  for(size_t i=0;i<num;++i)
    res.at(i) = func(i);
  return res;
}

#endif // SERIAL_GENERATE_HPP
