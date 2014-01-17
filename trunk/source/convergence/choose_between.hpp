// Pure virtual
#ifndef CHOOSE_BETWEEN_HPP
#define CHOOSE_BETWEEN_HPP 1

template<class T> T& choose_between(bool flag,
				    T& tres,
				    T& fres)
{
  if(flag)
    return tres;
  else
    return fres;
}

#endif // CHOOSE_BETWEEN_HPP
