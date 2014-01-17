#include <cmath>
#include "bisection.hpp"
#include "universal_error.hpp"

double find_upper_bracket(Func1Var const& f,
			  double xl)
{
  const int max_iter = 10;
  int iter = 0;

  double fl = f.eval(xl);
  double xr = xl;
  // TODO: change to fr = fl
  double fr = f.eval(xr);
  while(fl*fr>=0){
    xr = 2*xr;
    fr = f.eval(xr);
    if(iter>max_iter)
      throw UniversalError("Error in bisection::find_upper_bracket. Max number of iterations exceeded");
    else
      ++iter;
  }
  return 2*xr;
}

double bisection(Func1Var const& f,
		 double xl,
		 double xr,
		 double tol)
{
  const int max_iter = 100;

  int iter = 0;

  double fl = f.eval(xl);
  if(fl==0)
    return xl;
  
  double fr = f.eval(xr);
  if(fr==0)
    return xr;

  if(fl*fr>0){
    UniversalError eo("Error in bisection: root is not bracketed");
    eo.AddEntry("xl",xl);
    eo.AddEntry("xr",xr);
    eo.AddEntry("fl",fl);
    eo.AddEntry("fr",fr);
    throw eo;
  }

  //  double xm = 0.5*(xl+xm);
  while(abs(xl-xr)>tol){
    double xm = 0.5*(xl+xr);
    double fm = f.eval(xm);
    if(fm==0)
      return xm;
    if(fm*fl>0)
      xl = xm;
    else if(fm*fr>0)
      xr = xm;
    else
      throw UniversalError("Something bad happened. Probably a NaN");

    if(iter>max_iter)
      throw UniversalError("Error in Bisection: Maximum number of iterations exceeded");
    else
      ++iter;
  }

  return 0.5*(xl+xr);
}
