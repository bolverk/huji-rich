#include "ReducedPrimitiveGradient2D.hpp"
using namespace std;

ReducedPrimitiveGradient2D::ReducedPrimitiveGradient2D(void):
density(),pressure(),xvelocity(),yvelocity(),tracers(vector<Vector2D> ()) {}

ReducedPrimitiveGradient2D::ReducedPrimitiveGradient2D(Vector2D const& d,
	Vector2D const& p,Vector2D const& vx,Vector2D const& vy,vector<Vector2D>
	const& trace=vector<Vector2D>()):
density(d),pressure(p),xvelocity(vx),yvelocity(vy),tracers(trace) {}

ReducedPrimitiveGradient2D& ReducedPrimitiveGradient2D::operator+=
	(ReducedPrimitiveGradient2D const& source)
{
	density += source.density;
	pressure += source.pressure;
	xvelocity += source.xvelocity;
	yvelocity += source.yvelocity;

	if(source.tracers.empty()&&!tracers.empty())
	  {
	    throw "Something bad happened in ReducedPrimitiveGradient2D::operator+=";
	  }
	if(!tracers.empty())
	  {
	    if(tracers.size()!=source.tracers.size())
	      throw "Error in ReducedPrimitiveGradient2D::operator+ tracers have different size";
	  }
	if(tracers.empty()&&!source.tracers.empty()){
	  tracers.resize(source.tracers.size());
	}

	if(!tracers.empty())
	{
	  int n=(int)tracers.size();
	  for(int i=0;i<n;++i)
	    tracers[i]+=source.tracers[i];
	}
	return *this;
}

ReducedPrimitiveGradient2D& ReducedPrimitiveGradient2D::operator*=
	(double s)
{
	density *= s;
	pressure *= s;
	xvelocity *= s;
	yvelocity *= s;

	if(!tracers.empty())
	{
		int n=(int)tracers.size();
		for(int i=0;i<n;++i)
			tracers[i]*=s;
	}
	return *this;
}
