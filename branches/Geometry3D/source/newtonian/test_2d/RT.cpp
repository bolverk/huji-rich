#include "RT.hpp"

RT_velocityY::RT_velocityY()
{
}

RT_velocityY::~RT_velocityY()
{
}

double RT_velocityY::operator()(Vector2D const& r) const
{
	return 0.0025*(1-cos(12.566370614359*r.x))*
		(1-cos(12.566370614359*r.y/3));
}

RT_Pressure::RT_Pressure(double g,double rhoup,double rhodown):
  g_(g),
  rhou_(rhoup),
  rhod_(rhodown)
{
  //	g_=g;
  //	rhou_=rhoup;
  //	rhod_=rhodown;
}

RT_Pressure::~RT_Pressure()
{
}

double RT_Pressure::operator()(Vector2D const& r) const
{
	if(r.y>0.75)
		return 2.5+g_*(r.y-0.75)*rhou_;
	else
		return 2.5+g_*(r.y-0.75)*rhod_;
}
