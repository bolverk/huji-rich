#include "polygon_overlap_area.hpp"

PolygonOverlap::PolygonOverlap(void):gen(0){}

double PolygonOverlap::PolyArea(vector<Vector2D> const& polygon)
{
  int ntriangles=static_cast<int>(polygon.size())-2;
	if(ntriangles<1)
		return 0;
	double res=0;
	for(size_t i=0;i<static_cast<size_t>(ntriangles);++i)
		res+=fabs(CrossProduct(polygon[i+1]-polygon[0],polygon[i+2]
			-polygon[0]));
	return 0.5*res;
}

double PolygonOverlap::polygon_overlap_area(vector<Vector2D> const& ch1,
	vector<Vector2D> const& ch2,double R0,double R1)
{
	vector<Vector2D> p(ch1),q(ch2);
	double r[2];
	int n=static_cast<int>(p.size());
	boost::random::uniform_real_distribution<> dist(-1,1);
	for(size_t i=0;i<static_cast<size_t>(n);++i)
	{
		r[0]=2*dist(gen)-1;
		r[1]=2*dist(gen)-1;
		p[i].Set(p[i].x+r[0]*R0,p[i].y+r[1]*R0);
	}
	n=static_cast<int>(q.size());
	for(size_t i=0;i<static_cast<size_t>(n);++i)
	{
		r[0]=2*dist(gen)-1;
		r[1]=2*dist(gen)-1;
		q[i].Set(q[i].x+r[0]*R1,q[i].y+r[1]*R1);
	}
	vector<Vector2D> poly=ConvexIntersect(p,q);
	return PolyArea(poly);
}
