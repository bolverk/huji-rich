#include "polygon_overlap_area.hpp"

PolygonOverlap::PolygonOverlap(void):gen(0){}

double PolygonOverlap::PolyArea(vector<Vector2D> const& polygon)
{
	int ntriangles=(int)polygon.size()-2;
	if(ntriangles<1)
		return 0;
	double res=0;
	for(int i=0;i<ntriangles;++i)
		res+=abs(CrossProduct(polygon[i+1]-polygon[0],polygon[i+2]
			-polygon[0]));
	return 0.5*res;
}

double PolygonOverlap::polygon_overlap_area(vector<Vector2D> const& ch1,
	vector<Vector2D> const& ch2,double R0,double R1)
{
	vector<Vector2D> p(ch1),q(ch2);
	double r[2];
	int n=(int)p.size();
	boost::random::uniform_real_distribution<> dist(-1,1);
	for(int i=0;i<n;++i)
	{
		r[0]=dist(gen);
		r[1]=dist(gen);
		p[i].Set(p[i].x+r[0]*R0,p[i].y+r[1]*R0);
	}
	n=(int)q.size();
	for(int i=0;i<n;++i)
	{
		r[0]=dist(gen);
		r[1]=dist(gen);
		q[i].Set(q[i].x+r[0]*R1,q[i].y+r[1]*R1);
	}
	vector<Vector2D> poly=ConvexIntersect(p,q);
	return PolyArea(poly);
}
