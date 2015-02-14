#include "PolyIntersect.hpp"

using std::vector;
using std::max;

IntersectFlags SegmentIntersection(Vector2D const& p0,Vector2D const& p1,
	Vector2D const& q0,Vector2D const& q1,Vector2D &Intersection)
{
	/*int n=static_cast<int>(p.size());
	int m=static_cast<int>(q.size());
	Vector2D p0(p[ploc]);
	Vector2D p1(p[(ploc+1)%n]);
	Vector2D q0(q[qloc]);
	Vector2D q1(q[(qloc+1)%m]);*/
	if(min(p0.x,p1.x)>max(q0.x,q1.x)||min(q1.x,q0.x)>max(p0.x,p1.x)||
		min(p0.y,p1.y)>max(q1.y,q0.y)||min(q1.y,q0.y)>max(p0.y,p1.y))
		return False;
	boost::array<Vector2D,3> points;
	points[0]=p0;
	points[1]=p1;
	points[2]=q0;
	double temp=orient2d(TripleConstRef<Vector2D>(p0,p1,q0));
	points[2]=q1;
	if(temp*orient2d(TripleConstRef<Vector2D>(p0,p1,q1))>0)
		return False;
	points[0]=q0;
	points[1]=q1;
	points[2]=p0;
	temp=orient2d(TripleConstRef<Vector2D>(q0,q1,p0));
	points[2]=p1;
	if(temp*orient2d(TripleConstRef<Vector2D>(q0,q1,p1))>0)
		return False;
	points[0]=Vector2D(0,0);
	points[1]=p1-p0;
	points[2]=q1-q0;
	double d=orient2d(TripleConstRef<Vector2D>(Vector2D(0,0),
						   p1-p0,
						   q1-q0));
	double eps = 1e-7*abs(p1 - p0)*abs(q1 - q0);
	if(fabs(d)<0)
		return Par;
	double xi=((q0.x-q1.x)*(p0.x*p1.y-p1.x*p0.y)-(p0.x-p1.x)*
		(q0.x*q1.y-q1.x*q0.y))/d;
	double yi=((q0.y-q1.y)*(p0.x*p1.y-p1.x*p0.y)-(p0.y-p1.y)*
		(q0.x*q1.y-q1.x*q0.y))/d;
	Intersection.Set(xi,yi);
	eps=1e-7*sqrt(abs(p1-p0)*abs(q1-q0));
	if((xi+eps)<min(p0.x,p1.x)||(xi-eps)>max(p0.x,p1.x))
		return Par;
	if((xi+eps)<min(q0.x,q1.x)||(xi-eps)>max(q0.x,q1.x))
		return Par;
	if((yi+eps)<min(p0.y,p1.y)||(yi-eps)>max(p0.y,p1.y))
		return Par;
	if((yi+eps)<min(q0.y,q1.y)||(yi-eps)>max(q0.y,q1.y))
		return Par;
	return True;
}

/*
vector<Vector2D> GetParEdge(Vector2D const& p0,Vector2D const& p1,
	Vector2D const& q0,Vector2D const& q1)
{
	vector<Vector2D> res(2);
	if((p1.x-p0.x)*(p1.y-p0.y)>0)
	{
		res[0].Set(min(max(p0.x,p1.x),max(q0.x,q1.x)),
			min(max(p0.y,p1.y),max(q0.y,q1.y)));
		res[1].Set(max(min(p0.x,p1.x),min(q0.x,q1.x)),
			max(min(p0.y,p1.y),min(q0.y,q1.y)));
	}
	else
	{
		res[0].Set(max(min(p0.x,p1.x),min(q0.x,q1.x)),
			min(max(p0.y,p1.y),max(q0.y,q1.y)));
		res[1].Set(min(max(p0.x,p1.x),max(q0.x,q1.x)),
			max(min(p0.y,p1.y),min(q0.y,q1.y)));
	}
	return res;
}
*/

vector<Vector2D> ConvexIntersect(vector<Vector2D> const& poly0,vector<Vector2D>
	const& poly1)
{
	// poly0 is a poly1 is b
	vector<Vector2D> res;
	InFlags flag=UnKnown;
	bool FirsPoint=true;
	const int n=static_cast<int>(poly0.size());
	const int m=static_cast<int>(poly1.size());
	int p0index=0,p1index=0;
	boost::array<Vector2D,3> AreaSignCheck;
	int p0counter=0,p1counter=0;
	do
	{
	  AreaSignCheck[0]=poly0[static_cast<size_t>((p0index+n-1)%n)];
	  AreaSignCheck[1]=poly0[static_cast<size_t>(p0index)];
	  AreaSignCheck[2]=poly1[static_cast<size_t>(p1index)];
	  double bLeftOfa=orient2d
	    (TripleConstRef<Vector2D>
	     (poly0[static_cast<size_t>((p0index+n-1)%n)],
	      poly0[static_cast<size_t>(p0index)],
	      poly1[static_cast<size_t>(p1index)]));
	  AreaSignCheck[0]=poly1[static_cast<size_t>((p1index+m-1)%m)];
	  AreaSignCheck[1]=poly1[static_cast<size_t>(p1index)];
	  AreaSignCheck[2]=poly0[static_cast<size_t>(p0index)];
		double aLeftOfb=orient2d
		  (TripleConstRef<Vector2D>
		   (poly1[static_cast<size_t>((p1index+m-1)%m)],
		    poly1[static_cast<size_t>(p1index)],
		    poly0[static_cast<size_t>(p0index)]));

		Vector2D intersect;
		IntersectFlags inter=SegmentIntersection(poly0[static_cast<size_t>((p0index+n-1)%n)],
							 poly0[static_cast<size_t>(p0index)],
							 poly1[static_cast<size_t>(p1index)],
							 poly1[static_cast<size_t>((p1index+m-1)%m)],intersect);
		if(inter==True)
		{
			res.push_back(intersect);
			if(aLeftOfb>0)
				flag=Pi;
			else
				flag=Qi;
			if(FirsPoint)
			{
				FirsPoint=false;
				p0counter=0;
				p1counter=0;
			}
		}

		AreaSignCheck[0]=Vector2D(0,0);
		AreaSignCheck[1]=poly0[static_cast<size_t>(p0index)]-poly0[static_cast<size_t>((p0index+n-1)%n)];
		AreaSignCheck[2]=poly1[static_cast<size_t>(p1index)]-poly1[static_cast<size_t>((p1index+m-1)%m)];
		double areasign=orient2d
		  (TripleConstRef<Vector2D>
		   (Vector2D(0,0),
		    poly0[static_cast<size_t>(p0index)]-poly0[static_cast<size_t>((p0index+n-1)%n)],
		    poly1[static_cast<size_t>(p1index)]-poly1[static_cast<size_t>((p1index+m-1)%m)]));
		if(inter==Par)
		{
		  if(ScalarProd(poly0[static_cast<size_t>(p0index)]-poly0[static_cast<size_t>((p0index+n-1)%n)],
				poly1[static_cast<size_t>(p1index)]-poly1[static_cast<size_t>((p1index+m-1)%m)])<0)
			{
				// Do we need this? The area will be zero
				/*res=GetParEdge(poly0[p0index],poly0[(p0index+n-1)%n],
					poly1[p1index],poly1[(p1index+m-1)%m]);*/
				return vector<Vector2D> ();
			}
		}
		if((fabs(areasign)<1e-9)&&(aLeftOfb<0)&&(bLeftOfa<0))
		  return vector<Vector2D> ();
		if((fabs(areasign)<1e-9)&&(fabs(aLeftOfb)<1e-9)&&(fabs(bLeftOfa)<1e-9))
		{
			if(flag==Pi)
			{
				++p1counter;
				p1index=(p1index+1)%m;
			}
			else
			{
				++p0counter;
				p0index=(p0index+1)%n;
			}
		}
		if(areasign>=0)
		{
			if(bLeftOfa>0)
			{
				if(flag==Pi)
				  res.push_back(poly0[static_cast<size_t>(p0index)]);
				++p0counter;
				p0index=(p0index+1)%n;
			}
			else
			{
				if(flag==Qi)
				  res.push_back(poly1[static_cast<size_t>(p1index)]);
				++p1counter;
				p1index=(p1index+1)%m;
			}
		}
		else
		{
			if(aLeftOfb>0)
			{
				if(flag==Qi)
				  res.push_back(poly1[static_cast<size_t>(p1index)]);
				++p1counter;
				p1index=(p1index+1)%m;
			}
			else
			{
				if(flag==Pi)
				  res.push_back(poly0[static_cast<size_t>(p0index)]);
				++p0counter;
				p0index=(p0index+1)%n;
			}
		}
	}
	while(((p0counter<n)||(p1counter<m))&&(p0counter<2*n)&&(p1counter<2*m));
	return res;
}
