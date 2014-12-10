#include "Delaunay.hpp"
#include <vector>
#include <cmath>
#include "../misc/triplet.hpp"

namespace {
pair<int,int> find_diff(const facet& f1,const facet& f2)
{
  if(f1.vertices.first!=f2.vertices.first &&
     f1.vertices.first!=f2.vertices.second &&
     f1.vertices.first!=f2.vertices.third)
    return pair<int,int>(f1.vertices.first,0);
  else if(f1.vertices.second!=f2.vertices.first &&
	  f1.vertices.second!=f2.vertices.second &&
	  f1.vertices.second!=f2.vertices.third)
    return pair<int,int>(f1.vertices.second,1);
  else if(f1.vertices.third!=f2.vertices.first &&
	  f1.vertices.third!=f2.vertices.second &&
	  f1.vertices.third!=f2.vertices.third)
    return pair<int,int>(f1.vertices.third,2);
  else
    throw UniversalError("Delaunay, Couldn't find difference bewteen two facets");
}
}

Delaunay::DataOnlyForBuild::DataOnlyForBuild():insert_order(vector<int> ()),
	copied(vector<vector<char> > ())
{}

Delaunay::DataOnlyForBuild::DataOnlyForBuild(DataOnlyForBuild const& other):
insert_order(other.insert_order),copied(other.copied){}

Delaunay::DataOnlyForBuild& Delaunay::DataOnlyForBuild::operator=
	(DataOnlyForBuild const& other)
{
	if(this != &other)
	{
		copied=other.copied;
		insert_order=other.insert_order;
	}
	return *this;
}

Delaunay::Delaunay(void):
lastFacet(0),CalcRadius(false),
	radius(vector<double>()),cell_points(vector<Vector2D> ()),
	PointWasAdded(false),
	last_facet_added(0),
	f(vector<facet>()),
	cor(vector<Vector2D>()),
	length(0),
	olength(0),location_pointer(0), last_loc(0),
	logger(0) {}

Delaunay::Delaunay(Delaunay const& other):
lastFacet(other.lastFacet),
	CalcRadius(other.CalcRadius),
	radius(other.radius),cell_points(other.cell_points),
	PointWasAdded(other.PointWasAdded),
	last_facet_added(other.last_facet_added),
	f(other.f),
	cor(other.cor),
	length(other.length),
	olength(other.olength),
	location_pointer(other.location_pointer),
	last_loc(other.last_loc),
	logger(other.logger) {}

Delaunay::~Delaunay(void)
{
	cor.clear();
	f.clear();
	cell_points.clear();
}

namespace
{
	// Checks if a point is inside a triangle
	bool InTriangle(const TripleConstRef<Vector2D>& tri,
			const Vector2D& point)
	{
	  return (orient2d(TripleConstRef<Vector2D>(tri.first,
						   tri.second,
						    point))>0) &&
	    (orient2d(TripleConstRef<Vector2D>(tri.second,
					       tri.third,
					       point))>0) &&
	    (orient2d(TripleConstRef<Vector2D>(tri.third,
					       tri.first,
					       point))>0);
	}

	// Assume cell is orederd in convexhull counterclockwise
	bool InCell(vector<Vector2D> const& points,Vector2D const& p)
	{
		int n=(int)points.size();
		for(int i=0;i<n;++i)
		{
			if(CrossProduct(points[(size_t)i]-p,points[(size_t)((i+1)%n)]-p)<0)
				return false;
		}
		return true;
	}

	vector<double> CellSize(vector<Vector2D> const& points)
	{
		int n=(int)points.size();
		double minx=points[0].x;
		double miny=points[0].y;
		double maxx=minx;
		double maxy=miny;
		for(int i=1;i<n;++i)
		{
			minx=min(points[(size_t)i].x,minx);
			miny=min(points[(size_t)i].y,miny);
			maxx=max(points[(size_t)i].x,maxx);
			maxy=max(points[(size_t)i].y,maxy);
		}
		vector<double> res(4);
		res[0]=minx;
		res[1]=maxx;
		res[2]=miny;
		res[3]=maxy;
		return res;
	}

	int find_index(facet const& fc, int i)
	{
		for(int j=0;j<3;++j)
		{
			if(fc.neighbors[(size_t)j]==i)
				return j;
		}
		throw UniversalError("Error in find_index: Index not found");
	}
}

void Delaunay::add_point(size_t index)
{
	// Check if point is inside big triangle
  assert(InTriangle(TripleConstRef<Vector2D>(cor[olength],
					     cor[olength+1],
					     cor[olength+2]),
		    cor[index]));
  /*
  const TripleConstRef<Vector2D> tocheck(cor[olength],
					 cor[olength+1],
					 cor[olength+2]);
	if(!InTriangle(tocheck,cor[(size_t)index]))
	{
		UniversalError eo("Point not inside large triangle of Delaunay");
		eo.AddEntry("Point x",cor[(size_t)index].x);
		eo.AddEntry("Point y",cor[(size_t)index].y);
		eo.AddEntry("Big tirangle point 1 x",tocheck.first.x);
		eo.AddEntry("Big tirangle point 1 y",tocheck.first.y);
		eo.AddEntry("Big tirangle point 2 x",tocheck.second.x);
		eo.AddEntry("Big tirangle point 2 y",tocheck.second.y);
		eo.AddEntry("Big tirangle point 3 x",tocheck.third.x);
		eo.AddEntry("Big tirangle point 3 y",tocheck.third.y);
		throw eo;
	}
  */
	const size_t triangle=Walk(index);
	facet f_temp;
	const Triplet<int> outer(f[triangle].vertices);
	const Triplet<int> temp_friends(f[triangle].neighbors);
	// create and _update the new facets
	//	f.push_back(f_temp);
	//	f.push_back(f_temp);
	f[triangle].vertices.set(outer.third,outer.first,index);
	f[triangle].neighbors.set(temp_friends.third,
				  location_pointer+1,
				  location_pointer+2);
	f.push_back(facet(TripleConstRef<int>(outer.first,
					      outer.second,
					      index),
			  TripleConstRef<int>(temp_friends.first,
					      location_pointer+2,
					      triangle)));
	f.push_back(facet(TripleConstRef<int>(outer.second,
					      outer.third,
					      index),
			  TripleConstRef<int>(temp_friends.second,
					      triangle,
					      location_pointer+1)));
	// _update the friends list of the friends
	if(temp_friends.second!=last_loc)
	{
		const int i=find_index(f[(size_t)temp_friends.second],triangle);
		f[(size_t)temp_friends.second].neighbors[(size_t)i] = location_pointer+2;
	}
	if(temp_friends.first!=last_loc)
	{
	  const int i=find_index(f[(size_t)temp_friends.first],triangle);
		f[(size_t)temp_friends.first].neighbors[(size_t)i] = location_pointer+1;
	}
	// Calculate radius if needed
	if(CalcRadius)
	{
		radius[(size_t)triangle]=CalculateRadius(triangle);
		int n=int(f.size());
		int m=int(radius.size());
		if(n>m-1)
		{
			radius.push_back(CalculateRadius(location_pointer+1));
			radius.push_back(CalculateRadius(location_pointer+2));
		}
		else
			if(n>m)
			{
				radius[(size_t)location_pointer+1]=CalculateRadius(location_pointer+1);
				radius.push_back(CalculateRadius(location_pointer+2));
			}
			else
			{
				radius[(size_t)location_pointer+1]=CalculateRadius(location_pointer+1);
				radius[(size_t)location_pointer+2]=CalculateRadius(location_pointer+2);
			}
	}

	// check if flipping is needed
	flip(triangle,temp_friends.third);
	flip(location_pointer+1,temp_friends.first);
	flip(location_pointer+2,temp_friends.second);

	// _update number of facets
	location_pointer+=2;
}

void Delaunay::flip(size_t i, size_t j)
{
  stack<std::pair<size_t,size_t> > flip_stack
    (std::deque<std::pair<size_t,size_t> >(1,std::pair<size_t,size_t>(i,j)));
  //  flip_stack.push(std::pair<int,int>(i,j));
	while(!flip_stack.empty())
	{
	  if(flip_stack.top().second==(size_t)last_loc)
			flip_stack.pop();
		else
		{
		  const pair<size_t,size_t> indexes = flip_stack.top();
			// Returns the index to the point to check in coordinates and the index of the point in the facet
		  const pair<int,int> check = find_diff(f[indexes.second],
					  f[indexes.first]);
		  const pair<int,int> other = find_diff(f[indexes.first],
					  f[indexes.second]);

		  facet& prefetch_1 = f[indexes.first];
			if(incircle(cor[(size_t)prefetch_1.vertices.first],
				    cor[(size_t)prefetch_1.vertices.second],
				    cor[(size_t)prefetch_1.vertices.third],
				    cor[(size_t)check.first])>0)
			{
				//The point is in a circle change the facets and their friends
				const int v1=prefetch_1.vertices[(size_t)(other.second+1)%3];
				const int f1=prefetch_1.neighbors[(size_t)other.second];
				const int f12=prefetch_1.neighbors[(size_t)(other.second+2)%3];
				facet& prefetch_2 = f[indexes.second];
				const int v2=prefetch_2.vertices[(size_t)(check.second+1)%3];
				const int f2=prefetch_2.neighbors[(size_t)(check.second+2)%3];
				const int f22=prefetch_2.neighbors[(size_t)check.second];
				prefetch_1.vertices.set(other.first,v1,check.first);
				prefetch_2.vertices.set(check.first,v2,other.first);
				prefetch_1.neighbors.set(f1,f2,indexes.second);
				prefetch_2.neighbors.set(f22,f12,indexes.first);
				// change the friends of the friends if needed
				if(f2!=last_loc)
				{
					f[(size_t)f2].neighbors[(size_t)find_index(f[(size_t)f2],indexes.second)] = indexes.first;
				}
				if(f12!=last_loc)
				{
					f[(size_t)f12].neighbors[(size_t)find_index(f[(size_t)f12],indexes.first)] = indexes.second;
				}
				// Calculate the new radius if needed
				if(CalcRadius)
				{
				  radius[indexes.first]=CalculateRadius(indexes.first);
				  radius[indexes.second]=CalculateRadius(indexes.second);
				}
				// clear the checked facets
				flip_stack.pop();
				// push into the stack the new facets to check
				flip_stack.push(std::pair<size_t,size_t>(indexes.second,
									 prefetch_2.neighbors.first));
				flip_stack.push(std::pair<size_t,size_t>(indexes.first,
									 prefetch_1.neighbors.second));
			}
			else
			{
				// clear the checked facets
				flip_stack.pop();
			}
		}
	}
}

void Delaunay::build_delaunay(vector<Vector2D>const& vp,vector<Vector2D> const& cpoints)
{
	cell_points=cpoints;
	DataOnlyForBuild data;
	lastFacet=0;
	CalcRadius=false;
	length=int(vp.size()+3);
	int len=length-3;
	olength=len;
	f.clear();
	cor.clear();
	f.reserve((size_t)(2*length+1+(int)(17*sqrt(1.*length))));
	cor.reserve((size_t)(length+9*(int)sqrt(1.*length)));
	last_loc=INT_MAX;
	for(int i=0;i<len;i++)
	{
		cor.push_back(vp[(size_t)i]);
	}
	// Check point input
	CheckInput();

	data.insert_order=HilbertOrder(cor,olength,0);

	// add the 3 extreme points
	Vector2D p_temp;
	vector<double> cellsize=CellSize(cell_points);
	double width=cellsize[1]-cellsize[0];
	double height=cellsize[3]-cellsize[2];
	width=max(width,height);
	height=max(width,height);
	p_temp.x = cellsize[0]-100*width;
	p_temp.y = cellsize[2]-100*height;
	cor.push_back(p_temp);
	p_temp.x = cellsize[1]+100*width;
	p_temp.y = cellsize[2]-100*height;
	cor.push_back(p_temp);
	p_temp.x = (cellsize[0]+cellsize[1])/2.0;
	p_temp.y = cellsize[3]+100*height;
	cor.push_back(p_temp);
	// Create the big triangle, and assign friends
	facet f_temp;
	f.push_back(f_temp);
	f[0].vertices[0] = len;
	f[0].vertices[1] = len+1;
	f[0].vertices[2] = len+2;
	for(size_t i=0;i<3;i++)
		f[0].neighbors[i] = last_loc;
	location_pointer=0;
	// add the points
	for(int i=0;i<length-3;i++)
	{
		add_point(data.insert_order[(size_t)i]);
	}
	// Calculate radius
	radius.resize(f.size());
	int n=int(f.size());
	for(int i=0;i<n;++i)
		radius[(size_t)i]=CalculateRadius(i);
	CalcRadius=true;
}

double Delaunay::triangle_area(int index)
{
  const TripleConstRef<Vector2D> p
    (cor[(size_t)f[(size_t)index].vertices.first],
     cor[(size_t)f[(size_t)index].vertices.second],
     cor[(size_t)f[(size_t)index].vertices.third]);
  const double x1=p.third.x-p.first.x;
	const double x2=p.second.x-p.first.x;
	const double y1=p.third.y-p.first.y;
	const double y2=p.second.y-p.first.y;
	return -0.5*(x1*y2-x2*y1);
}

void Delaunay::update(const vector<Vector2D>& points,vector<Vector2D>
	const& cpoints)
{
	if(logger)
		logger->output(cor,f);
	build_delaunay(points,cpoints);
}

namespace {

  int Triplet<int>::* walk_condition(const vector<Vector2D>& cor,
				     const Triplet<int>& vertices,
				     size_t point)
  {
    if(orient2d(TripleConstRef<Vector2D>
		(cor[(size_t)vertices.first],
		 cor[(size_t)vertices.second],
		 cor[point]))<0)
      return &Triplet<int>::first;
    else if(orient2d(TripleConstRef<Vector2D>
		     (cor[(size_t)vertices.second],
		      cor[(size_t)vertices.third],
		      cor[point]))<0)
      return &Triplet<int>::second;
    else if(orient2d(TripleConstRef<Vector2D>
		     (cor[(size_t)vertices.third],
		      cor[(size_t)vertices.first],
		      cor[point]))<0)
      return &Triplet<int>::third;
    else
      return 0;
  }

  class WalkBookkeeper
  {
  public:

  private:
    
  };

  size_t find_new_facet(const vector<Vector2D>& cor,
			const vector<facet>& f,
			size_t point,
			size_t last_facet)
  {
    size_t res = last_facet;
    int Triplet<int>::* next = walk_condition(cor,
					 f[res].vertices,
					 point);
    while(next){
      res = f[res].neighbors.*next;
      next = walk_condition(cor,
			    f[res].vertices,
			    point);
    }
    return res;
  }
}

size_t Delaunay::Walk(size_t point)
{
  lastFacet = find_new_facet(cor,f,point,lastFacet);
  return lastFacet;
}

vector<int> Delaunay::FindContainingTetras(int StartTetra, int point)
{
	vector<int> res;
	FindContainingTetras(StartTetra, point, res);
	return res;
}

double Delaunay::FindMaxRadius(int point)
{
	const vector<int> vec = FindContainingTetras(Walk(point),point);
	double r=0;
	for(size_t i=0;i<vec.size();++i)
		r = max(r,radius[(size_t)vec[(size_t)i]]);
	return 2*r;
}

void Delaunay::FindContainingTetras(int StartFacet,int point,vector<int> &result)
{
	result.clear();
	int PointLocation=FindPointInFacet(StartFacet,point);
	int NextFacet=f[(size_t)StartFacet].neighbors[(size_t)PointLocation];
	result.reserve(12);
	result.push_back(NextFacet);
	while(NextFacet!=StartFacet)
	{
		PointLocation=FindPointInFacet(NextFacet,point);
		NextFacet=f[(size_t)NextFacet].neighbors[(size_t)PointLocation];
		result.push_back(NextFacet);
	}
}

int Delaunay::FindPointInFacet(int facet,int point)
{
	for(int i=0;i<3;++i)
		if(f[(size_t)facet].vertices[(size_t)i]==point)
			return i;
	UniversalError eo("Error in Delaunay, FindPointInFacet");
	eo.AddEntry("Facet number",facet);
	eo.AddEntry("Point number",point);
	throw eo;
}

bool Delaunay::IsOuterFacet(int facet)const
{
	//int PointNum=length-1;
	for(int i=0;i<3;++i)
		for(size_t j=0;j<3;++j)
		  if(f[(size_t)facet].vertices[(size_t)i]==(int)(olength+j))
				return true;
	return false;
}

double Delaunay::CalculateRadius(int facet)
{
	const double big=1e10;
	const double a=cor[(size_t)f[(size_t)facet].vertices[0]].distance(cor[(size_t)f[(size_t)facet].vertices[1]]);
	const double b=cor[(size_t)f[(size_t)facet].vertices[0]].distance(cor[(size_t)f[(size_t)facet].vertices[2]]);
	const double c=cor[(size_t)f[(size_t)facet].vertices[2]].distance(cor[(size_t)f[(size_t)facet].vertices[1]]);
	const double temp1=b+c-a;
	if(temp1<=0)
	{
		if(a>big*b||a>big*c) // Do we have a small edge?
			return 0.5*a;
		else
			return 0.5*(b+c); // we have 3 points on a line
	}
	const double temp2=c+a-b;
	if(temp2<=0)
	{
		if(b>big*a||b>big*c) // Do we have a small edge?
			return 0.5*b;
		else
			return 0.5*(a+c); // we have 3 points on a line
	}
	const double temp3=b-c+a;
	if(temp3<=0)
	{
		if(c>big*b||c>big*a) // Do we have a small edge?
			return 0.5*c;
		else
			return 0.5*(b+a); // we have 3 points on a line
	}
	return a*b*c/sqrt((a+b+c)*temp1*temp2*temp3);
}

void Delaunay::CheckInput()
{
	int n=(int)cor.size();
	for(int i=0;i<n;++i)
	{
		if(!InCell(cell_points,cor[(size_t)i]))
		{
			UniversalError eo("Mesh point outside cell");
			eo.AddEntry("Point number",i);
			eo.AddEntry("X location",cor[(size_t)i].x);
			eo.AddEntry("Y location",cor[(size_t)i].y);
			throw eo;
		}
	}
}

int Delaunay::GetOriginalIndex(int NewPoint) const
{
	return NewPoint;
}

double Delaunay::GetFacetRadius(int facet) const
{
	return radius[(size_t)facet];
}

void Delaunay::ChangeOlength(int n)
{
	olength=n;
}

void Delaunay::Changelength(int n)
{
	length=n+3;
}

vector<Vector2D>& Delaunay::ChangeCor(void)
{
	return cor;
}

const vector<Vector2D>& Delaunay::getCor(void) const
{
  return cor;
}

facet* Delaunay::get_facet(int index)
{
	return &f[(size_t)index];
}

double Delaunay::get_facet_coordinate(int Facet,int vertice, int dim)
{
	if(dim==0)
		return cor[(size_t)f[(size_t)Facet].vertices[(size_t)vertice]].x;
	else
		return cor[(size_t)f[(size_t)Facet].vertices[(size_t)vertice]].y;
}

Vector2D Delaunay::get_point(size_t index) const
{
	return cor[index];
}

double Delaunay::get_cor(int index, int dim) const
{
	if(dim==0)
		return cor[(size_t)index].x;
	else if(dim==1)
		return cor[(size_t)index].y;
	else
		throw UniversalError("Error in Delaunay::get_cor. Invalid index");
}

int Delaunay::get_num_facet(void)
{
	return (int)f.size();
}

int Delaunay::get_length(void) const
{
	return length-3;
}

int Delaunay::get_last_loc(void) const
{
	return last_loc;
}

void Delaunay::set_point(int index, Vector2D p)
{
	cor[(size_t)index]=p;
}

int Delaunay::GetOriginalLength(void) const
{
	return olength;
}

vector<Vector2D>& Delaunay::GetMeshPoints(void)
{
	return cor;
}

int Delaunay::GetTotalLength(void)
{
	return (int)cor.size();
}

void Delaunay::AddBoundaryPoints(vector<Vector2D> const& points)
{
	int n=(int)points.size();
	//	vector<int> order=HilbertOrder(points,n);
	for(int i=0;i<n;++i)
	{
		cor.push_back(points[(size_t)i]);
		add_point((int)cor.size()-1);
	}
}

void Delaunay::AddAditionalPoint(Vector2D const& vec)
{
	cor.push_back(vec);
}

int Delaunay::GetCorSize(void)const
{
	return (int)cor.size();
}

bool Delaunay::IsTripleOut(int index) const
{
	int counter=0;
	for(size_t i=0;i<3;++i)
		if(IsOuterFacet(f[(size_t)index].neighbors[(size_t)i]))
			++counter;
	if(counter>1)
		return true;
	else
		return false;
}

int Delaunay::FindTripleLoc(facet const& fct)const
{
	for(size_t i=0;i<3;++i)
		if(!IsOuterFacet(fct.neighbors[(size_t)i]))
			return (i+1)%3;
	throw UniversalError("Trouble in constructing boundary triangles. No inner neighbor");
}

namespace
{
	bool IsOuterQuick(facet const& f,int olength)
	{
		for(size_t i=0;i<3;++i)
			if(f.vertices[(size_t)i]>=olength)
				return true;
		return false;
	}

	bool IsEdgeFacet(vector<facet> const& facets,facet const& f,int olength)
	{
		int counter=0;
		for(size_t i=0;i<3;++i)
		{
			if(f.vertices[(size_t)i]>=olength)
				return false;
			if(IsOuterQuick(facets[(size_t)f.neighbors[(size_t)i]],olength))
				++counter;
		}
		if(counter>0)
			return true;
		else
			return false;
	}

	bool CircleSegmentIntersect(Edge const& edge,Vector2D const& center,double R)
	{
		Vector2D AC=center-edge.vertices.first;
		Vector2D AB=edge.vertices.second-edge.vertices.first;
		double d=ScalarProd(AC,AB);
		if(d<0)
		{
			if(abs(AC)>R)
				return false;
			else
				return true;
		}
		double LAB=abs(AB);
		if(d>LAB*LAB)
		{
			if(abs(center-edge.vertices.second)>R)
				return false;
			else
				return true;
		}
		Vector2D closest=edge.vertices.first+AB*d/(LAB*LAB);
		if(abs(center-closest)>R)
			return false;
		else
			return true;
	}
}

vector<int> Delaunay::GetOuterFacets(int start_facet,int real_point,int olength2)
{
	int cur_facet=start_facet;
	vector<int> f_temp,containing_facets;
	f_temp.reserve((size_t)(10*sqrt(1.0*olength2)));
	int point_index=FindPointInFacet(cur_facet,real_point);
	if(IsOuterQuick(f[(size_t)f[(size_t)cur_facet].neighbors[(size_t)point_index]],olength2))
	{
		point_index=(point_index+1)%3;
		real_point=f[(size_t)cur_facet].vertices[(size_t)point_index];
	}
	if(IsOuterQuick(f[(size_t)f[(size_t)cur_facet].neighbors[(size_t)point_index]],olength2))
	{
		point_index=(point_index+1)%3;
		real_point=f[(size_t)cur_facet].vertices[(size_t)point_index];
	}
	do
	{
		FindContainingTetras(cur_facet,real_point,containing_facets);
		int old_current=cur_facet;
		for(size_t i=0;i<containing_facets.size();++i)
		{
			if(IsEdgeFacet(f,f[(size_t)containing_facets[(size_t)i]],olength2)&&
				containing_facets[(size_t)i]!=old_current)
				cur_facet=containing_facets[(size_t)i];
			if(!IsOuterQuick(f[(size_t)containing_facets[(size_t)i]],olength2))
				f_temp.push_back(containing_facets[(size_t)i]);
		}
		point_index=(1+FindPointInFacet(cur_facet,real_point))%3;
		if(IsTripleOut(cur_facet))
			point_index=(point_index+1)%3;
		real_point=f[(size_t)cur_facet].vertices[(size_t)point_index];
	}while(start_facet!=cur_facet);
	sort(f_temp.begin(),f_temp.end());
	f_temp=unique(f_temp);
	return f_temp;
}

#ifdef RICH_MPI
vector<vector<int> > Delaunay::FindOuterPointsMPI(OuterBoundary const* obc,
						  vector<Edge> const& edges,
						  Tessellation const& tproc,
						  vector<vector<int> > &Nghost,
						  vector<int> &proclist)
{
	// We add the points in a counter clockwise fashion
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	vector<vector<int> > res(edges.size());
	vector<int> res_temp, outer_points, f_temp, f_add(f.size(), 0);
	if (olength>100)
	{
		res_temp.reserve((size_t)(20 * sqrt(1.0*olength)));
		outer_points.reserve((size_t)(10 * sqrt(1.0*olength)));
		// Walk to an outer point
		int cur_facet = Walk(olength);
		// Find the real point
		int real_point = 0;
		for (int i = 0; i < 3; ++i)
		{
			if (f[(size_t)cur_facet].vertices[(size_t)i] < olength)
			{
				real_point = f[(size_t)cur_facet].vertices[(size_t)i];
				break;
			}
		}
		vector<int> containing_facets;
		FindContainingTetras(cur_facet, real_point, containing_facets);
		for (size_t i = 0; i < containing_facets.size(); ++i)
		{
			if (IsEdgeFacet(f, f[(size_t)containing_facets[(size_t)i]], olength))
			{
				cur_facet = containing_facets[(size_t)i];
				break;
			}
		}

		cur_facet = Walk(real_point);
		FindContainingTetras(cur_facet, real_point, containing_facets);
		for (size_t i = 0; i < containing_facets.size(); ++i)
		{
			if (IsEdgeFacet(f, f[(size_t)containing_facets[(size_t)i]], olength))
			{
				cur_facet = containing_facets[(size_t)i];
				break;
			}
		}
		f_temp = GetOuterFacets(cur_facet, real_point, olength);
	}
	else
	{
		for (size_t i = 0; i < f.size(); ++i)
		{
			if (!IsOuterQuick(f[i], olength))
				f_temp.push_back(i);
		}
	}
	// sort all the points and facets with the correct neighbors
	vector<int> allpoints;
	for(size_t i=0;i<f_temp.size();++i)
		for(size_t j=0;j<3;++j)
			allpoints.push_back(f[f_temp[i]].vertices[j]);
	sort(allpoints.begin(),allpoints.end());
	allpoints=unique(allpoints);
	allpoints=VectorValues(allpoints,HilbertOrder(VectorValues(cor,allpoints),(int)allpoints.size()));
	for(size_t i=0;i<allpoints.size();++i)
	{
		vector<int> neigh2=FindContainingTetras(Walk(allpoints[i]),allpoints[i]);
		for(size_t k=0;k<neigh2.size();++k)
		{
			const Vector2D center=GetCircleCenter(neigh2[k]);
			for(size_t j=0;j<edges.size();++j)
				if(CircleSegmentIntersect(edges[j],center,radius[neigh2[k]]))
					res[j].push_back(allpoints[i]);	
		}
	}
	for(size_t i=0;i<res.size();++i)
	{
		sort(res[i].begin(),res[i].end());
		res[i]=unique(res[i]);
		res[i]=VectorValues(res[i],HilbertOrder(VectorValues(cor,res[i]),(int)res[i].size()));
	}
	// Send/Recv the data
	vector<vector<Vector2D> > tosend;
	vector<vector<int> > ownsend,oldres;
	vector<Edge> ownedge;
	vector<int> neigh;
	for(size_t i=0;i<edges.size();++i)
	{
		const int temp=(edges[i].neighbors.first==rank) ? edges[i].neighbors.second : edges[i].neighbors.first;
		if(temp<0)
		{
			ownsend.push_back(res[i]);
			ownedge.push_back(edges[i]);
		}
		else
		{
			neigh.push_back(temp);
			tosend.push_back(VectorValues(cor,res[i]));
			oldres.push_back(res[i]);
		}
	}
	try
	{
		// Do own send
		AddRigid(obc, ownedge, ownsend);
		// Send/Recv from neighbors
		SendRecvFirstBatch(tosend, neigh, Nghost);
	}
	catch (UniversalError &eo)
	{
		eo.AddEntry("Error occured in first batch of MPI boundary creation", 0);
		throw;
	}
	vector<int> orgneigh = neigh;
	// Now find all other candidates
	// Recursively look for more points
	vector<bool> checked((size_t)olength,false);
	vector<vector<int> > toduplicate;
	toduplicate.resize(neigh.size());

	//cur_facet=Walk(real_point);
	//FindContainingTetras(cur_facet,real_point,containing_facets);
	for(size_t i=0;i<allpoints.size();++i)
		AddOuterFacetsMPI(allpoints[i],toduplicate,neigh,checked,tproc);
	for(size_t i=0;i<neigh.size();++i)
	{
		sort(toduplicate[i].begin(),toduplicate[i].end());
		toduplicate[i]=unique(toduplicate[i]);
	}
	// Add the points
	// Make sure not to send twice
	for(size_t i=0;i<oldres.size();++i)
	{
		vector<int>  vtemp(oldres[i]);
		sort(vtemp.begin(),vtemp.end());
		toduplicate[i]=RemoveList(toduplicate[i],vtemp);
	}
	// Make sure not to send rigid walls twice, and find if point intersect more than one wall
	vector<Edge> outeredges=obc->GetBoxEdges();
	vector<Edge> etemp=ownedge;
	// Rearrange the edges to agree with cell edges
	for(size_t i=0;i<outeredges.size();++i)
	{
		bool good=true;
		const double l1=outeredges[i].GetLength();
		for(size_t j=0;j<ownedge.size();++j)
		{
			const double l2=ownedge[j].GetLength();
			// are the edges parallel?
			const double parvalue=std::abs(ScalarProd(outeredges[i].vertices.first-outeredges[i].vertices.second,
				ownedge[j].vertices.first-ownedge[j].vertices.second));
			if(parvalue>0.99*l1*l2)
				if(DistanceToEdge(ownedge[j].vertices.first,outeredges[i])<l2*1e-8)
				{
					good=false;
					break;
				}
		}
		if(good)
			etemp.push_back(outeredges[i]);
	}
	outeredges=etemp;
	vector<vector<int> > extradd(outeredges.size());
	if(toduplicate.size()>oldres.size())
	{
		for(size_t i=oldres.size();i<toduplicate.size();++i)
		{
			if(neigh[i]==-1)
			{
				for(size_t j=0;j<outeredges.size();++j)
				{
					vector<int> tempvec;
					if(j<ownedge.size()&&!ownsend[j].empty())
					{
						tempvec=ownsend[j];
						sort(tempvec.begin(),tempvec.end());
					}
					for(size_t k=0;k<toduplicate[i].size();++k)
					{
						if(tempvec.empty()||(!std::binary_search(tempvec.begin(),tempvec.end(),toduplicate[i][k])))
						{
							if(CircleSegmentIntersect(outeredges[j],cor[toduplicate[i][k]],
								2*GetMaxRadius(toduplicate[i][k],Walk(toduplicate[i][k]))))
								extradd[j].push_back(toduplicate[i][k]);
						}
					}
				}
			}
		}
		// Remove the rigid points
		vector<vector<int> > dtemp;
		vector<int> ntemp;
		for(size_t i=0;i<toduplicate.size();++i)
		{
			if(neigh[i]>=0)
			{
				dtemp.push_back(toduplicate[i]);
				ntemp.push_back(neigh[i]);
			}
		}
		neigh=ntemp;
		toduplicate=dtemp;
	}

	sort(orgneigh.begin(), orgneigh.end());
	// Make sure we have that all cpu talk with the relevent other and remove uneeded talks
	vector<int> sendnumber(get_mpi_size(),0),scounts(get_mpi_size(),1);
	int nrecv;
	const double R = tproc.GetWidth(rank);
	vector<int> neightemp2;
	for (size_t i = 0; i < neigh.size(); ++i)
	{
		if (tproc.GetMeshPoint(neigh[i]).distance(tproc.GetMeshPoint(rank)) <( 15 * R) || std::binary_search(orgneigh.begin(), orgneigh.end(), neigh[i]))
		{
			sendnumber[neigh[i]] = 1;
			neightemp2.push_back(neigh[i]);
		}
	}
	sort(neightemp2.begin(), neightemp2.end());
	MPI_Reduce_scatter(&sendnumber[0],&nrecv,&scounts[0],MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	char ctemp;
	vector<MPI_Request> req(neightemp2.size());
	for (size_t i = 0; i<neightemp2.size(); ++i)
		MPI_Isend(&ctemp, 1, MPI_CHAR, neightemp2[i], 2, MPI_COMM_WORLD, &req[i]);
	vector<int> sentme;
	for(size_t i=0;i<(size_t)nrecv;++i)
	{
		MPI_Status status;
		MPI_Probe(MPI_ANY_SOURCE,2,MPI_COMM_WORLD,&status);
		MPI_Recv(&ctemp,1,MPI_CHAR,status.MPI_SOURCE,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		sentme.push_back(status.MPI_SOURCE);
	}
	sort(sentme.begin(),sentme.end());
	MPI_Barrier(MPI_COMM_WORLD);
	vector<vector<int> > toduptemp;
	sendnumber.clear();
	for(size_t i=0;i<neigh.size();++i)
	{
		if (std::binary_search(sentme.begin(), sentme.end(), neigh[i]) && std::binary_search(neightemp2.begin(), neightemp2.end(), neigh[i]))
		{
			toduptemp.push_back(toduplicate[i]);
			sendnumber.push_back(neigh[i]);
		}
	}
	vector<int> oldneigh = neigh;
	neigh=sendnumber;
	toduplicate=toduptemp;
	vector<int> neightemp;
	tosend.clear();
	for(size_t i=0;i<neigh.size();++i)
	{
		if(neigh[i]>=0)
		{
			neightemp.push_back(neigh[i]);
			tosend.push_back(VectorValues(cor,toduplicate[i]));
		}
	}
	neigh=neightemp;

	try
	{
		SendRecvFirstBatch(tosend, neigh, Nghost);
	}
		catch (UniversalError &eo)
		{
			eo.AddEntry("Error occured in second batch of MPI boundary creation", 0);
			for (size_t i = 0; i < neigh.size(); ++i)
				eo.AddEntry("Talked to cpu=", neigh[i]);
			for (size_t i = 0; i < oldneigh.size(); ++i)
				eo.AddEntry("Wanted to talk with cpu=", oldneigh[i]);
			delaunay_loggers::BinaryLogger log("del.bin");
			log.output(cor, f);
			throw;
	}
	
	for(size_t i=0;i<oldres.size();++i)
		if(!oldres[i].empty())
			toduplicate[i].insert(toduplicate[i].begin(),oldres[i].begin(),oldres[i].end());
	proclist=neigh;

	// Add the second batch of rigid points, first make sure we need all of the points, we ddo not want points that are no longer needed
	try
	{
		vector<vector<int> > newadd(outeredges.size());
		for (size_t i = 0; i < extradd.size(); ++i)
		{
			for (size_t j = 0; j < extradd[i].size(); ++j)
			{
				if (CircleSegmentIntersect(outeredges[i], cor[extradd[i][j]],
					2*GetMaxRadius(extradd[i][j], Walk(extradd[i][j]))))
					newadd[i].push_back(extradd[i][j]);
			}
		}
		AddRigid(obc, outeredges, newadd);
	}
	catch (UniversalError &eo)
	{
		eo.AddEntry("Error occured in second batch of MPI boundary creation of rigid points", 0);
		for (size_t i = 0; i < neigh.size(); ++i)
			eo.AddEntry("Talked to cpu=", neigh[i]);
		for (size_t i = 0; i < oldneigh.size(); ++i)
			eo.AddEntry("Wanted to talk with cpu=", oldneigh[i]);
		for (size_t i = 0; i < outeredges.size(); ++i)
		{
			eo.AddEntry("Edge x0=", outeredges[i].vertices.first.x);
			eo.AddEntry("Edge y0=", outeredges[i].vertices.first.y);
			eo.AddEntry("Edge x1=", outeredges[i].vertices.second.x);
			eo.AddEntry("Edge y1=", outeredges[i].vertices.second.y);
			eo.AddEntry("Edge n0=", outeredges[i].neighbors.first);
			eo.AddEntry("Edge n0=", outeredges[i].neighbors.second);
			for (size_t j = 0; j < extradd[i].size(); ++j)
				eo.AddEntry("Added point", extradd[i][j]);
		}
		delaunay_loggers::BinaryLogger log("del.bin");
		log.output(cor, f);
		throw;
	}
	return toduplicate;
}
#endif

#ifdef RICH_MPI
void Delaunay::SendRecvFirstBatch(vector<vector<Vector2D> > &tosend,
	vector<int> const& proclist,vector<vector<int> > &Nghost)
{
	const int rank = get_mpi_rank();
	const int ws = get_mpi_size();
	vector<int> procorder=GetProcOrder(rank,ws);
	int nlist=(int)proclist.size();
	Nghost.resize(proclist.size());
	for(size_t i=0;i<procorder.size();++i)
	{
		const int index=find(proclist.begin(),proclist.end(),procorder[i])-proclist.begin();
		// Do we talk with this processor?
		if(index<nlist)
		{
			// Create send data
			vector<double> send;
			ConvertVector2DToDouble(tosend[index],send);
			// Send/Recv data
			MPI_Status status;
			vector<double> recv;
			int nrecv;
			if(rank<procorder[(size_t)i])
			{
				if(!send.empty())
					MPI_Send(&send[0],(int)send.size(),MPI_DOUBLE,procorder[i],0,MPI_COMM_WORLD);
				else
				{
					double temp=0;
					MPI_Send(&temp,1,MPI_DOUBLE,procorder[i],1,MPI_COMM_WORLD);
				}
				MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
				MPI_Get_count(&status,MPI_DOUBLE,&nrecv);
				recv.resize(nrecv);
				int rtag=status.MPI_TAG;
				if(rtag==0)
					MPI_Recv(&recv[0],nrecv,MPI_DOUBLE,procorder[i],0,MPI_COMM_WORLD,&status);
				else
				{
					double temp=0;
					MPI_Recv(&temp,1,MPI_DOUBLE,procorder[i],1,MPI_COMM_WORLD,&status);
				}
			}
			else
			{
				MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
				MPI_Get_count(&status,MPI_DOUBLE,&nrecv);
				recv.resize(nrecv);
				if(status.MPI_TAG==0)
					MPI_Recv(&recv[0],nrecv,MPI_DOUBLE,procorder[i],0,MPI_COMM_WORLD,&status);
				else
				{
					double temp=0;
					MPI_Recv(&temp,1,MPI_DOUBLE,procorder[i],1,MPI_COMM_WORLD,&status);
				}
				if(!send.empty())
					MPI_Send(&send[0],(int)send.size(),MPI_DOUBLE,procorder[i],0,MPI_COMM_WORLD);
				else
				{
					double temp=0;
					MPI_Send(&temp,1,MPI_DOUBLE,procorder[i],1,MPI_COMM_WORLD);
				}
			}
			vector<Vector2D> toadd;
			ConvertDoubleToVector2D(toadd,recv);
			// Add the points
			try
			{

				if (!toadd.empty())
				{
					for (size_t i = 0; i < toadd.size(); ++i)
						Nghost[index].push_back(cor.size() + i);
					AddBoundaryPoints(toadd);
				}
			}
			catch (UniversalError &eo)
			{
				eo.AddEntry("Error in SendRecvFirstBatch", 0);
				throw;
			}
		}
	}
}
#endif

vector<vector<int> > Delaunay::FindOuterPoints(vector<Edge> const& edges)
{
	// We add the points in a counter clockwise fashion
	vector<vector<int> > res(edges.size());
	if(olength<100)
	{
		for(size_t j=0;j<edges.size();++j)
		{
			res[(size_t)j].resize((size_t)olength);
			for(size_t i=0;i<olength;++i)
			  res[(size_t)j][i]=(int)i;
		}
		return res;
	}
	vector<int> res_temp,outer_points,f_temp,f_add(f.size(),0);
	res_temp.reserve((size_t)(20*sqrt(1.0*olength)));
	f_temp.reserve((size_t)(10*sqrt(1.0*olength)));
	outer_points.reserve((size_t)(10*sqrt(1.0*olength)));
	// Walk to an outer point
	int cur_facet=Walk(olength);
	// Find the real point
	int real_point = 0;
	for(int i=0;i<3;++i)
	{
	  if(f[(size_t)cur_facet].vertices[(size_t)i]<(int)olength)
		{
			real_point=f[(size_t)cur_facet].vertices[(size_t)i];
			break;
		}
	}
	vector<int> containing_facets;
	FindContainingTetras(cur_facet,real_point,containing_facets);
	for(size_t i=0;i<containing_facets.size();++i)
	{
		if(IsEdgeFacet(f,f[(size_t)containing_facets[(size_t)i]],olength))
		{
			cur_facet=containing_facets[(size_t)i];
			break;
		}
	}
	int start_facet=cur_facet;
	int point_index=FindPointInFacet(cur_facet,real_point);
	// Make sure we are at the right location in the triangle, we want to be in the last outer point
	if(IsOuterQuick(f[(size_t)f[(size_t)cur_facet].neighbors[(size_t)point_index]],olength))
	{
		point_index=(point_index+1)%3;
		real_point=f[(size_t)cur_facet].vertices[(size_t)point_index];
	}
	if(IsOuterQuick(f[(size_t)f[(size_t)cur_facet].neighbors[(size_t)point_index]],olength))
	{
		point_index=(point_index+1)%3;
		real_point=f[(size_t)cur_facet].vertices[(size_t)point_index];
	}
	do
	{
		FindContainingTetras(cur_facet,real_point,containing_facets);
		int old_current=cur_facet;
		for(size_t i=0;i<containing_facets.size();++i)
		{
			if(IsEdgeFacet(f,f[(size_t)containing_facets[(size_t)i]],olength)&&
				containing_facets[(size_t)i]!=old_current)
				cur_facet=containing_facets[(size_t)i];
			if(!IsOuterQuick(f[(size_t)containing_facets[(size_t)i]],olength))
				f_temp.push_back(containing_facets[(size_t)i]);
		}
		point_index=(1+FindPointInFacet(cur_facet,real_point))%3;
		if(IsTripleOut(cur_facet))
			point_index=(point_index+1)%3;
		real_point=f[(size_t)cur_facet].vertices[(size_t)point_index];
	}while(start_facet!=cur_facet);
	sort(f_temp.begin(),f_temp.end());
	f_temp=unique(f_temp);
	// Found all initial outer facets

	//Find the points in the outer facets
	vector<vector<int> > toduplicate(edges.size());
	// Recursively look for more points
	vector<bool> checked((size_t)olength,false);
	for(size_t i=0;i<f_temp.size();++i)
		AddOuterFacets(f_temp[i],toduplicate,edges,checked);
	for(size_t i=0;i<edges.size();++i)
	{
		sort(toduplicate[(size_t)i].begin(),toduplicate[(size_t)i].end());
		toduplicate[(size_t)i]=unique(toduplicate[(size_t)i]);
	}
	return toduplicate;
}

void Delaunay::AddRigid(OuterBoundary const* /*obc*/,vector<Edge> const& edges,
	vector<vector<int> > &toduplicate)
{
	vector<int> toremove;
	for(size_t i=0;i<edges.size();++i)
	{
		toremove.clear();
		if(toduplicate[i].empty())
			continue;
		vector<Vector2D> toadd;
		toadd.reserve(toduplicate[i].size());
		Vector2D par(Parallel(edges[i]));
		par=par/abs(par);
		for(size_t j=0;j<toduplicate[i].size();++j)
		{
			Vector2D temp=cor[(size_t)toduplicate[i][j]]-edges[i].vertices.first;
			temp=2*par*ScalarProd(par,temp)-temp+edges[i].vertices.first;
			if (InTriangle(TripleConstRef<Vector2D>
				       (cor[static_cast<size_t>(olength)],
					cor[static_cast<size_t>(olength+1)],
					cor[static_cast<size_t>(olength+2)]), 
				       temp))
				toadd.push_back(temp);
			else
				toremove.push_back((int)j);
		}
		RemoveVector(toduplicate[i], toremove);
		vector<int> order=HilbertOrder(toadd,(int)toadd.size());
		ReArrangeVector(toadd,order);
		try
		{
			AddBoundaryPoints(toadd);
		}
		catch (UniversalError &eo)
		{
			eo.AddEntry("Error in AddRigid", 0);
			throw;
		}
		ReArrangeVector(toduplicate[i],order);
	}
}

namespace
{
	/*
	vector<vector<int> > FindCorners(vector<vector<int> > const& toduplicate)
	{
	vector<vector<int> > res(toduplicate.size());
	for(size_t i=0;i<toduplicate.size();++i)
	{
	for(size_t j=0;j<toduplicate[(size_t)i].size();++j)
	{
	if(binary_search(toduplicate[(size_t)(i+1)%4].begin(),toduplicate[(size_t)(i+1)%4].end(),
	toduplicate[(size_t)i][(size_t)j]))
	res[(size_t)i].push_back(toduplicate[(size_t)i][(size_t)j]);
	}
	}
	return res;
	}
	*/
}
namespace
{
	vector<Edge> GetCornerEdges(OuterBoundary const* obc)
	{
		const double dx=obc->GetGridBoundary(Right)-obc->GetGridBoundary(Left);
		const double dy=obc->GetGridBoundary(Up)-obc->GetGridBoundary(Down);
		vector<Edge> res;
		const Vector2D RU(obc->GetGridBoundary(Right),obc->GetGridBoundary(Up));
		const Vector2D LU(obc->GetGridBoundary(Left),obc->GetGridBoundary(Up));
		const Vector2D LD(obc->GetGridBoundary(Left),obc->GetGridBoundary(Down));
		const Vector2D RD(obc->GetGridBoundary(Right),obc->GetGridBoundary(Down));
		res.push_back(Edge(RU,Vector2D(dx,0)+RU,0,0));
		res.push_back(Edge(RU,Vector2D(0,dy)+RU,0,0));
		res.push_back(Edge(LU,Vector2D(0,dy)+LU,0,0));
		res.push_back(Edge(LU,Vector2D(-dx,0)+LU,0,0));
		res.push_back(Edge(LD,Vector2D(-dx,0)+LD,0,0));
		res.push_back(Edge(LD,Vector2D(0,-dy)+LD,0,0));
		res.push_back(Edge(RD,Vector2D(0,-dy)+RD,0,0));
		res.push_back(Edge(RD,Vector2D(dx,0)+RD,0,0));
		return res;
	}
}

vector<vector<int> > Delaunay::AddPeriodic(OuterBoundary const* obc,vector<Edge> const& edges,
	vector<vector<int> > &toduplicate)
{
	const double dx=obc->GetGridBoundary(Right)-obc->GetGridBoundary(Left);
	const double dy=obc->GetGridBoundary(Up)-obc->GetGridBoundary(Down);
	for(size_t i=0;i<edges.size();++i)
	{
		if(toduplicate[(size_t)i].empty())
			continue;
		Vector2D change;
		switch(i)
		{
		case(0):
			change.x=-dx;
			break;
		case(1):
			change.y=-dy;
			break;
		case(2):
			change.x=dx;
			break;
		case(3):
			change.y=dy;
			break;
		}
		vector<Vector2D> toadd;
		toadd.reserve(toduplicate[(size_t)i].size());
		//vector<int> pointstemp(toduplicate[(size_t)i].size());
		for(size_t j=0;j<toduplicate[(size_t)i].size();++j)
		{
			toadd.push_back(cor[(size_t)toduplicate[(size_t)i][(size_t)j]]+change);
			//	pointstemp[j]=j;
		}
		vector<int> order=HilbertOrder(toadd,(int)toadd.size());
		ReArrangeVector(toadd,order);
		AddBoundaryPoints(toadd);
		ReArrangeVector(toduplicate[(size_t)i],order);
		//toduplicate[i]=pointstemp;
	}
	// Done with sides do corners now
	vector<Edge> corneredges=GetCornerEdges(obc);
	vector<vector<int> > corners(toduplicate.size());
	for(size_t i=0;i<toduplicate.size();++i)
	{
		for(size_t j=0;j<toduplicate[(size_t)i].size();++j)
		{
			const int facet_loc=Walk(toduplicate[(size_t)i][(size_t)j]);
			const Vector2D center=cor[(size_t)toduplicate[(size_t)i][(size_t)j]];
			const double R=2*GetMaxRadius(toduplicate[(size_t)i][(size_t)j],facet_loc);
			if(CircleSegmentIntersect(corneredges[(size_t)2*i],center,R))
				corners[(size_t)i].push_back(toduplicate[(size_t)i][(size_t)j]);
			if(CircleSegmentIntersect(corneredges[(size_t)(2*i+7)%8],center,R))
				corners[(size_t)(i+3)%4].push_back(toduplicate[(size_t)i][(size_t)j]);
		}
	}
	for(size_t i=0;i<corners.size();++i)
	{
		if(corners[(size_t)i].empty())
			continue;
		sort(corners[(size_t)i].begin(),corners[(size_t)i].end());
		corners[(size_t)i]=unique(corners[(size_t)i]);
		Vector2D change;
		switch(i)
		{
		case(0):
			change.x=-dx;
			change.y=-dy;
			break;
		case(1):
			change.y=-dy;
			change.x=dx;
			break;
		case(2):
			change.x=dx;
			change.y=dy;
			break;
		case(3):
			change.y=dy;
			change.x=-dx;
			break;
		}
		vector<Vector2D> toadd;
		toadd.reserve(corners[(size_t)i].size());
		//		vector<int> pointstemp(corners[i].size());
		for(size_t j=0;j<corners[(size_t)i].size();++j)
		{
			toadd.push_back(cor[(size_t)corners[(size_t)i][(size_t)j]]+change);
			//		pointstemp[j]=j;
		}
		vector<int> order=HilbertOrder(toadd,(int)toadd.size());
		ReArrangeVector(toadd,order);
		AddBoundaryPoints(toadd);
		ReArrangeVector(corners[(size_t)i],order);
		//	corners[i]=pointstemp;
	}
	return corners;
}

void Delaunay::AddHalfPeriodic(OuterBoundary const* obc,vector<Edge> const& edges,
	vector<vector<int> > &toduplicate)
{
	const double dx=obc->GetGridBoundary(Right)-obc->GetGridBoundary(Left);
	//	const double dy=obc->GetGridBoundary(Up)-obc->GetGridBoundary(Down);
	for(size_t i=0;i<edges.size();++i)
	{
		if(toduplicate[(size_t)i].empty())
			continue;
		Vector2D change;
		switch(i)
		{
		case(0):
			change.x=-dx;
			break;
		case(1):
			break;
		case(2):
			change.x=dx;
			break;
		case(3):
			break;
		}
		vector<Vector2D> toadd;
		toadd.reserve(toduplicate[(size_t)i].size());
		//vector<int> pointstemp(toduplicate[(size_t)i].size());
		Vector2D par(Parallel(edges[(size_t)i]));
		par=par/abs(par);
		for(size_t j=0;j<toduplicate[(size_t)i].size();++j)
		{
			Vector2D temp=cor[(size_t)toduplicate[(size_t)i][(size_t)j]];
			if(i%2==1)
			{
				temp-=edges[(size_t)i].vertices.first;
				temp=2*par*ScalarProd(par,temp)-temp+edges[(size_t)i].vertices.first;
			}
			toadd.push_back(temp+change);
			//pointstemp[j]=j;
		}
		vector<int> order=HilbertOrder(toadd,(int)toadd.size());
		ReArrangeVector(toadd,order);
		AddBoundaryPoints(toadd);
		ReArrangeVector(toduplicate[(size_t)i],order);
		//toduplicate[i]=pointstemp;
	}
}

#ifdef RICH_MPI
vector<vector<int> > Delaunay::BuildBoundary(OuterBoundary const* obc,
	Tessellation const& tproc,vector<vector<int> > &Nghost,vector<int> &proclist)
{
	vector<Edge> edges;
	vector<int> edge_index=tproc.GetCellEdges(get_mpi_rank());
	for(size_t i=0;i<edge_index.size();++i)
		edges.push_back(tproc.GetEdge(edge_index[i]));
//	delaunay_loggers::BinaryLogger log("del"+int2str(get_mpi_rank())+".bin");
//	log.output(cor,f);
	return FindOuterPointsMPI(obc,edges,tproc,Nghost,proclist);
}
#endif

vector<vector<int> > Delaunay::BuildBoundary(OuterBoundary const* obc,vector<Edge> const& edges)
{
//	delaunay_loggers::BinaryLogger log("c:\\del.bin");
//	log.output(cor,f);
	vector<vector<int> > toduplicate=FindOuterPoints(edges);
	if(obc->GetBoundaryType()==Rectengular)
	{
		AddRigid(obc,edges,toduplicate);
	}
	else
	{
		if(obc->GetBoundaryType()==Periodic)
		{
			vector<vector<int> > corners=AddPeriodic(obc,edges,toduplicate);
			for(size_t i=0;i<4;++i)
				toduplicate.push_back(corners[(size_t)i]);
		}
		else
		{
			AddHalfPeriodic(obc,edges,toduplicate);
		}
	}
	return toduplicate;
}

Vector2D Delaunay::GetCircleCenter(int index)const
{
	Vector2D center;
	facet const& F=f[(size_t)index];
	double x1=cor[(size_t)F.vertices[0]].x;
	double x2=cor[(size_t)F.vertices[1]].x;
	double x3=cor[(size_t)F.vertices[2]].x;
	double y1=cor[(size_t)F.vertices[0]].y;
	double y2=cor[(size_t)F.vertices[1]].y;
	double y3=cor[(size_t)F.vertices[2]].y;
	// Do we have a case where two point are very close compared to the third?
	double d12=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
	double d23=(x3-x2)*(x3-x2)+(y3-y2)*(y3-y2);
	double d13=(x1-x3)*(x1-x3)+(y1-y3)*(y1-y3);
	int scenario=0;
	if(d12<0.1*(d23+d13))
		scenario=1;
	else
		if(d23<0.1*(d13+d12))
			scenario=3;
		else
			if(d13<0.1*(d23+d12))
				scenario=2;
	switch(scenario)
	{
	case(0):
	case(1):
	case(2):
		{
			x2-=x1;
			x3-=x1;
			y2-=y1;
			y3-=y1;
			double d_inv=1/(2*(x2*y3-y2*x3));
			center.Set((y3*(x2*x2+y2*y2)-y2*(x3*x3+y3*y3))*d_inv+x1,
				(-x3*(x2*x2+y2*y2)+x2*(x3*x3+y3*y3))*d_inv+y1);
			break;
		}
	case(3):
		{
			x1-=x2;
			x3-=x2;
			y1-=y2;
			y3-=y2;
			double d_inv=1/(2*(x3*y1-y3*x1));
			center.Set((y1*(x3*x3+y3*y3)-y3*(x1*x1+y1*y1))*d_inv+x2,
				(x3*(x1*x1+y1*y1)-x1*(x3*x3+y3*y3))*d_inv+y2);
			break;
		}
	default:
		throw UniversalError("Unhandled case in switch statement VoronoiMesh::get_center");
	}
	return center;
}

double Delaunay::GetMaxRadius(int point,int startfacet)
{
	double res=0;
	vector<int> neigh=FindContainingTetras(startfacet,point);
	for(size_t i=0;i<neigh.size();++i)
		res=max(res,radius[(size_t)neigh[(size_t)i]]);
	return res;
}

void Delaunay::AddOuterFacets(int tri,vector<vector<int> > &toduplicate,
	vector<Edge> const& edges,vector<bool> &checked)
{
	stack<int> tocheck;
	tocheck.push(tri);
	while(!tocheck.empty())
	{
		int cur_facet=tocheck.top();
		tocheck.pop();
		for(size_t i=0;i<3;++i)
		{
			bool added=false;
			if(checked[(size_t)f[(size_t)cur_facet].vertices[(size_t)i]]||(f[(size_t)cur_facet].vertices[(size_t)i]>=(int)olength))
				continue;
			vector<int> neigh=FindContainingTetras(cur_facet,f[(size_t)cur_facet].vertices[(size_t)i]);
			for(size_t k=0;k<neigh.size();++k)
			{
				Vector2D center=GetCircleCenter(neigh[(size_t)k]);
				for(size_t l=0;l<edges.size();++l)
				{
					if(CircleSegmentIntersect(edges[(size_t)l],center,radius[(size_t)neigh[(size_t)k]]))
					{
						toduplicate[(size_t)l].push_back(f[(size_t)cur_facet].vertices[(size_t)i]);
						added=true;
					}
				}
			}
			checked[(size_t)f[(size_t)cur_facet].vertices[(size_t)i]]=true;
			if(added)
			{
				for(size_t j=0;j<neigh.size();++j)
				{
					if(!IsOuterQuick(f[(size_t)neigh[(size_t)j]],olength))
						tocheck.push(neigh[(size_t)j]);
				}
			}
		}
	}
}

#ifdef RICH_MPI
void Delaunay::AddOuterFacetsMPI(int point,vector<vector<int> > &toduplicate,
	vector<int> &neigh,vector<bool> &checked,Tessellation const &tproc)
{
	stack<int> tocheck;
	vector<int> neightemp=FindContainingTetras(Walk(point),point);
	for(size_t i=0;i<neightemp.size();++i)
		tocheck.push(neightemp[i]);
	const int rank=get_mpi_rank();
	while(!tocheck.empty())
	{
		int cur_facet=tocheck.top();
		tocheck.pop();
		for(size_t i=0;i<3;++i)
		{
			bool added=false;
			if(f[(size_t)cur_facet].vertices[i]>=olength)
				continue;
			if(checked[(size_t)f[(size_t)cur_facet].vertices[i]])
				continue;
			vector<int> neighs=FindContainingTetras(cur_facet,f[(size_t)cur_facet].vertices[(size_t)i]);
			for(size_t k=0;k<neighs.size();++k)
			{
				Circle circ(GetCircleCenter(neighs[k]),radius[(size_t)neighs[k]]);
				vector<int> cputosendto;
				find_affected_cells(tproc,rank,circ,cputosendto);
				sort(cputosendto.begin(),cputosendto.end());
				cputosendto=unique(cputosendto);
				//				vector<int> toremove;
				RemoveVal(cputosendto,rank);
				if(!cputosendto.empty())
				{
					added=true;
					for(size_t j=0;j<cputosendto.size();++j)
					{
						size_t index=find(neigh.begin(),neigh.end(),cputosendto[j])
							-neigh.begin();
						if(index<neigh.size())
							toduplicate[index].push_back(f[(size_t)cur_facet].vertices[i]);
						else
						{
							neigh.push_back(cputosendto[j]);
							vector<int> tempvec;
							tempvec.push_back(f[(size_t)cur_facet].vertices[i]);
							toduplicate.push_back(tempvec);
						}
					}
				}
			}
			checked[(size_t)f[(size_t)cur_facet].vertices[i]]=true;
			if(added)
			{
				for(size_t j=0;j<neighs.size();++j)
				{
					if(!IsOuterQuick(f[(size_t)neighs[j]],olength))
						tocheck.push(neighs[j]);
				}
			}
		}
	}
}
#endif
