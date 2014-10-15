#include "Delaunay.hpp"
#include <vector>
#include <cmath>

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
	bool InTriangle(boost::array<Vector2D,3> const& tri,Vector2D const& point)
	{
		boost::array<Vector2D,3> tocheck;
		tocheck[2]=point;
		for(int i=0;i<3;++i)
		{
			tocheck[0]=tri[i];
			tocheck[1]=tri[(i+1)%3];
			if(orient2d(tocheck)<0)
				return false;
		}
		return true;
	}

	// Assume cell is orederd in convexhull counterclockwise
	bool InCell(vector<Vector2D> const& points,Vector2D const& p)
	{
		int n=(int)points.size();
		for(int i=0;i<n;++i)
		{
			if(CrossProduct(points[i]-p,points[(i+1)%n]-p)<0)
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
			minx=min(points[i].x,minx);
			miny=min(points[i].y,miny);
			maxx=max(points[i].x,maxx);
			maxy=max(points[i].y,maxy);
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
			if(fc.neighbors[j]==i)
				return j;
		}
		throw UniversalError("Error in find_index: Index not found");
	}
}

void Delaunay::add_point(int index)
{
	// Check if point is inside big triangle
	boost::array<Vector2D,3> tocheck;
	tocheck[0]=cor[olength];
	tocheck[1]=cor[olength+1];
	tocheck[2]=cor[olength+2];
	if(!InTriangle(tocheck,cor[index]))
	{
		UniversalError eo("Point not inside large triangle of Delaunay");
		eo.AddEntry("Point x",cor[index].x);
		eo.AddEntry("Point y",cor[index].y);
		eo.AddEntry("Big tirangle point 1 x",tocheck[0].x);
		eo.AddEntry("Big tirangle point 1 y",tocheck[0].y);
		eo.AddEntry("Big tirangle point 2 x",tocheck[1].x);
		eo.AddEntry("Big tirangle point 2 y",tocheck[1].y);
		eo.AddEntry("Big tirangle point 3 x",tocheck[2].x);
		eo.AddEntry("Big tirangle point 3 y",tocheck[2].y);
		throw eo;
	}
	int triangle=Walk(index);
	boost::array<int,3> outer,temp_friends;
	facet f_temp;
	int i=0;
	for(i=0;i<3;++i)
	{
		outer[i]=f[triangle].vertices[i];
		temp_friends[i]=f[triangle].neighbors[i];
	}
	// create and _update the new facets
	f.push_back(f_temp);
	f.push_back(f_temp);
	f[triangle].vertices[0] = outer[2];
	f[triangle].vertices[1] = outer[0];
	f[triangle].vertices[2] = index;
	f[location_pointer+1].vertices[0] = outer[0];
	f[location_pointer+1].vertices[1] = outer[1];
	f[location_pointer+1].vertices[2] = index;
	f[location_pointer+2].vertices[0] = outer[1];
	f[location_pointer+2].vertices[1] = outer[2];
	f[location_pointer+2].vertices[2] = index;

	f[triangle].neighbors[0] = temp_friends[2];
	f[triangle].neighbors[1] = location_pointer+1;;
	f[triangle].neighbors[2] = location_pointer+2;
	f[location_pointer+1].neighbors[0] = temp_friends[0];
	f[location_pointer+1].neighbors[1] = location_pointer+2;
	f[location_pointer+1].neighbors[2] = triangle;
	f[location_pointer+2].neighbors[0] = temp_friends[1];
	f[location_pointer+2].neighbors[1] = triangle;
	f[location_pointer+2].neighbors[2] = location_pointer+1;
	// _update the friends list of the friends
	if(temp_friends[1]!=last_loc)
	{
		i=find_index(f[temp_friends[1]],triangle);
		f[temp_friends[1]].neighbors[i] = location_pointer+2;
	}
	if(temp_friends[0]!=last_loc)
	{
		i=find_index(f[temp_friends[0]],triangle);
		f[temp_friends[0]].neighbors[i] = location_pointer+1;
	}
	// Calculate radius if needed
	if(CalcRadius)
	{
		radius[triangle]=CalculateRadius(triangle);
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
				radius[location_pointer+1]=CalculateRadius(location_pointer+1);
				radius.push_back(CalculateRadius(location_pointer+2));
			}
			else
			{
				radius[location_pointer+1]=CalculateRadius(location_pointer+1);
				radius[location_pointer+2]=CalculateRadius(location_pointer+2);
			}
	}

	// check if flipping is needed
	flip(triangle,temp_friends[2]);
	flip(location_pointer+1,temp_friends[0]);
	flip(location_pointer+2,temp_friends[1]);

	// _update number of facets
	location_pointer+=2;
}

void Delaunay::flip(int i, int j)
{
	stack<boost::array<int,2> > flip_stack;
	boost::array<int,2> array_temp={{i,j}};
	flip_stack.push(array_temp);
	boost::array<int,2> check;
	boost::array<int,2> other;
	boost::array<int,2> indexes;
	boost::array<Vector2D,4> circle_test;
	while(!flip_stack.empty())
	{
		if(flip_stack.top()[1]==last_loc)
			flip_stack.pop();
		else
		{
			indexes[0]=flip_stack.top()[0];
			indexes[1]=flip_stack.top()[1];
			// Returns the index to the point to check in coordinates and the index of the point in the facet
			find_diff(&f[indexes[1]],&f[indexes[0]],&check[0]);
			find_diff(&f[indexes[0]],&f[indexes[1]],&other[0]);

			for(int k=0;k<3;++k)
			{
				circle_test[k]=cor[f[indexes[0]].vertices[k]];
			}
			circle_test[3]=cor[check[0]];

			if(incircle(circle_test)>0)
			{
				//The point is in a circle change the facets and their friends
				const int v1=f[indexes[0]].vertices[(other[1]+1)%3];
				const int f1=f[indexes[0]].neighbors[other[1]];
				const int f12=f[indexes[0]].neighbors[(other[1]+2)%3];
				const int f13=f[indexes[0]].neighbors[(other[1]+2)%3];
				const int v2=f[indexes[1]].vertices[(check[1]+1)%3];
				const int f2=f[indexes[1]].neighbors[(check[1]+2)%3];
				const int f22=f[indexes[1]].neighbors[check[1]];
				const int f23=f[indexes[1]].neighbors[(check[1]+2)%3];
				f[indexes[0]].vertices[0] = other[0];
				f[indexes[0]].vertices[1] = v1;
				f[indexes[0]].vertices[2] = check[0];
				f[indexes[1]].vertices[0] = check[0];
				f[indexes[1]].vertices[1] = v2;
				f[indexes[1]].vertices[2] = other[0];
				f[indexes[0]].neighbors[0] = f1;
				f[indexes[0]].neighbors[1] = f2;
				f[indexes[0]].neighbors[2] = indexes[1];
				f[indexes[1]].neighbors[0] = f22;
				f[indexes[1]].neighbors[1] = f12;
				f[indexes[1]].neighbors[2] = indexes[0];
				// change the friends of the friends if needed
				if(f23!=last_loc)
				{
					f[f23].neighbors[find_index(f[f23],indexes[1])] = indexes[0];
				}
				if(f13!=last_loc)
				{
					f[f13].neighbors[find_index(f[f13],indexes[0])] = indexes[1];
				}
				// Calculate the new radius if needed
				if(CalcRadius)
				{
					radius[indexes[0]]=CalculateRadius(indexes[0]);
					radius[indexes[1]]=CalculateRadius(indexes[1]);
				}
				// clear the checked facets
				flip_stack.pop();
				// push into the stack the new facets to check
				array_temp[0]=indexes[1];
				array_temp[1]=f[indexes[1]].neighbors[0];
				flip_stack.push(array_temp);
				array_temp[0]=indexes[0];
				array_temp[1]=f[indexes[0]].neighbors[1];
				flip_stack.push(array_temp);
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
	f.reserve(2*length+1+(int)(17*sqrt(1.*length)));
	cor.reserve(length+9*(int)sqrt(1.*length));
	last_loc=INT_MAX;
	for(int i=0;i<len;i++)
	{
		cor.push_back(vp[i]);
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
	for(int i=0;i<3;i++)
		f[0].neighbors[i] = last_loc;
	location_pointer=0;
	// add the points
	for(int i=0;i<length-3;i++)
	{
		add_point(data.insert_order[i]);
	}
	// Calculate radius
	radius.resize(f.size());
	int n=int(f.size());
	for(int i=0;i<n;++i)
		radius[i]=CalculateRadius(i);
	CalcRadius=true;
}

void Delaunay::find_diff(facet *f1,facet *f2,int *p) const
{
	int i=0;
	for(;i<3;++i)
	{
		bool counter=false;
		for(int j=0;j<3;++j)
		{
			if(f1->vertices[i]==f2->vertices[j])
			{
				counter=true;
				break;
			}
		}
		if(counter==false)
		{
			p[0]=f1->vertices[i];
			p[1]=i;
			return;
		}
	}
	UniversalError eo("Delaunay, Couldn't find difference bewteen two facets");
	eo.AddEntry("Facet 1",double(f1-&f[0]));
	eo.AddEntry("Facet 2",double(f2-&f[0]));
	throw eo;
}

double Delaunay::triangle_area(int index)
{
	boost::array<Vector2D,3> p;
	p[0]=cor[f[index].vertices[0]];
	p[1]=cor[f[index].vertices[1]];
	p[2]=cor[f[index].vertices[2]];
	double x1=p[2].x-p[0].x;
	double x2=p[1].x-p[0].x;
	double y1=p[2].y-p[0].y;
	double y2=p[1].y-p[0].y;
	return -0.5*(x1*y2-x2*y1);
}

void Delaunay::update(const vector<Vector2D>& points,vector<Vector2D>
	const& cpoints)
{
	if(logger)
		logger->output(cor,f);
	build_delaunay(points,cpoints);
}

int Delaunay::Walk(int point)
{ //This method starts from the last added facet and walks to the facet that contains point
	//returns the index of the facet that contains point;
	int cur_facet=lastFacet;
	int finish=0;
	boost::array<Vector2D,3> points;
	points[2]=cor[point];
	while(finish==0)
	{
		finish=1;
		//Test friends
		for(int i=0;i<3;++i)
		{
			points[0]=cor[f[cur_facet].vertices[i]];
			points[1]=cor[f[cur_facet].vertices[(i+1)%3]];
			if(orient2d(points)<0)
			{
				finish=0;
				cur_facet=f[cur_facet].neighbors[i];
				break;
			}
		}
	}
	lastFacet=cur_facet;
	return cur_facet;
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
		r = max(r,radius[vec[i]]);
	return 2*r;
}

void Delaunay::FindContainingTetras(int StartFacet,int point,vector<int> &result)
{
	result.clear();
	int PointLocation=FindPointInFacet(StartFacet,point);
	int NextFacet=f[StartFacet].neighbors[PointLocation];
	result.reserve(12);
	result.push_back(NextFacet);
	while(NextFacet!=StartFacet)
	{
		PointLocation=FindPointInFacet(NextFacet,point);
		NextFacet=f[NextFacet].neighbors[PointLocation];
		result.push_back(NextFacet);
	}
}

int Delaunay::FindPointInFacet(int facet,int point)
{
	for(int i=0;i<3;++i)
		if(f[facet].vertices[i]==point)
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
		for(int j=0;j<3;++j)
			if(f[facet].vertices[i]==(olength+j))
				return true;
	return false;
}

double Delaunay::CalculateRadius(int facet)
{
	const double big=1e10;
	const double a=cor[f[facet].vertices[0]].distance(cor[f[facet].vertices[1]]);
	const double b=cor[f[facet].vertices[0]].distance(cor[f[facet].vertices[2]]);
	const double c=cor[f[facet].vertices[2]].distance(cor[f[facet].vertices[1]]);
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
		if(!InCell(cell_points,cor[i]))
		{
			UniversalError eo("Mesh point outside cell");
			eo.AddEntry("Point number",i);
			eo.AddEntry("X location",cor[i].x);
			eo.AddEntry("Y location",cor[i].y);
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
	return radius[facet];
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

facet* Delaunay::get_facet(int index)
{
	return &f[index];
}

double Delaunay::get_facet_coordinate(int Facet,int vertice, int dim)
{
	if(dim==0)
		return cor[f[Facet].vertices[vertice]].x;
	else
		return cor[f[Facet].vertices[vertice]].y;
}

Vector2D Delaunay::get_point(int index) const
{
	return cor[index];
}

double Delaunay::get_cor(int index, int dim) const
{
	if(dim==0)
		return cor[index].x;
	else if(dim==1)
		return cor[index].y;
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
	cor[index]=p;
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
		cor.push_back(points[i]);
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
		if(IsOuterFacet(f[index].neighbors[i]))
			++counter;
	if(counter>1)
		return true;
	else
		return false;
}

int Delaunay::FindTripleLoc(facet const& fct)const
{
	for(size_t i=0;i<3;++i)
		if(!IsOuterFacet(fct.neighbors[i]))
			return (i+1)%3;
	throw UniversalError("Trouble in constructing boundary triangles. No inner neighbor");
}

namespace
{
	bool IsOuterQuick(facet const& f,int olength)
	{
		for(size_t i=0;i<3;++i)
			if(f.vertices[i]>=olength)
				return true;
		return false;
	}

	bool IsEdgeFacet(vector<facet> const& facets,facet const& f,int olength)
	{
		int counter=0;
		for(size_t i=0;i<3;++i)
		{
			if(f.vertices[i]>=olength)
				return false;
			if(IsOuterQuick(facets[f.neighbors[i]],olength))
				++counter;
		}
		if(counter>0)
			return true;
		else
			return false;
	}

	void AddAllFacets(vector<facet> const& f,int olength,vector<int> const& tocheck,
		vector<int> &res,vector<int> &res_facet)
	{
		vector<int> temp;
		vector<int> ftemp;
		temp.reserve(30);
		for(size_t i=0;i<tocheck.size();++i)
		{
			for(size_t j=0;j<3;++j)
				if(f[tocheck[i]].vertices[j]<olength)
				{
					temp.push_back(f[tocheck[i]].vertices[j]);
					ftemp.push_back(tocheck[i]);
				}
		}
		sort(temp.begin(),temp.end());
		sort(ftemp.begin(),ftemp.end());
		temp=unique(temp);
		ftemp=unique(ftemp);
		res.insert(res.end(),temp.begin(),temp.end());
		res_facet.insert(res_facet.end(),ftemp.begin(),ftemp.end());
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

vector<vector<int> > Delaunay::FindOuterPoints(vector<Edge> const& edges)
{
	// We add the points in a counter clockwise fashion
	vector<vector<int> > res(edges.size());
	if(olength<100)
	{
		for(size_t j=0;j<edges.size();++j)
		{
			res[j].resize(olength);
			for(int i=0;i<olength;++i)
				res[j][i]=i;
		}
		return res;
	}
	vector<int> res_temp,outer_points,f_temp,f_add(f.size(),0);
	res_temp.reserve((int)(20*sqrt(1.0*olength)));
	f_temp.reserve((int)(10*sqrt(1.0*olength)));
	outer_points.reserve((int)(10*sqrt(1.0*olength)));
	// Walk to an outer point
	int cur_facet=Walk(olength);
	// Find the real point
	int real_point = 0;
	for(int i=0;i<3;++i)
	{
		if(f[cur_facet].vertices[i]<olength)
		{
			real_point=f[cur_facet].vertices[i];
			break;
		}
	}
	vector<int> containing_facets;
	FindContainingTetras(cur_facet,real_point,containing_facets);
	for(size_t i=0;i<containing_facets.size();++i)
	{
		if(IsEdgeFacet(f,f[containing_facets[i]],olength))
		{
			cur_facet=containing_facets[i];
			break;
		}
	}
	int start_facet=cur_facet;
	int point_index=FindPointInFacet(cur_facet,real_point);
	if(IsOuterQuick(f[f[cur_facet].neighbors[point_index]],olength))
	{
		point_index=(point_index+1)%3;
		real_point=f[cur_facet].vertices[point_index];
	}
	if(IsOuterQuick(f[f[cur_facet].neighbors[point_index]],olength))
	{
		point_index=(point_index+1)%3;
		real_point=f[cur_facet].vertices[point_index];
	}
	do
	{
		FindContainingTetras(cur_facet,real_point,containing_facets);
		int old_current=cur_facet;
		for(size_t i=0;i<containing_facets.size();++i)
		{
			if(IsEdgeFacet(f,f[containing_facets[i]],olength)&&
				containing_facets[i]!=old_current)
				cur_facet=containing_facets[i];
			if(!IsOuterQuick(f[containing_facets[i]],olength))
				f_temp.push_back(containing_facets[i]);
		}
		point_index=(1+FindPointInFacet(cur_facet,real_point))%3;
		if(IsTripleOut(cur_facet))
			point_index=(point_index+1)%3;
		real_point=f[cur_facet].vertices[point_index];
	}while(start_facet!=cur_facet);
	sort(f_temp.begin(),f_temp.end());
	f_temp=unique(f_temp);
	// Found all initial outer facets

	//Find the points in the outer facets
	vector<vector<int> > toduplicate(edges.size());
	// Recursively look for more points
	vector<bool> checked(olength,false);
	for(size_t i=0;i<f_temp.size();++i)
		AddOuterFacets(f_temp[i],toduplicate,edges,checked);
	for(size_t i=0;i<edges.size();++i)
	{
		sort(toduplicate[i].begin(),toduplicate[i].end());
		toduplicate[i]=unique(toduplicate[i]);
	}
	return toduplicate;
}

void Delaunay::AddRigid(OuterBoundary const* /*obc*/,vector<Edge> const& edges,
	vector<vector<int> > &toduplicate)
{
	for(size_t i=0;i<edges.size();++i)
	{
		if(toduplicate[i].empty())
			continue;
		vector<Vector2D> toadd;
		toadd.reserve(toduplicate[i].size());
		//vector<int> pointstemp(toduplicate[i].size());
		Vector2D par(Parallel(edges[i]));
		par=par/abs(par);
		for(size_t j=0;j<toduplicate[i].size();++j)
		{
			Vector2D temp=cor[toduplicate[i][j]]-edges[i].vertices.first;
			temp=2*par*ScalarProd(par,temp)-temp+edges[i].vertices.first;
			toadd.push_back(temp);
			//	pointstemp[j]=j;
		}
		vector<int> order=HilbertOrder(toadd,(int)toadd.size());
		ReArrangeVector(toadd,order);
		AddBoundaryPoints(toadd);
		ReArrangeVector(toduplicate[i],order);
		//toduplicate[i]=pointstemp;
	}
}

namespace
{
	vector<vector<int> > FindCorners(vector<vector<int> > const& toduplicate)
	{
		vector<vector<int> > res(toduplicate.size());
		for(size_t i=0;i<toduplicate.size();++i)
		{
			for(size_t j=0;j<toduplicate[i].size();++j)
			{
				if(binary_search(toduplicate[(i+1)%4].begin(),toduplicate[(i+1)%4].end(),
					toduplicate[i][j]))
					res[i].push_back(toduplicate[i][j]);
			}
		}
		return res;
	}
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
		if(toduplicate[i].empty())
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
		toadd.reserve(toduplicate[i].size());
		//vector<int> pointstemp(toduplicate[i].size());
		for(size_t j=0;j<toduplicate[i].size();++j)
		{
			toadd.push_back(cor[toduplicate[i][j]]+change);
			//	pointstemp[j]=j;
		}
		vector<int> order=HilbertOrder(toadd,(int)toadd.size());
		ReArrangeVector(toadd,order);
		AddBoundaryPoints(toadd);
		ReArrangeVector(toduplicate[i],order);
		//toduplicate[i]=pointstemp;
	}
	// Done with sides do corners now
	vector<Edge> corneredges=GetCornerEdges(obc);
	vector<vector<int> > corners(toduplicate.size());
	for(size_t i=0;i<toduplicate.size();++i)
	{
		for(size_t j=0;j<toduplicate[i].size();++j)
		{
			const int facet_loc=Walk(toduplicate[i][j]);
			const Vector2D center=cor[toduplicate[i][j]];
			const double R=2*GetMaxRadius(toduplicate[i][j],facet_loc);
			if(CircleSegmentIntersect(corneredges[2*i],center,R))
				corners[i].push_back(toduplicate[i][j]);
			if(CircleSegmentIntersect(corneredges[(2*i+7)%8],center,R))
				corners[(i+3)%4].push_back(toduplicate[i][j]);
		}
	}
	for(size_t i=0;i<corners.size();++i)
	{
		if(corners[i].empty())
			continue;
		sort(corners[i].begin(),corners[i].end());
		corners[i]=unique(corners[i]);
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
		toadd.reserve(corners[i].size());
		//		vector<int> pointstemp(corners[i].size());
		for(size_t j=0;j<corners[i].size();++j)
		{
			toadd.push_back(cor[corners[i][j]]+change);
			//		pointstemp[j]=j;
		}
		vector<int> order=HilbertOrder(toadd,(int)toadd.size());
		ReArrangeVector(toadd,order);
		AddBoundaryPoints(toadd);
		ReArrangeVector(corners[i],order);
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
		if(toduplicate[i].empty())
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
		toadd.reserve(toduplicate[i].size());
		//vector<int> pointstemp(toduplicate[i].size());
		Vector2D par(Parallel(edges[i]));
		par=par/abs(par);
		for(size_t j=0;j<toduplicate[i].size();++j)
		{
			Vector2D temp=cor[toduplicate[i][j]];
			if(i%2==1)
			{
				temp-=edges[i].vertices.first;
				temp=2*par*ScalarProd(par,temp)-temp+edges[i].vertices.first;
			}
			toadd.push_back(temp+change);
			//pointstemp[j]=j;
		}
		vector<int> order=HilbertOrder(toadd,(int)toadd.size());
		ReArrangeVector(toadd,order);
		AddBoundaryPoints(toadd);
		ReArrangeVector(toduplicate[i],order);
		//toduplicate[i]=pointstemp;
	}
}

vector<vector<int> > Delaunay::BuildBoundary(OuterBoundary const* obc,vector<Edge> const& edges)
{
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
				toduplicate.push_back(corners[i]);
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
	facet const& F=f[index];
	double x1=cor[F.vertices[0]].x;
	double x2=cor[F.vertices[1]].x;
	double x3=cor[F.vertices[2]].x;
	double y1=cor[F.vertices[0]].y;
	double y2=cor[F.vertices[1]].y;
	double y3=cor[F.vertices[2]].y;
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
		res=max(res,radius[neigh[i]]);
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
			if(checked[f[cur_facet].vertices[i]]||(f[cur_facet].vertices[i]>=olength))
				continue;
			vector<int> neigh=FindContainingTetras(cur_facet,f[cur_facet].vertices[i]);
			for(size_t k=0;k<neigh.size();++k)
			{
				Vector2D center=GetCircleCenter(neigh[k]);
				for(size_t l=0;l<edges.size();++l)
				{
					if(CircleSegmentIntersect(edges[l],center,radius[neigh[k]]))
					{
						toduplicate[l].push_back(f[cur_facet].vertices[i]);
						added=true;
					}
				}
			}
			checked[f[cur_facet].vertices[i]]=true;
			if(added)
			{
				for(size_t j=0;j<neigh.size();++j)
				{
					if(!IsOuterQuick(f[neigh[j]],olength))
						tocheck.push(neigh[j]);
				}
			}
		}
	}
}
