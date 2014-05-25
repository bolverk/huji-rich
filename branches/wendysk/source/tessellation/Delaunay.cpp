#include "Delaunay.hpp"
#include <iostream>
#include <vector>
#include <cmath>

Delaunay::DataOnlyForBuild::DataOnlyForBuild():insert_order(vector<int> ()),
	copied(vector<vector<char> > ()),BoundaryCandidates(vector<int> ())
{}

Delaunay::DataOnlyForBuild::DataOnlyForBuild(DataOnlyForBuild const& other):
insert_order(other.insert_order),copied(other.copied),BoundaryCandidates
(other.BoundaryCandidates){}

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
  tree(0),
  treePoints(0),
  NewPointIndex(vector<int>()),
  lastFacet(0),CalcRadius(false),
  radius(vector<double>()),
  PointWasAdded(false),
  BIG(0), last_facet_added(0),
  f(vector<facet>()),
  cor(vector<Vector2D>()),
  totalCor(vector<Vector2D>()),
  length(0),
  olength(0),location_pointer(0), last_loc(0),bc(0)
{}

Delaunay::Delaunay(Delaunay const& other):
  tree(0),
  treePoints(0),
  NewPointIndex(other.NewPointIndex),
  lastFacet(other.lastFacet),
  CalcRadius(other.CalcRadius),
  radius(other.radius),
  PointWasAdded(other.PointWasAdded),
  BIG(other.BIG),
  last_facet_added(other.last_facet_added),
  f(other.f),
  cor(other.cor),
  totalCor(other.totalCor),
  length(other.length),
  olength(other.olength),
  location_pointer(other.location_pointer),
  last_loc(other.last_loc),
  bc(other.bc) {}

Delaunay::~Delaunay(void)
{
	//	cout<<"Entered d destruct"<<endl;
	cor.clear();
	f.clear();
	if(treePoints!=0)
	{
		annDeallocPts(treePoints);
		delete tree;
		annClose();	
	}
	NewPointIndex.clear();
	totalCor.clear();
}

namespace{
  int find_index(facet const& fc, int i)
  {
    for(int j=0;j<3;++j)
      {
	if(fc.get_friend(j)==i)
	  return j;
      }
    cout<<"Couldn't find number "<<i<<" in facet";
    return 0;
  }
}

void Delaunay::add_point(int index)
{
	int triangle=Walk(index);
	boost::array<int,3> outer,temp_friends;
	facet f_temp;
	int i=0;
	for(i=0;i<3;++i)
	{
		outer[i]=f[triangle].get_vertice(i);
		temp_friends[i]=f[triangle].get_friend(i);
	}
	// create and _update the new facets
	f.push_back(f_temp);
	f.push_back(f_temp);
	f[triangle].set_vertice(outer[2],0);
	f[triangle].set_vertice(outer[0],1);
	f[triangle].set_vertice(index,2);
	f[location_pointer+1].set_vertice(outer[0],0);
	f[location_pointer+1].set_vertice(outer[1],1);
	f[location_pointer+1].set_vertice(index,2);
	f[location_pointer+2].set_vertice(outer[1],0);
	f[location_pointer+2].set_vertice(outer[2],1);
	f[location_pointer+2].set_vertice(index,2);

	f[triangle].set_friend(temp_friends[2],0);
	f[triangle].set_friend(location_pointer+1,1);
	f[triangle].set_friend(location_pointer+2,2);
	f[location_pointer+1].set_friend(temp_friends[0],0);
	f[location_pointer+1].set_friend(location_pointer+2,1);
	f[location_pointer+1].set_friend(triangle,2);
	f[location_pointer+2].set_friend(temp_friends[1],0);
	f[location_pointer+2].set_friend(triangle,1);
	f[location_pointer+2].set_friend(location_pointer+1,2);
	// _update the friends list of the friends
	if(temp_friends[1]!=last_loc)
	{
		i=find_index(f[temp_friends[1]],triangle);
		f[temp_friends[1]].set_friend(location_pointer+2,i);
	}
	if(temp_friends[0]!=last_loc)
	{
		i=find_index(f[temp_friends[0]],triangle);
		f[temp_friends[0]].set_friend(location_pointer+1,i);
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
	while(flip_stack.size()>0)
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
				circle_test[k]=cor[f[indexes[0]].get_vertice(k)];
			}
			circle_test[3]=cor[check[0]];

			if(incircle(circle_test)>0)
			{
				//The point is in a circle change the facets and their friends
				const int v1=f[indexes[0]].get_vertice((other[1]+1)%3);
				const int f1=f[indexes[0]].get_friend(other[1]);
				const int f12=f[indexes[0]].get_friend((other[1]+2)%3);
				const int f13=f[indexes[0]].get_friend((other[1]+2)%3);
				const int v2=f[indexes[1]].get_vertice((check[1]+1)%3);
				const int f2=f[indexes[1]].get_friend((check[1]+2)%3);
				const int f22=f[indexes[1]].get_friend(check[1]);
				const int f23=f[indexes[1]].get_friend((check[1]+2)%3);
				f[indexes[0]].set_vertice(other[0],0);
				f[indexes[0]].set_vertice(v1,1);
				f[indexes[0]].set_vertice(check[0],2);
				f[indexes[1]].set_vertice(check[0],0);
				f[indexes[1]].set_vertice(v2,1);
				f[indexes[1]].set_vertice(other[0],2);
				f[indexes[0]].set_friend(f1,0);
				f[indexes[0]].set_friend(f2,1);
				f[indexes[0]].set_friend(indexes[1],2);
				f[indexes[1]].set_friend(f22,0);
				f[indexes[1]].set_friend(f12,1);
				f[indexes[1]].set_friend(indexes[0],2);
				// change the friends of the friends if needed
				if(f23!=last_loc)
				{
					f[f23].set_friend(indexes[0],find_index(f[f23],indexes[1]));
				}
				if(f13!=last_loc)
				{
					f[f13].set_friend(indexes[1],find_index(f[f13],indexes[0]));
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
				array_temp[1]=f[indexes[1]].get_friend(0);
				flip_stack.push(array_temp);
				array_temp[0]=indexes[0];
				array_temp[1]=f[indexes[0]].get_friend(1);
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

void Delaunay::build_delaunay(vector<Vector2D>const& vp,OuterBoundary const* bc1)
{
	DataOnlyForBuild data;
	// copy the boundary conditions
	bc=bc1;
	lastFacet=0;
	CalcRadius=false;
	double maxV=(bc->GetGridBoundary(Right)-bc->GetGridBoundary(Left)),
		maxH=(bc->GetGridBoundary(Up)-bc->GetGridBoundary(Down));
	BIG=pow(max(maxV,maxH),2);
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
	double width=bc->GetGridBoundary(Right)-bc->GetGridBoundary(Left);
	double height=bc->GetGridBoundary(Up)-bc->GetGridBoundary(Down);
	p_temp.set_x(bc->GetGridBoundary(Left)-5.5*width);
	p_temp.set_y(bc->GetGridBoundary(Down)-height);
	cor.push_back(p_temp);
	p_temp.set_x(bc->GetGridBoundary(Right)+5.5*width);
	p_temp.set_y(bc->GetGridBoundary(Down)-height);
	cor.push_back(p_temp);
	p_temp.set_x((bc->GetGridBoundary(Left)+bc->GetGridBoundary(Right))/2.0);
	p_temp.set_y(bc->GetGridBoundary(Up)+2*height);
	cor.push_back(p_temp);
	// Create the big triangle, and assign friends
	facet f_temp;
	f.push_back(f_temp);
	f[0].set_vertice(len,0);
	f[0].set_vertice(len+1,1);
	f[0].set_vertice(len+2,2);
	for(int i=0;i<3;i++)
		f[0].set_friend(last_loc,i);
	location_pointer=0;
	// add the points
	for(int i=0;i<length-3;i++)
	{
		add_point(data.insert_order[i]);
	}
	BuildBoundary(data);

}

void Delaunay::find_diff(facet *f1,facet *f2,int *p) const
{
	bool counter;
	int i=0;
	for(;i<3;++i)
	{
		counter=false;
		for(int j=0;j<3;++j)
		{
			if(f1->get_vertice(i)==f2->get_vertice(j))
			{
				counter=true;
				break;
			}
		}
		if(counter==false)
		{
			p[0]=f1->get_vertice(i);
			p[1]=i;
			return;
		}
	}
	UniversalError eo("Delaunay, Couldn't find difference bewteen two facets");
	eo.AddEntry("Facet 1",double(f1-&f[0]));
	eo.AddEntry("Facet 2",double(f2-&f[0]));
	throw eo;
}

void Delaunay::output()
{	
	fstream myFile ("c:\\Delaunay.bin",ios::out | ios::binary);
	length-=3;
	int int_temp=int(cor.size());
	myFile.write ((char*)&int_temp, sizeof (int));
	length+=3;
	int_temp=int(f.size());
	myFile.write ((char*)&(int_temp), sizeof (int));
	double temp;
	for(int i=0;i<(int)cor.size();i++)
	{
		temp=cor[i].x;
		myFile.write ((char*)&temp, sizeof (double));
		temp=cor[i].y;
		myFile.write ((char*)&temp, sizeof (double));
	}
	int temp2,j;

	for(int i=0;i<(int)f.size();i++)
	{
		{
			for(j=0;j<3;j++)
			{
				temp2=f[i].get_vertice(j);
				myFile.write ((char*)&temp2, sizeof (int));
			}
		}
	}
	for(int i=0;i<(int)f.size();i++)
	{
		{
			for(j=0;j<3;j++)
			{
				temp2=f[i].get_friend(j);
				myFile.write ((char*)&temp2, sizeof (int));
			}
		}
	}
	myFile.close();
}

double Delaunay::triangle_area(int index)
{
	boost::array<Vector2D,3> p;
	p[0]=cor[f[index].get_vertice(0)];
	p[1]=cor[f[index].get_vertice(1)];
	p[2]=cor[f[index].get_vertice(2)];
	double x1=p[2].x-p[0].x;
	double x2=p[1].x-p[0].x;
	double y1=p[2].y-p[0].y;
	double y2=p[1].y-p[0].y;
	return -0.5*(x1*y2-x2*y1);
}

void Delaunay::update(const vector<Vector2D>& points)
{
	build_delaunay(points,bc);
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
			points[0]=cor[f[cur_facet].get_vertice(i)];
			points[1]=cor[f[cur_facet].get_vertice((i+1)%3)];
			if(orient2d(points)<0)
			{
				finish=0;
				cur_facet=f[cur_facet].get_friend(i);
				break;
			}
		}
	}
	lastFacet=cur_facet;
	return cur_facet;
}

void Delaunay::BuildTree(DataOnlyForBuild &data)
{
	// Find the points to add
	data.BoundaryCandidates.reserve(int(20*sqrt(1.0*(double)cor.size())));
	int N=olength;
	for(int i=0;i<N;++i)
	{
		double R=FindMaxRadius(data.insert_order[i]);
		if(BoundaryCandidate(R,data.insert_order[i]))
			data.BoundaryCandidates.push_back(data.insert_order[i]);
	}
	// Build kd-tree
	N=(int)data.BoundaryCandidates.size();
	treePoints=annAllocPts(N,2);
	for(int i=0;i<N;++i)
	{
		treePoints[i][0]=cor[data.BoundaryCandidates[i]].x;
		treePoints[i][1]=cor[data.BoundaryCandidates[i]].y;
	}
	tree=new ANNkd_tree(treePoints,N,2,1,ANN_KD_SUGGEST);	
}

void Delaunay::BuildBoundary(DataOnlyForBuild &data)
{
	NewPointIndex.clear();
	NewPointIndex.reserve(8*(int)sqrt(double(cor.size())*1.2));

	DuplicatePoints(data);

	annDeallocPts(treePoints);
	delete tree;
	annClose();	
}

double Delaunay::FindMaxRadius(int point)
{
	int startFacet=Walk(point);
	vector<int> vec;
	FindContainingTetras(startFacet,point,vec);
	double r=0,temp;
	int n=(int)vec.size();
	for(int i=0;i<n;++i)
	{
		temp=radius[vec[i]];
		if(temp>r)
			r=temp;
	}
	return 2*r;
}

void Delaunay::FindContainingTetras(int StartFacet,int point,vector<int> &result)
{
	int PointLocation=FindPointInFacet(StartFacet,point);
	int NextFacet=f[StartFacet].get_friend(PointLocation);
	result.reserve(8);
	result.push_back(NextFacet);
	while(NextFacet!=StartFacet)
	{
		PointLocation=FindPointInFacet(NextFacet,point);
		NextFacet=f[NextFacet].get_friend(PointLocation);
		result.push_back(NextFacet);
	}
}

int Delaunay::FindPointInFacet(int facet,int point)
{
	for(int i=0;i<3;++i)
		if(f[facet].get_vertice(i)==point)
			return i;
	UniversalError eo("Error in Delaunay, FindPointInFacet");
	eo.AddEntry("Facet number",facet);
	eo.AddEntry("Point number",point);
	throw eo;
}

bool Delaunay::IsOuterFacet(int facet)
{
	//int PointNum=length-1;
	for(int i=0;i<3;++i)
		for(int j=0;j<3;++j)
			if(f[facet].get_vertice(i)==(olength+j))
				return true;
	return false;
}

double Delaunay::CalculateRadius(int facet)
{
	double a=cor[f[facet].get_vertice(0)].distance(cor[f[facet].get_vertice(1)]);
	double b=cor[f[facet].get_vertice(0)].distance(cor[f[facet].get_vertice(2)]);
	double c=cor[f[facet].get_vertice(2)].distance(cor[f[facet].get_vertice(1)]);
	return a*b*c/sqrt((a+b+c)*(b+c-a)*(a-b+c)*(a+b-c));
}

double Delaunay::ReflectDistance(double SearchR,Sides side,int point,int maxPoints,int &closest)
{
	// Finds the distance from point point reflecting at side and looking at the 
	// numOfPoints closest neighbor. Also returns the index of the closest point
	Vector2D newPoint;
	ANNpoint queryPt;
	boost::array<double,2> pointToCheck;
	switch(side)
	{
	case(LU):
	  throw UniversalError("I do not know what to do in case LU in Delaunay::ReflectDistance");
	case(LD):
	  throw UniversalError("I do not know what to do in case LD in Delaunay::ReflectDistance");
	case(RU):
	  throw UniversalError("I do not know what to do in case RU in Delaunay::ReflectDistance");
	case(RD):
	  throw UniversalError("I do not know what to do in case RD in Delaunay::ReflectDistance");
	case(RIGHT):
		{
			pointToCheck[0]=-cor[point].x+2*bc->GetGridBoundary(Right);
			pointToCheck[1]=cor[point].y;
			if((pointToCheck[0]-SearchR)>bc->GetGridBoundary(Right))
			{
				closest=point;
				return 2*BIG;
			}
			break;
		}
	case(LEFT):
		{
			pointToCheck[0]=-cor[point].x+2*bc->GetGridBoundary(Left);
			pointToCheck[1]=cor[point].y;
			if((pointToCheck[0]+SearchR)<bc->GetGridBoundary(Left))
			{
				closest=point;
				return 2*BIG;
			}
			break;
		}
	case(DOWN):
		{
			pointToCheck[0]=cor[point].x;
			pointToCheck[1]=2*bc->GetGridBoundary(Down)-cor[point].y;
			if((pointToCheck[1]+SearchR)<bc->GetGridBoundary(Down))
			{
				closest=point;
				return 2*BIG;
			}
			break;
		}
	case(UP):
		{
			pointToCheck[0]=cor[point].x;
			pointToCheck[1]=2*bc->GetGridBoundary(Up)-cor[point].y;
			if((pointToCheck[1]-SearchR)>bc->GetGridBoundary(Up))
			{
				closest=point;
				return 2*BIG;
			}
			break;
		}
	default:
		UniversalError eo("Error in Delaunay::ReflectDistance: Expected UP, DOWN, LEFT or RIGHT, but got a different input");
		throw eo;
	}
	queryPt=annAllocPt(2);
	queryPt[0]=pointToCheck[0];
	queryPt[1]=pointToCheck[1];
	ANNidxArray nnIdx; // near neighbor indices
	ANNdistArray dists; // near neighbor distances
	nnIdx = new ANNidx[maxPoints]; // allocate near neigh indices
	dists = new ANNdist[maxPoints]; // allocate near neighbor dists
	tree->annkFRSearch(queryPt,SearchR*SearchR,maxPoints,nnIdx,dists,0);
	closest=nnIdx[maxPoints-1];
	double res=dists[maxPoints-1];

	//cleanup
	delete []nnIdx;
	delete []dists;
	annDeallocPt(queryPt);

	return res;
}

void Delaunay::AddBoundaryPoint(int pointToAdd,Sides side,bool partial,DataOnlyForBuild
	&data)
{
	bool AddedNewPoint=false;
	vector<Vector2D> *vecPoint;
	if(partial)
		vecPoint=&totalCor;
	else
		vecPoint=&cor;
	Vector2D newPoint;
	switch(side)
	{
	case(LD):
		{
			newPoint.Set(-bc->GetGridBoundary(Right)+bc->GetGridBoundary(Left)+
				vecPoint->at(pointToAdd).x,	vecPoint->at(pointToAdd).y+
				bc->GetGridBoundary(Down)-bc->GetGridBoundary(Up));
			NewPointIndex.push_back(pointToAdd);
			AddedNewPoint=true;
			break;
		}
	case(LU):
		{
			newPoint.Set(-bc->GetGridBoundary(Right)+bc->GetGridBoundary(Left)+
				vecPoint->at(pointToAdd).x,	vecPoint->at(pointToAdd).y-bc->GetGridBoundary(Down)+
				bc->GetGridBoundary(Up));
			NewPointIndex.push_back(pointToAdd);
			AddedNewPoint=true;
			break;
		}
	case(RD):
		{
			newPoint.Set(bc->GetGridBoundary(Right)-bc->GetGridBoundary(Left)+
				vecPoint->at(pointToAdd).x,	vecPoint->at(pointToAdd).y+bc->GetGridBoundary(Down)-
				bc->GetGridBoundary(Up));
			NewPointIndex.push_back(pointToAdd);
			AddedNewPoint=true;
			break;
		}
	case(RU):
		{
			newPoint.Set(bc->GetGridBoundary(Right)-bc->GetGridBoundary(Left)+
				vecPoint->at(pointToAdd).x,	vecPoint->at(pointToAdd).y-bc->GetGridBoundary(Down)+
				bc->GetGridBoundary(Up));
			NewPointIndex.push_back(pointToAdd);
			AddedNewPoint=true;
			break;
		}
	case(RIGHT):
		{
			if(bc->GetBoundaryType()==Rectengular)
				newPoint.Set(2*bc->GetGridBoundary(Right)-vecPoint->at(pointToAdd).x,
				vecPoint->at(pointToAdd).y);
			else
			{
				newPoint.Set(bc->GetGridBoundary(Right)-bc->GetGridBoundary(Left)+vecPoint->at(pointToAdd).x,
					vecPoint->at(pointToAdd).y);
				NewPointIndex.push_back(pointToAdd);
				AddedNewPoint=true;
			}
			break;
		}
	case(LEFT):
		{
			if(bc->GetBoundaryType()==Rectengular)
				newPoint.Set(2*bc->GetGridBoundary(Left)-vecPoint->at(pointToAdd).x,
				vecPoint->at(pointToAdd).y);
			else
			{
				newPoint.Set(bc->GetGridBoundary(Left)-bc->GetGridBoundary(Right)+vecPoint->at(pointToAdd).x,
					vecPoint->at(pointToAdd).y);
				NewPointIndex.push_back(pointToAdd);
				AddedNewPoint=true;
			}
			break;
		}
	case(DOWN):
		{
			if(bc->GetBoundaryType()!=Periodic)
			{
				newPoint.Set(vecPoint->at(pointToAdd).x,
					2*bc->GetGridBoundary(Down)-vecPoint->at(pointToAdd).y);
				if(bc->GetBoundaryType()==HalfPeriodic)
				{
					NewPointIndex.push_back(pointToAdd);
					AddedNewPoint=true;
				}
			}
			else
			{
				newPoint.Set(vecPoint->at(pointToAdd).x,
					bc->GetGridBoundary(Down)-bc->GetGridBoundary(Up)+vecPoint->at(pointToAdd).y);
				NewPointIndex.push_back(pointToAdd);
				AddedNewPoint=true;
			}
			break;
		}
	case(UP):
		{
			if(bc->GetBoundaryType()!=Periodic)
			{
				newPoint.Set(vecPoint->at(pointToAdd).x,
					2*bc->GetGridBoundary(Up)-vecPoint->at(pointToAdd).y);
				if(bc->GetBoundaryType()==HalfPeriodic)
				{
					NewPointIndex.push_back(pointToAdd);
					AddedNewPoint=true;
				}
			}
			else
			{
				newPoint.Set(vecPoint->at(pointToAdd).x,
					-bc->GetGridBoundary(Down)+bc->GetGridBoundary(Up)+vecPoint->at(pointToAdd).y);
				NewPointIndex.push_back(pointToAdd);
				AddedNewPoint=true;
			}
			break;
		}
	default:
		UniversalError eo("Error in Delaunay::AddBoundaryPoint: got unexpected input");
		throw eo;
	}	
	if(data.copied[side][pointToAdd]==0)
	{
		cor.push_back(newPoint);
		add_point(int(cor.size())-1);
		data.copied[side][pointToAdd]=1;
		PointWasAdded=true;
	}
	else
		if(AddedNewPoint)
			NewPointIndex.pop_back();
}

int Delaunay::AddSphere(double minR,double maxR,ANNpoint &queryPt,Sides side,
	bool partial,DataOnlyForBuild &data)
{
	// Get how many point are in the sphere
	int NumberInSphere=tree->annkFRSearch(queryPt,maxR*maxR,0);
	//allocate
	ANNidxArray nnIdx; // near neighbor indices
	ANNdistArray dists; // near neighbor distances
	nnIdx = new ANNidx[NumberInSphere]; // allocate near neigh indices
	dists = new ANNdist[NumberInSphere]; // allocate near neighbor dists
	// Get the points
	tree->annkFRSearch(queryPt,maxR*maxR,NumberInSphere,nnIdx,dists,0);
	for(int i=0;i<NumberInSphere;++i)
	{
		if(dists[i]>(minR*minR))
			AddBoundaryPoint(data.BoundaryCandidates[nnIdx[i]],side,partial,data); // add the point
	}
	delete [] nnIdx;
	delete [] dists;
	return NumberInSphere;
}

bool Delaunay::GetQueryPointReflective(Vector2D &qpoint,Sides side,double SearchR,int point)
{
	switch(side)
	{
	case(LU):
	  throw UniversalError("I do not know what to do in case LU in Delaunay::GetQueryPointReflective");
	case(LD):
	  throw UniversalError("I do not know what to do in case LD in Delaunay::GetQueryPointReflective");
	case(RU):
	  throw UniversalError("I do not know what to do in case RU in Delaunay::GetQueryPointReflective");
	case(RD):
	  throw UniversalError("I do not know what to do in case RD in Delaunay::GetQueryPointReflective");
	case(RIGHT):
		{
			qpoint.x=-cor[point].x+2*bc->GetGridBoundary(Right);
			qpoint.y=cor[point].y;
			if((qpoint.x-SearchR)>bc->GetGridBoundary(Right))
			{
				return false;
			}
			break;
		}
	case(LEFT):
		{
			qpoint.x=-cor[point].x+2*bc->GetGridBoundary(Left);
			qpoint.y=cor[point].y;
			if((qpoint.x+SearchR)<bc->GetGridBoundary(Left))
			{
				return false;
			}
			break;
		}
	case(DOWN):
		{
			qpoint.x=cor[point].x;
			qpoint.y=2*bc->GetGridBoundary(Down)-cor[point].y;
			if((qpoint.y+SearchR)<bc->GetGridBoundary(Down))
			{
				return false;
			}
			break;
		}
	case(UP):
		{
			qpoint.x=cor[point].x;
			qpoint.y=2*bc->GetGridBoundary(Up)-cor[point].y;
			if((qpoint.y-SearchR)>bc->GetGridBoundary(Up))
			{
				return false;
			}
			break;
		}
	default:
		UniversalError eo("Error in Delaunay::GetQueryPoint got unexpected input");
		throw eo;
	}
	return true;
}

bool Delaunay::GetQueryPointPeriodic(Vector2D &qpoint,Sides side,double SearchR,int point)
{
	switch(side)
	{
	case(RIGHT):
		{
			qpoint.x=cor[point].x-bc->GetGridBoundary(Right)
				+bc->GetGridBoundary(Left);
			qpoint.y=cor[point].y;
			if((qpoint.x+SearchR)<bc->GetGridBoundary(Left))
			{
				return false;
			}
			break;
		}
	case(LEFT):
		{
			qpoint.x=cor[point].x+bc->GetGridBoundary(Right)
				-bc->GetGridBoundary(Left);
			qpoint.y=cor[point].y;
			if((qpoint.x-SearchR)>bc->GetGridBoundary(Right))
			{
				return false;
			}
			break;
		}
	case(DOWN):
		{
			qpoint.x=cor[point].x;
			qpoint.y=bc->GetGridBoundary(Up)-bc->GetGridBoundary(Down)
				+cor[point].y;
			if((qpoint.y-SearchR)>bc->GetGridBoundary(Up))
			{
				return false;
			}
			break;
		}
	case(UP):
		{
			qpoint.x=cor[point].x;
			qpoint.y=bc->GetGridBoundary(Down)-bc->GetGridBoundary(Up)
				+cor[point].y;
			if((qpoint.y+SearchR)<bc->GetGridBoundary(Down))
			{
				return false;
			}
			break;
		}
	case(LU):
		{
			qpoint.x=cor[point].x+bc->GetGridBoundary(Right)
				-bc->GetGridBoundary(Left);
			qpoint.y=bc->GetGridBoundary(Down)-bc->GetGridBoundary(Up)
				+cor[point].y;
			if(((qpoint.y+SearchR)<bc->GetGridBoundary(Down))||
				((qpoint.x-SearchR)>bc->GetGridBoundary(Right)))
			{
				return false;
			}
			break;
		}
	case(LD):
		{
			qpoint.x=cor[point].x+bc->GetGridBoundary(Right)
				-bc->GetGridBoundary(Left);
			qpoint.y=-bc->GetGridBoundary(Down)+bc->GetGridBoundary(Up)
				+cor[point].y;
			if(((qpoint.y-SearchR)>bc->GetGridBoundary(Up))||
				((qpoint.x-SearchR)>bc->GetGridBoundary(Right)))
			{
				return false;
			}
			break;
		}
	case(RD):
		{
			qpoint.x=cor[point].x-bc->GetGridBoundary(Right)
				+bc->GetGridBoundary(Left);
			qpoint.y=-bc->GetGridBoundary(Down)+bc->GetGridBoundary(Up)
				+cor[point].y;
			if(((qpoint.y-SearchR)>bc->GetGridBoundary(Up))||
				((qpoint.x+SearchR)<bc->GetGridBoundary(Left)))
			{
				return false;
			}
			break;
		}
	case(RU):
		{
			qpoint.x=cor[point].x-bc->GetGridBoundary(Right)
				+bc->GetGridBoundary(Left);
			qpoint.y=bc->GetGridBoundary(Down)-bc->GetGridBoundary(Up)
				+cor[point].y;
			if(((qpoint.y+SearchR)<bc->GetGridBoundary(Down))||
				((qpoint.x+SearchR)<bc->GetGridBoundary(Left)))
			{
				return false;
			}
			break;
		}
	default:
		UniversalError eo("Error in Delaunay::GetQueryPoint got unexpected input");
		throw eo;
	}
	return true;
}

void Delaunay::CheckInput(void)
{
	try
	{
		if(bc->GetBoundaryType()==Rectengular)
		{
			for(int i=0;i<(int)cor.size();++i)
			{
				if((cor[i].x>bc->GetGridBoundary(Right))||(cor[i].x<bc->GetGridBoundary(Left)))
				{	
					throw UniversalError("Mesh point to outside the x axis");
				}
				if((cor[i].y>bc->GetGridBoundary(Up))||(cor[i].y<bc->GetGridBoundary(Down)))
				{
					UniversalError eo("Mesh point to outside the y axis");
					eo.AddEntry("Point number",i);
					eo.AddEntry("X coordinate",cor[i].x);
					eo.AddEntry("Y coordinate",cor[i].y);
					throw eo;
				}
			}
		}
		if(bc->GetBoundaryType()==Periodic)
		{
			for(int i=0;i<(int)cor.size();++i)
			{
				if(cor[i].x>bc->GetGridBoundary(Right))
					cor[i].x=cor[i].x-bc->GetGridBoundary(Right)+bc->GetGridBoundary(Left);
				if(cor[i].x<bc->GetGridBoundary(Left))
					cor[i].x=cor[i].x+bc->GetGridBoundary(Right)-bc->GetGridBoundary(Left);
				if(cor[i].y>bc->GetGridBoundary(Up))
					cor[i].y=cor[i].y-bc->GetGridBoundary(Up)+bc->GetGridBoundary(Down);
				if(cor[i].y<bc->GetGridBoundary(Down))
					cor[i].y=cor[i].y+bc->GetGridBoundary(Up)-bc->GetGridBoundary(Down);
			}
		}
		if(bc->GetBoundaryType()==HalfPeriodic)
		{
			for(int i=0;i<(int)cor.size();++i)
			{
				if(cor[i].x>bc->GetGridBoundary(Right))
					cor[i].x=cor[i].x-bc->GetGridBoundary(Right)+bc->GetGridBoundary(Left);
				if(cor[i].x<bc->GetGridBoundary(Left))
					cor[i].x=cor[i].x+bc->GetGridBoundary(Right)-bc->GetGridBoundary(Left);
				if((cor[i].y>bc->GetGridBoundary(Up))||(cor[i].y<bc->GetGridBoundary(Down)))
				{
					UniversalError eo("Mesh point to outside the y axis");
					eo.AddEntry("Point number",i);
					eo.AddEntry("X coordinate",cor[i].x);
					eo.AddEntry("Y coordinate",cor[i].y);
					throw eo;
				}
			}
		}
	}
	catch(UniversalError& eo)
	{
		throw;
	}
}

int Delaunay::GetOriginalIndex(int NewPoint) const
{
	if(!NewPointIndex.empty())
		if(NewPoint>=length)
			return NewPointIndex[NewPoint-length];
		else
			if(NewPoint>olength)
				return NewPointIndex[NewPoint-olength-3];
			else
				return NewPoint;
	else
		return NewPoint;
}

double Delaunay::GetInitialSearchR(int index,double &maxradius)
{
	int startFacet=Walk(index);
	vector<int> neigh;
	FindContainingTetras(startFacet,index,neigh);
	double res=1;
	maxradius=0;
	size_t n=neigh.size();
	for(size_t i=0;i<n;++i)
	{
		res*=1.5*radius[neigh[i]];
		maxradius=max(maxradius,radius[neigh[i]]);
	}
	maxradius*=2;
	return pow(res,1.0/(double)neigh.size());
}

bool Delaunay::BoundaryCandidate(double maxR,int point)
{
	if(cor[point].x-bc->GetGridBoundary(Left)<4*maxR)
		return true;
	if(cor[point].y-bc->GetGridBoundary(Down)<4*maxR)
		return true;
	if(bc->GetGridBoundary(Up)-cor[point].y<4*maxR)
		return true;
	if(bc->GetGridBoundary(Right)-cor[point].x<4*maxR)
		return true;
	return false;
}

void Delaunay::DuplicatePoints(DataOnlyForBuild &data)
{
	int sidesNumber,point;
	if(bc->GetBoundaryType()!=Periodic) // How many duplicates do I need to check
		sidesNumber=4;
	else
		sidesNumber=8;
	data.copied.clear();
	data.copied.resize(sidesNumber);
	for(int i=0;i<sidesNumber;++i)
		data.copied[i].assign(olength,0);
	vector<Sides> SidesToCheck;
	SidesToCheck.resize(sidesNumber);
	for(int i=0;i<sidesNumber;++i)
		SidesToCheck[i]=(Sides)i;
	Vector2D qpoint;
	ANNpoint queryPt;
	queryPt=annAllocPt(2);
	// Copy the outer points
	vector<vector<int> > outerpoints;
	GetOuterPoints2(outerpoints);
	int n;
	for(int i=0;i<4;++i)
	{
	  n=int(outerpoints[i].size());
		for(int j=0;j<n;++j)
		{
			if(bc->GetBoundaryType()==Rectengular)
				AddBoundaryPoint(outerpoints[i][j],(Sides)i,false,data);
			else
				if(bc->GetBoundaryType()==Periodic)
					AddBoundaryPoint(outerpoints[i][j],(Sides)((i+2)%4),false,data);
				else
					if(i==1||i==3)
						AddBoundaryPoint(outerpoints[i][j],(Sides)i,false,data);
					else
						AddBoundaryPoint(outerpoints[i][j],(Sides)((i+2)%4),false,
						data);
		}
	}
	// Calculate radius
	radius.resize(f.size());
	n=int(f.size());
	for(int i=0;i<n;++i)
		radius[i]=CalculateRadius(i);
	/*vector<double> searchr(olength);
	for(int i=0;i<olength;++i)
	{
		searchr[data.insert_order[i]]=GetInitialSearchR(data.insert_order[i]);
	}*/
	CalcRadius=true;
	BuildTree(data);
	// Check the integrity for all of the other points
	double searchR;
	for(int i=0;i<olength;++i)
	{
		point=data.insert_order[i];
		bool condition=true;
		//double maxR=FindMaxRadius(point);
		double maxR;
		searchR=GetInitialSearchR(data.insert_order[i],maxR);
		if(!BoundaryCandidate(maxR,point))
			continue;
		double tempR=0;
		while(condition)
		{
			PointWasAdded=false;
			for(int j=0;j<sidesNumber;++j)
			{
				if(bc->GetBoundaryType()==Rectengular||
					((bc->GetBoundaryType()==HalfPeriodic)&&(j==1||j==3)))
				{
					if(GetQueryPointReflective(qpoint,SidesToCheck[j],searchR,point))
					{
						queryPt[0]=qpoint.x;
						queryPt[1]=qpoint.y;
						AddSphere(tempR,searchR,queryPt,SidesToCheck[j],false,data);
					}
				}
				if(bc->GetBoundaryType()==Periodic||
					((bc->GetBoundaryType()==HalfPeriodic)&&(j!=1&&j!=3)))
				{
					if(GetQueryPointPeriodic(qpoint,SidesToCheck[j],searchR,point))
					{
						queryPt[0]=qpoint.x;
						queryPt[1]=qpoint.y;
						AddSphere(tempR,searchR,queryPt,SidesToCheck[j],false,data);
					}
				}
			}
			if(PointWasAdded)
				maxR=FindMaxRadius(point);
			if(maxR<=searchR)
			{
				condition=false;
			}
			else
			{
				tempR=searchR;
				searchR=tempR*1.4;
			}
		}
	}
	annDeallocPt(queryPt);
}

void Delaunay::GetOuterPoints2(vector<vector<int> > &OuterPoints)
{
	// Find the first point and triangle

	int cur_facet=Walk(olength);
	int startpoint=0;
	for(int i=0;i<3;++i)
	{
		if(f[cur_facet].get_vertice(i)<olength)
		{
			startpoint=f[cur_facet].get_vertice(i);
			break;
		}
	}
	vector<int> tetras;
	FindContainingTetras(cur_facet,startpoint,tetras);
	bool foundok=false;
	for(int i=0;i<(int)tetras.size();++i)
	{
		if(IsOuterFacet(f[tetras[i]].get_friend((FindPointInFacet(tetras[i],
			startpoint)+2)%3))&&!IsOuterFacet(tetras[i]))
		{
			cur_facet=tetras[i];
			foundok=true;
			break;
		}
	}
	if(!foundok)
		throw UniversalError("Error in locating first outer point in Delaunay.cpp");
	// Start adding the outer facets and allocate data
	if(OuterPoints.empty())
		OuterPoints.resize(4);
	vector<int> AllPoints;
	int n=(int)(2*sqrt(olength*1.0));
	AllPoints.reserve(4*n);
	if(olength<6)
	{
		for(int i=0;i<4;++i)
			for(int j=0;j<olength;++j)
				OuterPoints[i].push_back(j);
		return;
	}
	// Start the walking
	int point=startpoint;
	int nextpoint,nextfacet=cur_facet;
	AllPoints.push_back(point);
	nextfacet=FindNextOuterFacet(nextfacet,point,nextpoint);
	while(nextpoint!=startpoint)
	{
		AllPoints.push_back(nextpoint);
		point=nextpoint;
		nextfacet=FindNextOuterFacet(nextfacet,point,nextpoint);
	}
	// Find the containing tetras
	n=(int)AllPoints.size();
	vector<int> OuterTetras;
	OuterTetras.reserve(n*4);
	for(int i=0;i<n;++i)
	{
		int facet_index=Walk(AllPoints[i]);
		vector<int> facets;
		FindContainingTetras(facet_index,AllPoints[i],facets);
		for(int j=0;j<(int)facets.size();++j)
			OuterTetras.push_back(facets[j]);
	}
	// Sort and clean the vector
	sort(OuterTetras.begin(),OuterTetras.end());
	OuterTetras=unique(OuterTetras);
	// Find all of the points in the sorted facets
	AllPoints.clear();
	for(int i=0;i<(int)OuterTetras.size();++i)
		for(int j=0;j<3;++j)
			if(f[OuterTetras[i]].get_vertice(j)<olength)
				AllPoints.push_back(f[OuterTetras[i]].get_vertice(j));
	// Sort and clean the vector
	sort(AllPoints.begin(),AllPoints.end());
	AllPoints=unique(AllPoints);
	// Do Hilbert ordering
	vector<Vector2D> points_cor;
	n=(int)AllPoints.size();
	points_cor.resize(n);
	for(int i=0;i<n;++i)
		points_cor[i]=cor[AllPoints[i]];
	vector<int> order=HilbertOrder(points_cor,n,0);
	boost::array<Edge,4> edges;
	edges[1]=Edge(Vector2D(bc->GetGridBoundary(Right),bc->GetGridBoundary(Up)),
		Vector2D(bc->GetGridBoundary(Left),bc->GetGridBoundary(Up)),0,0);
	edges[2]=Edge(Vector2D(bc->GetGridBoundary(Left),bc->GetGridBoundary(Up)),
		Vector2D(bc->GetGridBoundary(Left),bc->GetGridBoundary(Down)),0,0);
	edges[3]=Edge(Vector2D(bc->GetGridBoundary(Left),bc->GetGridBoundary(Down)),
		Vector2D(bc->GetGridBoundary(Right),bc->GetGridBoundary(Down)),0,0);
	edges[0]=Edge(Vector2D(bc->GetGridBoundary(Right),bc->GetGridBoundary(Down)),
		Vector2D(bc->GetGridBoundary(Right),bc->GetGridBoundary(Up)),0,0);
	// Divide into sides
	for(int i=0;i<4;++i)
		OuterPoints[i].reserve((int)(n*0.26));
	for(int i=0;i<n;++i)
	{
		int min_index=0;
		double min_distance=DistanceToEdge(cor[AllPoints[order[i]]],edges[0]);
		for(int j=1;j<4;++j)
		{
			double temp=DistanceToEdge(cor[AllPoints[order[i]]],edges[j]);
			if(temp<min_distance)
			{
				min_distance=temp;
				min_index=j;
			}
		}
		OuterPoints[min_index].push_back(AllPoints[order[i]]);
	}
}

void Delaunay::GetOuterPoints(vector<vector<int> > &OuterPoints)
{
	// Find the first point and triangle

	int cur_facet=Walk(olength);
	int startpoint=0;
	for(int i=0;i<3;++i)
	{
		if(f[cur_facet].get_vertice(i)<olength)
		{
			startpoint=f[cur_facet].get_vertice(i);
			break;
		}
	}
	vector<int> tetras;
	FindContainingTetras(cur_facet,startpoint,tetras);
	bool foundok=false;
	for(int i=0;i<(int)tetras.size();++i)
	{
		if(IsOuterFacet(f[tetras[i]].get_friend((FindPointInFacet(tetras[i],
			startpoint)+2)%3))&&!IsOuterFacet(tetras[i]))
		{
			cur_facet=tetras[i];
			foundok=true;
			break;
		}
	}
	if(!foundok)
		throw UniversalError("Error in locating first outer point in Delaunay.cpp");
	// Allocate the space needed
	OuterPoints.clear();
	OuterPoints.resize(4);
	vector<int> AllPoints;
	int n=(int)(2*sqrt(olength*1.0));
	AllPoints.reserve(4*n);
	for(int i=0;i<4;++i)
		OuterPoints[i].reserve(n);
	if(olength<6)
	{
		for(int i=0;i<4;++i)
			for(int j=0;j<olength;++j)
				OuterPoints[i].push_back(j);
		return;
	}
	// Start the walking
	int point=startpoint;
	int nextpoint,nextfacet=cur_facet;
	AllPoints.push_back(point);
	nextfacet=FindNextOuterFacet(nextfacet,point,nextpoint);
	while(nextpoint!=startpoint)
	{
		AllPoints.push_back(nextpoint);
		point=nextpoint;
		nextfacet=FindNextOuterFacet(nextfacet,point,nextpoint);
	}
	// Divide the points into quadrants
	Vector2D direction;
	n=int(AllPoints.size());
	double angle;
	for(int i=0;i<n;++i)
	{
		direction=cor[AllPoints[(i+1)%n]]-cor[AllPoints[i]];
		angle=direction.y/direction.x;
		if(abs(angle)<0.2) 
		{
			// Horizontal
			if(direction.x<0)
				OuterPoints[1].push_back(AllPoints[i]);
			else
				OuterPoints[3].push_back(AllPoints[i]);
		}
		else
		{
			if(abs(angle)>0.8)
			{
				// Vertical
				if(direction.y>0)
					OuterPoints[0].push_back(AllPoints[i]);
				else
					OuterPoints[2].push_back(AllPoints[i]);
			}
			else
			{
				// We are at a corner probably
				if(angle>0)
				{
					if(direction.y>0)
					{
						OuterPoints[0].push_back(AllPoints[i]);
						OuterPoints[3].push_back(AllPoints[i]);
					}
					else
					{
						OuterPoints[1].push_back(AllPoints[i]);
						OuterPoints[2].push_back(AllPoints[i]);
					}
				}
				else
				{
					if(direction.y>0)
					{
						OuterPoints[0].push_back(AllPoints[i]);
						OuterPoints[1].push_back(AllPoints[i]);
					}
					else
					{
						OuterPoints[3].push_back(AllPoints[i]);
						OuterPoints[2].push_back(AllPoints[i]);
					}
				}
			}
		}
	}
}

int Delaunay::FindNextOuterFacet(int facet,int point,int &nextpoint)
{
	// Method gets an initial facet and point and returns
	// the next point on the outer edge and its facet
	int location=FindPointInFacet(facet,point);
	int next_facet=f[facet].get_friend(location);
	// Does facet have 3 points on edge?
	if(IsOuterFacet(next_facet))
	{
		nextpoint=f[facet].get_vertice((location+1)%3);
		return facet;
	}
	while(!IsOuterFacet(next_facet))
	{
		location=next_facet;
		next_facet=f[location].get_friend(FindPointInFacet(location,point));
	}
	nextpoint=f[location].get_vertice((FindPointInFacet(location,point)+1)%3);
	return location;
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

vector<int>& Delaunay::ChangeNewPointIndex(void)
{
	return NewPointIndex;
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
		return cor[f[Facet].get_vertice(vertice)].x;
	else 
		return cor[f[Facet].get_vertice(vertice)].y;
}

Vector2D Delaunay::get_point(int index) const
{
	return cor[index];
}

double Delaunay::get_cor(int index,int dim)
{
	if(dim==0) 
		return cor[index].x;
	else 
		return cor[index].y;
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

/* Reserved for future work
Vector2D Delaunay::CalCircleCenter(int facet)
{
	Vector2D B=cor[f[facet].get_vertice(1)]-cor[f[facet].get_vertice(0)];
	Vector2D C=cor[f[facet].get_vertice(2)]-cor[f[facet].get_vertice(0)];
	Vector2D res;
	double D_1=1/(2*(B.x*C.y-B.y*C.x));
	res.Set((C.y*(B.x*B.x+B.y*B.y)
		-B.y*(C.x*C.x+C.y*C.y))
		*D_1,(-C.x*(B.x*B.x+B.y*B.y)
		+B.x*(C.x*C.x+C.y*C.y))
		*D_1);
	return res;
}
*/

/* Reserved for future work
vector<double> Delaunay::GetNonLocalRadii(int point,vector<Vector2D> &centers)
{
	int start=Walk(point);
	vector<double> res;
	centers.clear();
	vector<int> facets;
	FindContainingTetras(start,point,facets);
	int n=facets.size();
	for(int i=0;i<n;++i)
	{
		//if(IsOuterFacet(facets[i]))
		{
			centers.push_back(CalCircleCenter(facets[i]));
			res.push_back(radius[facets[i]]);
		}
	}
	return res;
}
*/

/* Reserved for future work
void Delaunay::AddBigTirangle(int point)
{
	vector<double> radii;
	vector<Vector2D> centers;
	ANNpoint queryPt;
	queryPt=annAllocPt(2);
	double searchR;
	//bool condition=true;
	//condition=false;
	radii=GetNonLocalRadii(point,centers);
	int n=radii.size();
	for(int j=0;j<n;++j)
	{
		queryPt[0]=centers[j].x;
		queryPt[1]=centers[j].y;
		searchR=radii[j];
		ANNidxboost::array nnIdx; // near neighbor indices
		ANNdistboost::array dists; // near neighbor distances
		int NumberInSphere=tree->annkFRSearch(queryPt,searchR*searchR,0);
		if(NumberInSphere>5)
			continue;
		nnIdx = new ANNidx[NumberInSphere]; // allocate near neigh indices
		dists = new ANNdist[NumberInSphere]; // allocate near neighbor dists
		// Get the points
		tree->annkFRSearch(queryPt,searchR*searchR,NumberInSphere,nnIdx,dists,0);
		for(int i=0;i<NumberInSphere;++i)
		{
			cor.push_back(totalCor[nnIdx[i]]);
			int facetLoc=WalkRobust(cor.size()-1);
			if(facetLoc!=-1)
			{
				add_point(cor.size()-1);
				++length;
				NewPointIndex.push_back(nnIdx[i]);
			}
			else
				cor.pop_back();
		}
		delete [] nnIdx;
		delete [] dists;
	}
	annDeallocPt(queryPt);
}
*/


/* Reserve for future work
void Delaunay::build_limited(vector<Vector2D>const& vp,vector<int> 
	_pointIndex,OuterBoundary const *_bc)
{// copy the boundary conditions
	bc=_bc;
	lastFacet=0;
	double d_temp;
	double maxV=(bc->GetGridBoundary(Right)-bc->GetGridBoundary(Left)),
		maxH=(bc->GetGridBoundary(Up)-bc->GetGridBoundary(Down));
	BIG=pow(max(maxV,maxH),2);
	eps=max(abs(bc->GetGridBoundary(Up)),abs(bc->GetGridBoundary(Down)));
	d_temp=max(abs(bc->GetGridBoundary(Right)),abs(bc->GetGridBoundary(Left)));
	eps=max(eps,d_temp);
	eps=eps*1E-12;
	length=_pointIndex.size()+3;
	int len=length-3;
	olength=len;
	f.clear();
	cor.clear();
	f.reserve(2*length+1);
	cor.reserve(length);
	last_loc=100*length;
	for(int i=0;i<(int)_pointIndex.size();i++)
	{
		cor.push_back(vp[_pointIndex[i]]);
	}
	// Check point input
	CheckInput();

	insert_order=HilbertOrder(cor,cor.size());	

	totalCor.reserve(vp.size()+3);
	for(int i=0;i<(int)vp.size();i++)
	{
		totalCor.push_back(vp[i]);
	}


	// add the 3 extreme points
	Vector2D p_temp;
	p_temp.set_x(bc->GetGridBoundary(Left)-5*bc->GetGridBoundary(Right));
	p_temp.set_y(bc->GetGridBoundary(Down)-5*bc->GetGridBoundary(Up));
	cor.push_back(p_temp);
	p_temp.set_x(5*bc->GetGridBoundary(Right));
	p_temp.set_y(bc->GetGridBoundary(Down)-5*bc->GetGridBoundary(Up));
	cor.push_back(p_temp);
	p_temp.set_x((bc->GetGridBoundary(Left)+bc->GetGridBoundary(Right))/2.0);
	p_temp.set_y(5*bc->GetGridBoundary(Up));
	cor.push_back(p_temp);
	// Create the big triangle, and assign friends
	facet f_temp;
	f.push_back(f_temp);
	f[0].set_vertice(len,0);
	f[0].set_vertice(len+1,1);
	f[0].set_vertice(len+2,2);
	for(int i=0;i<3;i++)
		f[0].set_friend(last_loc,i);
	location_pointer=0;
	// add the points
	for(int i=0;i<length-3;i++)
	{// Need to fix here lack of hilber order
		add_point(insert_order[i]);
	}
	BuildBoundaryLimited();
}
*/

/* Reserved for future work
void Delaunay::BuildBoundaryLimited(void)
{
	// Build kd-tree
	treePoints=annAllocPts(totalCor.size(),2);
	for(int i=0;i<(int)totalCor.size();i++)
	{
		treePoints[i][0]=totalCor[i].x;
		treePoints[i][1]=totalCor[i].y;
	}
	tree=new ANNkd_tree(treePoints,totalCor.size()-3,2,1,ANN_KD_SUGGEST);	
	NewPointIndex.clear();
	NewPointIndex.reserve(8*(int)sqrt(cor.size()*1.2));
	if(lastRadius.capacity()!=cor.size())
	{
		lastRadius.reserve(cor.size());
		for(int i=0;i<(int)cor.size();i++)
			lastRadius.push_back(0.3*BIG/sqrt((double)(totalCor.size())));
	}

	DuplicatePointsLimited();

	//annDeallocPts(treePoints);
	//delete tree;
	annClose();	
}
*/

/*
Reserved for future work
void Delaunay::DuplicatePointsLimited(void)
{
	int sidesNumber=0,point=0;
	if(bc->GetBoundaryType()==Rectengular) // How many duplicates do I need to check
		sidesNumber=4;
	else
		sidesNumber=8;
	vector<Sides> SidesToCheck;
	SidesToCheck.resize(sidesNumber);
	for(int i=0;i<sidesNumber;i++)
		SidesToCheck[i]=(Sides)i;

	double searchR,tempR;
	bool condition;
	Vector2D qpoint;
	ANNpoint queryPt;
	queryPt=annAllocPt(2);
	for(int i=0;i<(int)insert_order.size();i++)
	{
		point=insert_order[i];
		searchR=lastRadius[point];
		tempR=0;
		condition=true;
		while(condition)
		{
			queryPt[0]=cor[point].x;
			queryPt[1]=cor[point].y;
			int NumberInSphere=tree->annkFRSearch(queryPt,searchR*searchR,0);
			if(NumberInSphere>25)
			{
				AddBigTirangle(point);
				break;
			}
			ANNidxboost::array nnIdx; // near neighbor indices
			ANNdistboost::array dists; // near neighbor distances
			nnIdx = new ANNidx[NumberInSphere]; // allocate near neigh indices
			dists = new ANNdist[NumberInSphere]; // allocate near neighbor dists
			// Get the points
			tree->annkFRSearch(queryPt,searchR*searchR,NumberInSphere,nnIdx,dists,0);
			for(int i=0;i<NumberInSphere;i++)
			{
				if(dists[i]>(tempR*tempR))
				{
					cor.push_back(totalCor[nnIdx[i]]);
					int facetLoc=WalkRobust(cor.size()-1);
					if(facetLoc!=-1)
					{
						add_point(cor.size()-1);
						++length;
						NewPointIndex.push_back(nnIdx[i]);
					}
					else
						cor.pop_back();
				}
			}
			double maxR=FindMaxRadius(point);
			if(maxR<=searchR)
			{
				condition=false;
			}
			else
			{
				tempR=searchR;
				searchR=tempR*1.2;
			}
			delete [] nnIdx;
			delete [] dists;
		}
	}
	// Outer boundaries
	annDeallocPts(treePoints);
	delete tree;
	// Build kd-tree
	treePoints=annAllocPts(olength,2);
	for(int i=0;i<olength;i++)
	{
		treePoints[i][0]=cor[i].x;
		treePoints[i][1]=cor[i].y;
	}
	tree=new ANNkd_tree(treePoints,olength,2,1,ANN_KD_SUGGEST);
	int sidenum;
	for(int i=0;i<(int)insert_order.size();i++)
	{
		tempR=0;
		condition=true;
		searchR=lastRadius[point];
		sidenum=sidesNumber;
		while(condition)
		{
			for(int j=0;j<sidenum;j++)
			{
				if(bc->GetBoundaryType()==Rectengular)
				{
					if(GetQueryPointReflective(qpoint,SidesToCheck[j],searchR,point))
					{
						queryPt[0]=qpoint.x;
						queryPt[1]=qpoint.y;
						int ndup=AddSphere(tempR,searchR,queryPt,SidesToCheck[j],false);
						if(ndup>25)
						{
							SidesToCheck.erase(SidesToCheck.begin()+j);
							--j;
							--sidenum;
						}
					}
				}
				if(bc->GetBoundaryType()==Periodic)
				{
					if(GetQueryPointPeriodic(qpoint,SidesToCheck[j],searchR,point))
					{
						queryPt[0]=qpoint.x;
						queryPt[1]=qpoint.y;
						int ndup=AddSphere(tempR,searchR,queryPt,SidesToCheck[j],false);
						if(ndup>25)
						{
							SidesToCheck.erase(SidesToCheck.begin()+j);
							--j;
							--sidenum;
						}
					}
				}
			}
			double maxR=FindMaxRadius(point);
			if(maxR<=searchR)
			{
				condition=false;
				lastRadius[point]=maxR*1.2;
			}
			else
			{
				tempR=searchR;
				searchR=tempR*1.2;
			}
		}
	}
	annDeallocPt(queryPt);
	annDeallocPts(treePoints);
	delete tree;
	annClose();	
}
*/
