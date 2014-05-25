#include "VoronoiMesh.hpp"
#include <iostream>
#include <cmath>

namespace {
  UniversalError negative_volume_ratio(VoronoiMesh const& V,int index,
				       vector<int> const& ToRemove)
  {
    UniversalError eo("Negative volume ratio in cell refine");
    vector<int> bad_edgesindex=V.GetCellEdges(index);
    for(int jj=0;jj<(int)bad_edgesindex.size();++jj)
      {
	Edge e=V.GetEdge(bad_edgesindex[jj]);
	eo.AddEntry("Edge X cor",e.get_x(0));
	eo.AddEntry("Edge Y cor",e.get_y(0));
	eo.AddEntry("Edge X cor",e.get_x(1));
	eo.AddEntry("Edge Y cor",e.get_y(1));
	eo.AddEntry("Edge length",e.GetLength());
	eo.AddEntry("Edge neighbor 0",e.GetNeighbor(0));
	eo.AddEntry("Edge neighbor 1",e.GetNeighbor(1));
      }
    for(int i=0;i<(int)ToRemove.size();++i)
      {
	eo.AddEntry("ToRemove Index",ToRemove[i]);
	eo.AddEntry("ToRemove x cor",V.GetMeshPoint(ToRemove[i]).x);
	eo.AddEntry("ToRemove y cor",V.GetMeshPoint(ToRemove[i]).y);
      }
    return eo;
  }
}
  
Vector2D VoronoiMesh::CalcFaceVelocity(Vector2D wl, Vector2D wr,Vector2D rL, Vector2D rR,
	Vector2D f)const
{
	Vector2D wprime = ScalarProd(wl-wr,f-(rR+rL)/2)*	(rR-rL)/pow(abs(rR-rL),2);
	return 0.5*(wl+wr) + wprime;
}

vector<Vector2D> VoronoiMesh::calc_edge_velocities(HydroBoundaryConditions const* hbc,
	vector<Vector2D> const& point_velocities,double time)const
{
	vector<Vector2D> facevelocity;
	facevelocity.resize(edges.size());
	for(int i = 0; i < (int)edges.size(); ++i)
	{
		if(hbc->IsBoundary(edges[i],this))
		{
			// Boundary
			facevelocity[i] = hbc->CalcEdgeVelocity(this,point_velocities,edges[i],
				time);
		}
		else
		{
			// Bulk
			facevelocity[i] = CalcFaceVelocity(
				point_velocities[edges[i].GetNeighbor(0)],
				point_velocities[edges[i].GetNeighbor(1)],
				GetMeshPoint(edges[i].GetNeighbor(0)),
				GetMeshPoint(edges[i].GetNeighbor(1)),
				0.5*(edges[i].GetVertex(0)+edges[i].GetVertex(1)));
		}
	}
	return facevelocity;
}


bool VoronoiMesh::NearBoundary(int index) const
{
  int n=int(mesh_vertices[index].size());
	int n1,n0;
	for(int i=0;i<n;++i)
	{
		n0=edges[mesh_vertices[index][i]].GetNeighbor(0);
		n1=edges[mesh_vertices[index][i]].GetNeighbor(1);
		if(n0<0||n1<0||n0>n||n1>n)
			return true;
	}
	return false;
}

VoronoiMesh::~VoronoiMesh(void)
{
	//	cout<<"Entered v destruct"<<endl;
	delete Tri;
	edges.clear();
	mesh_vertices.clear();
}

int VoronoiMesh::GetOriginalIndex(int point) const
{
	return Tri->GetOriginalIndex(point);
}


VoronoiMesh::VoronoiMesh(): 
  bc(0),
  eps(1e-8),
  edges(vector<Edge>()),
  added(vector<vector<int> >()),
  CM(vector<Vector2D> ()),
  mesh_vertices(vector<vector<int> >()),
  Tri(0) {}

VoronoiMesh::VoronoiMesh(VoronoiMesh const& other):
  bc(other.bc),
  eps(other.eps),
  edges(other.edges),
  added(other.added),
  CM(other.CM),
  mesh_vertices(other.mesh_vertices),
  Tri(new Delaunay(*other.Tri)) {}

void VoronoiMesh::build_v()
{
	added.clear();
	added.resize(Tri->GetTotalLength());
	vector<int>::iterator it;
	Vector2D center,center_temp;
	int j;
	facet *to_check;
	Edge edge_temp;
	Vector2D p_temp;
	mesh_vertices.clear();
	mesh_vertices.resize(Tri->get_length());
	edges.reserve((int)(Tri->get_length()*3.5));
	int N=Tri->GetOriginalLength();
	for(int i=0;i<N;++i)
	{
		mesh_vertices[i].reserve(7);
		added[i].reserve(7);
	}

	for(int i=0;i<Tri->get_num_facet();i++)
	{
	  //		int NNN=edges.size();
		center=get_center(i);
		to_check=Tri->get_facet(i);
		for(j=0;j<3;j++)
		{
			if(to_check->get_friend(j)==Tri->get_last_loc())
				continue;
			//			double xx=center.x;
			center_temp=get_center(to_check->get_friend(j));
			{
				edge_temp.set_x(0,center.get_x());
				edge_temp.set_x(1,center_temp.get_x());
				edge_temp.set_y(0,center.get_y());
				edge_temp.set_y(1,center_temp.get_y());
				edge_temp.set_friend(0,to_check->get_vertice(j));
				edge_temp.set_friend(1,to_check->get_vertice((j+1)%3));
				if(legal_edge(&edge_temp))
				{
					// I added a change here, if edge has zero length I don't add it.
					if(edge_temp.GetLength()>eps*max(Tri->GetFacetRadius(edge_temp.GetNeighbor(0)),
						Tri->GetFacetRadius(edge_temp.GetNeighbor(1))))
					{
						it=find(added[edge_temp.GetNeighbor(1)].begin(),
							added[edge_temp.GetNeighbor(1)].end(),
							edge_temp.GetNeighbor(0));
						if(it==added[edge_temp.GetNeighbor(1)].end())
						{
							added[edge_temp.GetNeighbor(0)].push_back(edge_temp.GetNeighbor(1));
							if(bc->GetBoundaryType()==Periodic||
								(bc->GetBoundaryType()==HalfPeriodic&&
								(!bc->AreWeReflective(edge_temp))))
							{
								if(edge_temp.GetNeighbor(0)<=Tri->GetOriginalLength())
								  mesh_vertices[edge_temp.GetNeighbor(0)].push_back((int)edges.size());
								if(edge_temp.GetNeighbor(1)<=Tri->GetOriginalLength())
								  mesh_vertices[edge_temp.GetNeighbor(1)].push_back((int)edges.size());
								edges.push_back(edge_temp);
							}
							if(bc->GetBoundaryType()==Rectengular||
								(bc->GetBoundaryType()==HalfPeriodic&&
								(bc->AreWeReflective(edge_temp))))
							{
								if(edge_temp.GetNeighbor(0)>Tri->GetOriginalLength())
									edge_temp.set_friend(0,-1);
								else
								  mesh_vertices[edge_temp.GetNeighbor(0)].push_back((int)edges.size());
								if(edge_temp.GetNeighbor(1)>Tri->GetOriginalLength())
									edge_temp.set_friend(1,-1);
								else
								  mesh_vertices[edge_temp.GetNeighbor(1)].push_back((int)edges.size());
								edges.push_back(edge_temp);
							}
						}
					}
				}
			}
		}
	}
}

void VoronoiMesh::Initialise(vector<Vector2D>const& pv,OuterBoundary const* _bc)
{
	//Build the delaunay
	delete Tri;
	Tri=new Delaunay;
	Tri->build_delaunay(pv,_bc);
	// copy the boundary conditions
	bc=_bc;
	eps=1E-8;
	edges.clear();
	build_v();
	vector<vector<int> > ().swap(added);
	int n=GetPointNo();
	CM.resize(n);
	for(int i=0;i<n;++i)
		CM[i]=CalcCellCM(i);
}


Vector2D VoronoiMesh::get_center(int facet)
	// This function calculates the center of the circle surrounding the Facet
{
	Vector2D center;
	double x1=Tri->get_facet_coordinate(facet,0,0);
	double y1=Tri->get_facet_coordinate(facet,0,1);
	double x2=Tri->get_facet_coordinate(facet,1,0);
	double y2=Tri->get_facet_coordinate(facet,1,1);
	double x3=Tri->get_facet_coordinate(facet,2,0);
	double y3=Tri->get_facet_coordinate(facet,2,1);
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


void VoronoiMesh::output(string filename)
{	
  int temp2=int(edges.size());
	fstream myFile(filename.c_str(),ios::out | ios::binary);
	myFile.write ((char*)&temp2, sizeof (int));
	double temp;
	for(int i=0;i<temp2;i++)
	{
		//if(RemovedEdges.size()==0||!binary_search(RemovedEdges.begin(),
		//RemovedEdges.end(),i))
		{
			temp=edges[i].get_x(0);
			myFile.write ((char*)&temp, sizeof (double));
			temp=edges[i].get_x(1);
			myFile.write ((char*)&temp, sizeof (double));
		}
	}
	for(int i=0;i<temp2;i++)
	{
		//if(RemovedEdges.size()==0||!binary_search(RemovedEdges.begin(),
		//RemovedEdges.end(),i))
		{
			temp=edges[i].get_y(0);
			myFile.write ((char*)&temp, sizeof (double));
			temp=edges[i].get_y(1);
			myFile.write ((char*)&temp, sizeof (double));
		}
	}
	int j;
	for(int i=0;i<temp2;i++)
	{
		//if(RemovedEdges.size()==0||!binary_search(RemovedEdges.begin(),
		//RemovedEdges.end(),i))
		{
			j=edges[i].GetNeighbor(0);
			myFile.write ((char*)&j, sizeof (int));
			j=edges[i].GetNeighbor(1);
			myFile.write ((char*)&j, sizeof (int));
		}
	}
	temp2=Tri->get_length();
	myFile.write ((char*)&temp2, sizeof (int));
	for(int i=0;i<temp2;i++)
	{
		temp=Tri->get_cor(i,0);
		myFile.write ((char*)&temp, sizeof (double));
		temp=Tri->get_cor(i,1);
		myFile.write ((char*)&temp, sizeof (double));
	}
	vector<int> indexes;
	int length;
	for(int i=0;i<temp2;i++)
	{
		indexes=GetCellEdges(i);
		length=int(indexes.size());
		myFile.write((char*)&length,sizeof(int));
		for(int j2=0;j2<(int)indexes.size();j2++)
			myFile.write((char*)&indexes[j2],sizeof(int));
	}
	myFile.close();
	}


bool VoronoiMesh::legal_edge(Edge *e) //checks if both ends of the edge are outside the grid and that the edge doesn't cross the grid
{
	if((e->GetNeighbor(0)<Tri->get_length())||
		(e->GetNeighbor(1)<Tri->get_length()))
		return true;
	else
		return false;
}

double VoronoiMesh::GetWidth(int index)const
{
  return sqrt(GetVolume(index)/M_PI);
}

vector<int> const& VoronoiMesh::GetCellEdges(int index) const
{
  return mesh_vertices[index];
}

double VoronoiMesh::GetVolume(int index) const
{
	double x1,x2;
	double y1,y2;
	Vector2D center=Tri->get_point(index);
	double area=0;
	for (size_t i=0;i<mesh_vertices[index].size();++i)
	{
		x1=edges[mesh_vertices[index][i]].get_x(0)-center.get_x();
		x2=edges[mesh_vertices[index][i]].get_x(1)-center.get_x();
		y1=edges[mesh_vertices[index][i]].get_y(0)-center.get_y();
		y2=edges[mesh_vertices[index][i]].get_y(1)-center.get_y();
		area+=abs(x2*y1-y2*x1)*0.5;
	}
	return area;
}

Vector2D VoronoiMesh::CalcCellCM(int index) const
{
	double x1,x2;
	double y1,y2;
	double xc=0,yc=0;
	Vector2D center=Tri->get_point(index);
	double area=0,area_temp;
	for (int i=0;i<(int)mesh_vertices[index].size();i++)
	{
		int temp=mesh_vertices[index][i];
		x1=edges[temp].get_x(0)-center.get_x();
		x2=edges[mesh_vertices[index][i]].get_x(1)-center.get_x();
		y1=edges[mesh_vertices[index][i]].get_y(0)-center.get_y();
		y2=edges[mesh_vertices[index][i]].get_y(1)-center.get_y();
		area_temp=abs(x2*y1-y2*x1)*0.5;
		area+=area_temp;
		xc+=area_temp*(center.get_x()+edges[temp].get_x(0)+edges[temp].get_x(1))/3;
		yc+=area_temp*(center.get_y()+edges[temp].get_y(0)+edges[temp].get_y(1))/3;
	}
	area=1/area;
	Vector2D p(xc*area,yc*area);
	return p;
}

vector<Vector2D>& VoronoiMesh::GetMeshPoints(void)
{
  return Tri->GetMeshPoints();
}

void VoronoiMesh::Update(vector<Vector2D> const& points)
{
	// Clean_up last step
	edges.clear();	
	Tri->update(points);
	build_v();
	vector<vector<int> > ().swap(added);
	int n=GetPointNo();
	CM.resize(n);
	for(int i=0;i<n;++i)
		CM[i]=CalcCellCM(i);
}

vector<int> VoronoiMesh::GetNeighbors(int index)const
{
  int n=(int)mesh_vertices[index].size();
	int other;
	vector<int> res;
	res.reserve(n);
	for(int i=0;i<n;++i)
	{
		if((other=GetOriginalIndex(
			edges[mesh_vertices[index][i]].GetNeighbor(0)))!=index)
		{
			//if(other!=-1)
			res.push_back(other);
		}
		else
		{
			other=GetOriginalIndex(edges[mesh_vertices[index][i]]
			.GetNeighbor(1));
			//if(other!=-1)
			res.push_back(other);
		}	
	}
	return res;
}


Tessellation* VoronoiMesh::clone(void)const
{
	return new VoronoiMesh(*this);
}

void VoronoiMesh::RemoveCells(vector<int> &ToRemove,vector<vector<int> > &VolIndex,
	vector<vector<double> > &Volratio)
{
	Remove_Cells(*this,ToRemove,VolIndex,Volratio);
}

void Remove_Cells(VoronoiMesh &V,vector<int> &ToRemove,
	vector<vector<int> > &VolIndex,vector<vector<double> > &Volratio)
{
	VolIndex.clear();
	Volratio.clear();
	size_t n=ToRemove.size();
	VolIndex.resize(n);
	Volratio.resize(n);
	vector<int> RemovedEdges;
	for(size_t i=0;i<n;++i)
	{
		double vol_inverse=1/V.GetVolume(ToRemove[i]);
		// Get neighboring points
		vector<int> neigh=V.GetNeighbors(ToRemove[i]);
		vector<Vector2D> neighcor;
		//Remove The boundary points
		vector<int> real_neigh; //Keeps track of the real points I copied
		vector<int> period_points;// Keeps track of the points which were reflected
		real_neigh.reserve(10);
		period_points.reserve(10);
		for(size_t j=0;j<neigh.size();++j)
		{
			if(neigh[j]>-1)
			{
				real_neigh.push_back(V.Tri->GetOriginalIndex(
					neigh[j]));	
				period_points.push_back(neigh[j]);
			}
		}
		vector<int> AllPoints=real_neigh;
		for(size_t j=0;j<real_neigh.size();++j)
		{
			vector<int> temp_neigh=V.GetNeighbors(real_neigh[j]);
			for(size_t jj=0;jj<temp_neigh.size();++jj)
				if(temp_neigh[jj]>-1&&temp_neigh[jj]!=ToRemove[i])
					AllPoints.push_back(temp_neigh[jj]);
		}
		// save old volume
		vector<double> old_vol(real_neigh.size());
		sort(real_neigh.begin(),real_neigh.end());
		for(size_t j=0;j<real_neigh.size();++j)
			old_vol[j]=V.GetVolume(real_neigh[j]);
		sort(AllPoints.begin(),AllPoints.end());
		AllPoints=unique(AllPoints);
		AllPoints=RemoveList(AllPoints,real_neigh);
		AllPoints.insert(AllPoints.begin(),real_neigh.begin(),real_neigh.end());
		neighcor.reserve(AllPoints.size());
		for(size_t j=0;j<AllPoints.size();++j)
			neighcor.push_back(V.Tri->get_point(AllPoints[j]));
		// build new local tessellation
		VoronoiMesh Vlocal;
		Vlocal.Initialise(neighcor,V.bc);
		// save the volume change
		for(size_t j=0;j<real_neigh.size();++j)
		{
		  old_vol[(int)j]=(Vlocal.GetVolume((int)j)-old_vol[(int)j])*vol_inverse;
			if(old_vol[j]<(-1e-5))
			{
				cout<<"Negative volume ratio is "<<old_vol[j]<<endl;
				throw negative_volume_ratio(V,real_neigh[j],ToRemove);
			}
		}
		// Check we did ok
		double vol_sum=0;
		for(size_t j=0;j<real_neigh.size();++j)
		{
			vol_sum+=old_vol[j];
		}
		if(abs(1-vol_sum)>1e-4)
		{
			UniversalError eo("Bad total volume in Voronoi Remove_Cells");
			eo.AddEntry("Total Volume",vol_sum);
			eo.AddEntry("Removed cell",ToRemove[i]);
			throw eo;
		}
		Volratio[i]=old_vol;
		VolIndex[i]=real_neigh;
		// Fix the old Voronoi////////////////
		int Nedges=(int)Vlocal.edges.size();
		// Keep only the relevant edges
		vector<Edge> NewEdges;
		int N=(int)real_neigh.size();
		for(int j=0;j<Nedges;++j)
		{
			int temp0=Vlocal.edges[j].GetNeighbor(0);
			int temp1=Vlocal.edges[j].GetNeighbor(1);
			if(temp0<N)
			{
				if(temp1<N||temp1<0||Vlocal.GetOriginalIndex(temp1)<N)
				{
					NewEdges.push_back(Vlocal.edges[j]);
					continue;
				}
			}
			if(temp1<N)
				if(temp0<0||Vlocal.GetOriginalIndex(temp0)<N)
					NewEdges.push_back(Vlocal.edges[j]);
		}
		// Fix the new edges to correspond to the correct points in cor, use real_neigh
		for(size_t j=0;j<NewEdges.size();++j)
		{
			int temp=NewEdges[j].GetNeighbor(0);
			if(temp>N)
				NewEdges[j].set_friend(0,period_points[
					Vlocal.GetOriginalIndex(temp)]);
			else
				if(temp>-1)
					NewEdges[j].set_friend(0,real_neigh[temp]);
			temp=NewEdges[j].GetNeighbor(1);
			if(temp>N)
				NewEdges[j].set_friend(1,period_points[
					Vlocal.GetOriginalIndex(temp)]);
			else
				if(temp>-1)
					NewEdges[j].set_friend(1,real_neigh[temp]);
		}
		//Remove bad edge reference from mesh vertices and outer edges
		vector<int> oldedges=V.mesh_vertices[ToRemove[i]];
		for(int j=0;j<N;++j)
		{
			int temp1=real_neigh[j];
			int NN=(int)V.mesh_vertices[temp1].size();
			for(int jj=0;jj<NN;++jj)
			{
				Edge etemp=V.edges[V.mesh_vertices[temp1][jj]];
				if(etemp.GetNeighbor(0)==-1||etemp.GetNeighbor(0)==
					ToRemove[i])
				{
					if(etemp.GetNeighbor(0)<0)
						oldedges.push_back(V.mesh_vertices[temp1][jj]);
					V.mesh_vertices[temp1].erase(V.mesh_vertices[temp1].begin()
						+jj);
					--NN;
					--jj;
					continue;
				}
				if(etemp.GetNeighbor(1)==-1||etemp.GetNeighbor(1)==
					ToRemove[i])
				{
					if(etemp.GetNeighbor(1)<0)
						oldedges.push_back(V.mesh_vertices[temp1][jj]);
					V.mesh_vertices[temp1].erase(V.mesh_vertices[temp1].begin()
						+jj);
					--NN;
					--jj;
					continue;
				}
			}
		}

		// copy the new edges	
		Nedges=(int)NewEdges.size();	
		// Fix the mesh vertices to point to the correct edges
		int temp,other;
		for(int j=0;j<Nedges;++j)
		{
			for(int k=0;k<2;++k)
			{
				temp=V.GetOriginalIndex(NewEdges[j].GetNeighbor(k));
				if(temp==-1)
					continue;
				other=NewEdges[j].GetNeighbor((k+1)%2);
				size_t jj=0;			
				if(other!=-1)
				{
					for(jj=0;jj<V.mesh_vertices[temp].size();++jj)
					{
						if(V.edges[V.mesh_vertices[temp][jj]].GetNeighbor(0)==
							other)
						{
							V.edges[V.mesh_vertices[temp][jj]]=NewEdges[j];
							break;
						}
						if(V.edges[V.mesh_vertices[temp][jj]].GetNeighbor(1)==
							other)
						{
							V.edges[V.mesh_vertices[temp][jj]]=NewEdges[j];
							break;
						}
					}
				}
				else
				{
				  V.mesh_vertices[temp].push_back((int)V.edges.size());
					V.edges.push_back(NewEdges[j]);
				}
				// We made New neighbors
				if(jj==V.mesh_vertices[temp].size())
				{
				  V.mesh_vertices[temp].push_back((int)V.edges.size()-k);
					if(k==0)
						V.edges.push_back(NewEdges[j]);
				}
			}
		}
		RemovedEdges.insert(RemovedEdges.end(),oldedges.begin(),oldedges.end());
	}
	// Done with all of the points to remove now start fixing
	sort(RemovedEdges.begin(),RemovedEdges.end());
	RemovedEdges=unique(RemovedEdges);
	sort(ToRemove.begin(),ToRemove.end());

	//Fix edges
	size_t Nedges=V.edges.size();
	int toReduce;
	int temp;
	for(size_t i=0;i<Nedges;++i)
	{
		temp=V.edges[i].GetNeighbor(0);
		if(temp>-1)
		{
		  toReduce=int(lower_bound(ToRemove.begin(),ToRemove.end(),temp)-
			       ToRemove.begin());
			V.edges[i].set_friend(0,temp-toReduce);
		}
		temp=V.edges[i].GetNeighbor(1);
		if(temp>-1)
		{
		  toReduce=int(lower_bound(ToRemove.begin(),ToRemove.end(),temp)-
			       ToRemove.begin());
			V.edges[i].set_friend(1,temp-toReduce);
		}
	}
	// Remove bad edges
	RemoveVector(V.edges,RemovedEdges);
	// Fix mesh_vertices
	RemoveVector(V.mesh_vertices,ToRemove);	
	size_t N=V.mesh_vertices.size();
	// Fix edge number in mesh_vertices
	for(size_t i=0;i<N;++i)
	{
		for(size_t j=0;j<V.mesh_vertices[i].size();++j)
		{
		  toReduce=int(lower_bound(RemovedEdges.begin(),RemovedEdges.end(),
					   V.mesh_vertices[i][j])-RemovedEdges.begin());
			V.mesh_vertices[i][j]-=toReduce;
		}
	}
	// Fix cor in Tri
	RemoveVector(V.Tri->ChangeCor(),ToRemove);
	// Fix CM
	V.CM.resize(V.CM.size()-ToRemove.size());
	int NN=(int)V.CM.size();
	for(int i=0;i<NN;++i)
		V.CM[i]=V.CalcCellCM(i);
	// Fix point numer in Tri
	V.Tri->ChangeOlength(V.Tri->get_length()-(int)n);
	V.Tri->Changelength(V.Tri->get_length()-(int)n);
	// Fix NewIndex in Tri
	vector<int>& NewPoint=V.Tri->ChangeNewPointIndex();
	N=NewPoint.size();
	for(size_t i=0;i<N;++i)
	{
		temp=NewPoint[i];
		toReduce=int(lower_bound(ToRemove.begin(),ToRemove.end(),temp)-
			     ToRemove.begin());
		NewPoint[i]-=toReduce;
	}
	return;
	// Reset Tree if self gravity is needed
}


void VoronoiMesh::RefineCells(vector<int> const& ToRefine,
	vector<Vector2D> const& directions,double alpha)
{
	Refine_Cells(*this,ToRefine,alpha,directions,bc->GetBoundaryType()==Periodic);
	return;
}

void Refine_Cells(VoronoiMesh &V,vector<int> const& ToRefine,double alpha,
	vector<Vector2D> const& directions,bool PeriodicBoundary)
{
  int n=(int)ToRefine.size();
	if(n==0)
		return;
	int Npoints=V.GetPointNo();
	vector<Vector2D>& cor=V.Tri->ChangeCor();
	vector<Vector2D> cortemp(cor.size()-Npoints);
	// copy the ghost points into cortemp
	copy(cor.begin()+Npoints,cor.end(),cortemp.begin());
	int N=(int)cor.size();
	// Expand the cor vector to include the new points
	cor.resize(N+n);
	copy(cortemp.begin(),cortemp.end(),cor.begin()+Npoints+n);
	cortemp.clear();
	// Fix the boundary point refrences if needed
	if(PeriodicBoundary)
	{
		for(int i=0;i<(int)V.edges.size();++i)
		{
			int temp_neigh=V.edges[i].GetNeighbor(0);
			if(temp_neigh>Npoints)
			{
				V.edges[i].set_friend(0,temp_neigh+n);
				continue;
			}
			temp_neigh=V.edges[i].GetNeighbor(1);
			if(temp_neigh>Npoints)
				V.edges[i].set_friend(1,temp_neigh+n);
		}
	}
	// Update the lengths
	V.Tri->Changelength(Npoints+n);
	V.Tri->ChangeOlength(Npoints+n);
	// reserve space for mesh_vertices
	V.mesh_vertices.resize(Npoints+n);
	// Refine the points
	for(int i=0;i<n;++i)
	{
		// Get the edges
		vector<int> edge_index=V.GetCellEdges(ToRefine[i]);
		int nedges=(int)edge_index.size();
		vector<Edge> edges(nedges);
		for(int j=0;j<nedges;++j)
			edges[j]=V.GetEdge(edge_index[j]);
		// Make the new point by finding the closest edge
		Vector2D NewPoint=V.GetMeshPoint(ToRefine[i]);
		double R=V.GetWidth(ToRefine[i]);
		Vector2D normal;
		Vector2D slope;
		// Do we have a specified direction or do we do line conneting the CM or nearest edge
		if(!directions.empty())
		{
			slope=directions[i]/abs(directions[i]);
			normal.Set(slope.y,-slope.x);
		}
		else
			slope=FindBestSplit(&V,ToRefine[i],edges,R,normal);
		NewPoint+=alpha*R*slope;
		cor[Npoints+i]=NewPoint;
		Vector2D v(V.GetMeshPoint(ToRefine[i]));
		// Split edges and update neighbors
		vector<int> old_ref;
		vector<int> new_ref;
		// Calculate the splitting edge
		Edge splitedge;
		splitedge.set_friend(1,ToRefine[i]);
		splitedge.set_friend(0,Npoints+i);
		splitedge.set_x(0,0.5*(NewPoint.x+V.GetMeshPoint(ToRefine[i]).x));
		splitedge.set_y(0,0.5*(NewPoint.y+V.GetMeshPoint(ToRefine[i]).y));
		splitedge.set_x(1,splitedge.get_x(0)+normal.x*R*20);
		splitedge.set_y(1,splitedge.get_y(0)+normal.y*R*20);
		Vector2D intersect;
		// Find other intersecting segment
		int counter_edges=0;
		int coincide;
		for(int j=0;j<nedges;++j)
		{
			// Does the splitting edge cross the edge?
			SegmentIntersection(splitedge,edges[j],intersect);
			if(DistanceToEdge(intersect,edges[j])<1e-7*R)
			{
				splitedge.set_x(counter_edges,intersect.x);
				splitedge.set_y(counter_edges,intersect.y);
				++counter_edges;
				if(counter_edges==2)
					break;
			}
		}
		if(counter_edges!=2)
		{
			UniversalError eo("Error in splitting the cell");
			eo.AddEntry("Cell index",ToRefine[i]);
			throw eo;
		}
		// add the splitting edge
		old_ref.push_back((int)V.edges.size());
		new_ref.push_back((int)V.edges.size());
		V.edges.push_back(splitedge);
		// make the new edges
		for(int j=0;j<nedges;++j)
		{
			// Does the edge have a coinciding point with the splitting edge?
			if(splitedge.GetVertex(0).distance(edges[j].GetVertex(0))<
				1e-8*R||splitedge.GetVertex(0).distance(edges[j].
				GetVertex(1))<1e-8*R)
				coincide=0;
			else
				if(splitedge.GetVertex(1).distance(edges[j].GetVertex(0))<
					1e-8*R||splitedge.GetVertex(1).distance(edges[j].
					GetVertex(1))<1e-8*R)
					coincide=1;
				else
					coincide=2;
			if(coincide<2)
			{
				// No need to change vertices only neighbors
				if(DistanceToEdge(v,edges[j])>DistanceToEdge(NewPoint,edges[j]))
				{
					Edge NewEdge(edges[j]);
					if(NewEdge.GetNeighbor(0)==ToRefine[i])
						NewEdge.set_friend(0,Npoints+i);
					else
						NewEdge.set_friend(1,Npoints+i);
					V.edges[edge_index[j]]=NewEdge;
					new_ref.push_back(edge_index[j]);
				}
				else
					old_ref.push_back(edge_index[j]);
				continue;
			}
			// Does the splitting edge cross the edge?
			SegmentIntersection(splitedge,edges[j],intersect);
			if(DistanceToEdge(intersect,edges[j])<1e-8*R)
			{
				Edge NewEdge;
				NewEdge.set_friend(0,Npoints+i);
				if(edges[j].GetNeighbor(0)==ToRefine[i])
					NewEdge.set_friend(1,edges[j].GetNeighbor(1));
				else
					NewEdge.set_friend(1,edges[j].GetNeighbor(0));
				int index;
				if(NewPoint.distance(edges[j].GetVertex(0))<
					v.distance(edges[j].GetVertex(0)))
					index=0;
				else
					index=1;
				NewEdge.set_x(0,edges[j].get_x(index));
				NewEdge.set_y(0,edges[j].get_y(index));
				NewEdge.set_x(1,intersect.x);
				NewEdge.set_y(1,intersect.y);
				new_ref.push_back((int)V.edges.size());
				if((NewEdge.GetNeighbor(1)>-1)&&(NewEdge.GetNeighbor(1)<(n+Npoints)))
					V.mesh_vertices[NewEdge.GetNeighbor(1)].
					  push_back((int)V.edges.size());
				V.edges.push_back(NewEdge);
				// Do the other split
				NewEdge.set_x(0,edges[j].get_x((index+1)%2));
				NewEdge.set_y(0,edges[j].get_y((index+1)%2));
				NewEdge.set_friend(0,ToRefine[i]);
				V.edges[edge_index[j]]=NewEdge;
				old_ref.push_back(edge_index[j]);
				continue;
			}
			// No need to split the edge
			// check if update is needed
			if(NewPoint.distance(edges[j].GetVertex(0))<
				v.distance(edges[j].GetVertex(0)))
			{
				// Change edge neighbors
				if(edges[j].GetNeighbor(0)==ToRefine[i])
					V.edges[edge_index[j]].set_friend(0,Npoints+i);
				else
					V.edges[edge_index[j]].set_friend(1,Npoints+i);
				// add new reference
				new_ref.push_back(edge_index[j]);
			}
			else
				old_ref.push_back(edge_index[j]);
		}
		V.mesh_vertices[Npoints+i]=new_ref;
		V.mesh_vertices[ToRefine[i]]=old_ref;
	}
	// Calculate the new CM
	for(int i=0;i<(int)ToRefine.size();++i)
		V.CM.push_back(V.CalcCellCM(Npoints+i));
	for(int i=0;i<(int)ToRefine.size();++i)
		V.CM[ToRefine[i]]=V.CalcCellCM(ToRefine[i]);
	return;
	// Reset Tree if self gravity is needed
}

int VoronoiMesh::GetPointNo(void) const
{
  return Tri->get_length();
}

Vector2D VoronoiMesh::GetMeshPoint(int index) const
{
  return Tri->get_point(index);
}

int VoronoiMesh::GetTotalSidesNumber(void) const 
{
  return (int)edges.size();
}

Edge const& VoronoiMesh::GetEdge(int index) const 
{
  return edges[index];
}

Vector2D const& VoronoiMesh::GetCellCM(int index) const
{
  return CM[index];
}
