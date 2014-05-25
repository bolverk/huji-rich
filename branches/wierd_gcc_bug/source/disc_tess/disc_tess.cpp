#include "disc_tess.hpp"

DiscTess::DiscTess(double MaxR,double MinR,vector<double> R):_MaxR(MaxR),
	_MinR(MinR),_r(R){}

DiscTess::~DiscTess(void){}

void DiscTess::Initialise(vector<Vector2D> const& points,OuterBoundary const* /*bc*/)
{
	double r,dphi;
	int r_index;
	double eps=1e-8;
	pair<vector<double>::iterator,vector<double>::iterator> bounds;
	vector<vector<int> > point_indeces;
	vector<vector<double> > phi;
	point_indeces.resize(_r.size());
	phi.resize(_r.size());
	_mesh_vertices.resize(points.size());
	cor.resize(points.size());
	for(int i=0;i<(int)points.size();++i)
	{
		// Find the R index
		r=abs(points[i]);
		r_index=lower_bound(_r.begin(),_r.end(),r)-_r.begin();
		if(r_index==(int)_r.size())
			--r_index;
		else
			if(r_index>0)
				if(abs(r-_r[r_index])>abs(r-_r[r_index-1]))
					--r_index;
		if(abs(r-_r[r_index])>(0.1*r))
			throw("Bad initial setup, points are not on a circular ring");
		// Calculate the angle
		dphi=atan2(points[i].y,points[i].x);
		point_indeces[r_index].push_back(i);
		phi[r_index].push_back(dphi);
		cor[i]=points[i];
		// Fix radius if needed
		if(abs(r-_r[r_index])>(eps*r))
		{
			Vector2D temp(_r[r_index]*cos(dphi),_r[r_index]*cos(dphi));
			cor[i]=temp;
		}
	}

	// The edges phi
	vector<double> ephi,ephi_lower;
	vector<int> old_radial_e_index,radial_e_index;
	double phi_new;
	// Create the edges
	Edge edge;
	for(int i=0;i<(int)_r.size();++i)
	{
		//output("c:\\disco.bin");
		old_radial_e_index=radial_e_index;
		radial_e_index.clear();
		// copy the old radial edges angle
		ephi_lower=ephi;
		ephi.clear();
		int N=phi[i].size();
		// Find upper and lower annuli radii
		double Rlower,Rupper;
		if(i>0)
			Rlower=0.5*(_r[i]+_r[i-1]);
		else
			Rlower=0.5*(_r[i]+_MinR);
		if(i==(_r.size()-1))
			Rupper=0.5*(_MaxR+_r[i]);
		else
			Rupper=0.5*(_r[i+1]+_r[i]);

		// Make a local copy of the angles
		vector<double> phi_local=phi[i];
		vector<int> indeces;
		indeces.resize(phi_local.size());

		// sort the vector and get the sorted indeces
		sort_index(phi_local,indeces);
		sort(phi_local.begin(),phi_local.end());

		// Are we in the first radial bin?
		if(i==0)
		{
			for(int j=0;j<N;++j)
			{
				// Do radial edge
				phi_new=0.5*(phi_local[(j+1)%N]+phi_local[j]);
				if(phi_local[j]>0&&phi_local[(j+1)%N]<0)
					if(phi_new>0)
						phi_new-=PI;
					else
						phi_new+=PI;
				edge.set_x(0,Rlower*cos(phi_new));
				edge.set_y(0,Rlower*sin(phi_new));
				edge.set_x(1,Rupper*cos(phi_new));
				edge.set_y(1,Rupper*sin(phi_new));
				edge.set_friend(0,point_indeces[i][indeces[j]]);
				edge.set_friend(1,point_indeces[i][indeces[(j+1)%N]]);

				_mesh_vertices[edge.GetNeighbor(0)].push_back(_edges.size());
				_mesh_vertices[edge.GetNeighbor(1)].push_back(_edges.size());
				radial_e_index.push_back(_edges.size());
				ephi.push_back(phi_new);
				_edges.push_back(edge);
			}
			for(int j=0;j<N;++j)
			{
				// Do angular edge
				//phi_new=0.5*(phi_local[(j+1)%N]+phi_local[j]);
				phi_new=ephi[j];
				edge.set_x(0,Rlower*cos(phi_new));
				edge.set_y(0,Rlower*sin(phi_new));
				//phi_new=0.5*(phi_local[j]+phi_local[(j-1+N)%N]);
				phi_new=ephi[(j-1+N)%N];
				edge.set_x(1,Rlower*cos(phi_new));
				edge.set_y(1,Rlower*sin(phi_new));
				edge.set_friend(0,point_indeces[i][indeces[j]]);
				edge.set_friend(1,-1);

				_mesh_vertices[edge.GetNeighbor(0)].push_back(_edges.size());
				_edges.push_back(edge);
			}
			continue;
		}
		int Nedges=_edges.size();

		// Bulk annuli
		for(int j=0;j<N;++j)
		{
			// Do radial edge
			phi_new=0.5*(phi_local[(j+1)%N]+phi_local[j]);
			if(phi_local[j]>0&&phi_local[(j+1)%N]<0)
				if(phi_new>0)
					phi_new-=PI;
				else
					phi_new+=PI;
			edge.set_x(0,Rlower*cos(phi_new));
			edge.set_y(0,Rlower*sin(phi_new));
			edge.set_x(1,Rupper*cos(phi_new));
			edge.set_y(1,Rupper*sin(phi_new));
			edge.set_friend(0,point_indeces[i][indeces[j]]);
			edge.set_friend(1,point_indeces[i][indeces[(j+1)%N]]);

			_mesh_vertices[edge.GetNeighbor(0)].push_back(_edges.size());
			_mesh_vertices[edge.GetNeighbor(1)].push_back(_edges.size());
			radial_e_index.push_back(_edges.size());
			_edges.push_back(edge);
			ephi.push_back(phi_new);
		}
		//output("c:\\disco.bin");
		// Do angular edges
		// Sort the previous radial edges
		vector<int> indeces2;
		indeces2.resize(ephi_lower.size());

		// sort the vector and get the sorted indeces
		sort_index(ephi_lower,indeces2);
		sort(ephi_lower.begin(),ephi_lower.end());
		ReArrangeVector(old_radial_e_index,indeces2);
		for(int j=0;j<N;++j)
		{	
			//output("c:\\disco.bin");

			// Find the relevant radial edges
			// The right most edge
			int index_r=lower_bound(ephi_lower.begin(),ephi_lower.end(),ephi[j])
				-ephi_lower.begin();
			if(index_r==((int)ephi_lower.size()))
				index_r=0;
			else
				if(ephi_lower[index_r]==ephi[j])
					index_r=(index_r+1)%ephi_lower.size();

			// Are there no lower edges cutting this angular edge?
			if((ephi_lower[index_r]>=ephi[(j+1)%N]))
				//((ephi_lower[index_r]<=ephi[j])&&(ephi_lower[index_r]*ephi[j]>0)))
			{
				if(!((ephi[j]>0&&ephi[(j+1)%N]<0)&&ephi_lower[index_r]>0))
				{
					edge.set_x(0,Rlower*cos(ephi[j]));
					edge.set_y(0,Rlower*sin(ephi[j]));
					edge.set_x(1,Rlower*cos(ephi[(j+1)%N]));
					edge.set_y(1,Rlower*sin(ephi[(j+1)%N]));
					edge.set_friend(0,point_indeces[i][indeces[(j+1)%N]]);
					edge.set_friend(1,_edges[old_radial_e_index[index_r]].GetNeighbor(0));

					_mesh_vertices[edge.GetNeighbor(0)].push_back(_edges.size());
					_mesh_vertices[edge.GetNeighbor(1)].push_back(_edges.size());
					_edges.push_back(edge);
					continue;
				}
			}

			// The left most edge
			int index_l=lower_bound(ephi_lower.begin(),ephi_lower.end(),
				ephi[(j+1+N)%N])-ephi_lower.begin();
			if(index_l==((int)ephi_lower.size()))
				index_l=0;
			if((ephi_lower[index_l]>ephi[(j+1)%N])||
				((ephi_lower[index_l]<0)&&(ephi[(j+1)%N]>0)))
				index_l=(index_l-1+ephi_lower.size())%ephi_lower.size();

			// create the angular edges
			edge.set_x(0,Rlower*cos(ephi[j]));
			edge.set_y(0,Rlower*sin(ephi[j]));
			edge.set_x(1,_edges[old_radial_e_index[index_r]].get_x(1));
			edge.set_y(1,_edges[old_radial_e_index[index_r]].get_y(1));
			edge.set_friend(0,point_indeces[i][indeces[(j+1)%N]]);
			edge.set_friend(1,_edges[old_radial_e_index[index_r]].GetNeighbor(0));

			_mesh_vertices[edge.GetNeighbor(0)].push_back(_edges.size());
			_mesh_vertices[edge.GetNeighbor(1)].push_back(_edges.size());
			_edges.push_back(edge);

			edge.set_x(0,_edges[old_radial_e_index[index_l]].get_x(1));
			edge.set_y(0,_edges[old_radial_e_index[index_l]].get_y(1));
			edge.set_x(1,Rlower*cos(ephi[(j+1+N)%N]));
			edge.set_y(1,Rlower*sin(ephi[(j+1+N)%N]));
			edge.set_friend(1,_edges[old_radial_e_index[index_l]].GetNeighbor(1));

			_mesh_vertices[edge.GetNeighbor(0)].push_back(_edges.size());
			_mesh_vertices[edge.GetNeighbor(1)].push_back(_edges.size());
			_edges.push_back(edge);
			// Are there edges in between?
			if(index_l!=index_r)
			{
				//output("c:\\disco.bin");
				for(int jj=0;jj<((index_l-index_r+N)%N);++jj)
				{
					edge.set_x(1,edge.get_x(0));
					edge.set_y(1,edge.get_y(0));
					edge.set_x(0,_edges[old_radial_e_index[(index_l-1-jj+N)%N]].
						get_x(1));
					edge.set_y(0,_edges[old_radial_e_index[(index_l-1-jj+N)%N]].
						get_y(1));
					edge.set_friend(1,_edges[old_radial_e_index[(index_l-1-jj+N)%N]].
						GetNeighbor(1));
				}
			}
		}	

		// Are we in the last annuli?
		if(i==(_r.size()-1))
		{
			//output("c:\\disco.bin");
			// Do the additional angular edges
			for(int j=0;j<N;++j)
			{
				phi_new=0.5*(phi_local[(j+1)%N]+phi_local[j]);
				if(phi_local[j]>0&&phi_local[(j+1)%N]<0)
					if(phi_new>0)
						phi_new-=PI;
					else
						phi_new+=PI;
				edge.set_x(0,Rupper*cos(phi_new));
				edge.set_y(0,Rupper*sin(phi_new));
				phi_new=0.5*(phi_local[j]+phi_local[(j-1+N)%N]);
				if(phi_local[(j-1+N)%N]>0&&phi_local[j]<0)
					if(phi_new>0)
						phi_new-=PI;
					else
						phi_new+=PI;
				edge.set_x(1,Rupper*cos(phi_new));
				edge.set_y(1,Rupper*sin(phi_new));
				edge.set_friend(1,point_indeces[i][indeces[j]]);
				edge.set_friend(0,-1);

				_mesh_vertices[edge.GetNeighbor(1)].push_back(_edges.size());
				_edges.push_back(edge);
			}
		}
		// sort the indeces for next step
		indeces.resize(ephi.size());
		sort_index(ephi,indeces);
		ReArrangeVector(radial_e_index,indeces);
		ReArrangeVector(ephi,indeces);
	}
}

void DiscTess::Update(vector<Vector2D> const& points)
{
	Initialise(points);
}

int DiscTess::GetPointNo(void) const
{
	return cor.size();
}

Vector2D DiscTess::GetMeshPoint(int index) const
{
	return cor[index];
}

Vector2D DiscTess::GetCellCM(int index) const
{
	return cor[index];
}

int DiscTess::GetTotalSidesNumber(void) const
{
	return (int)_edges.size();
}

Edge const& DiscTess::GetEdge(int index) const
{
	return _edges[index];
}

double DiscTess::GetWidth(int index) const
{
	return sqrt(GetVolume(index)/PI);
}

double DiscTess::GetVolume(int index) const
{
	double phi;
	double phi_mesh=atan2(cor[index].y,cor[index].x);
	for(int i=0;i<(int)_mesh_vertices[index].size();++i)
	{
		phi=atan2(_edges[_mesh_vertices[index][i]].get_y(1),
			_edges[_mesh_vertices[index][i]].get_x(1));
		// Is it a radial edge?
		if(abs(phi-atan2(_edges[_mesh_vertices[index][i]].get_y(0),
			_edges[_mesh_vertices[index][i]].get_x(0)))<0.1)
		{
			// Same sign as phi_mesh?
			if(phi*phi_mesh>=0)
			{
				double rmin=abs(_edges[_mesh_vertices[index][i]].GetVertex(0));
				double rmax=abs(_edges[_mesh_vertices[index][i]].GetVertex(1));
				return 0.5*abs(phi-phi_mesh)*(rmax*rmax-rmin*rmin);
			}
		}
	}
	throw("Error in GetVolume, wrong signs");
}

vector<int>const& DiscTess::GetCellEdges(int index) const
{
	return _mesh_vertices[index];
}

int DiscTess::GetOriginalIndex(int point) const
{
	return point;
}

vector<Vector2D>& DiscTess::GetMeshPoints(void)
{
	return cor;
}

vector<int> DiscTess::GetNeighbors(int index)const
{
	vector<int> res;
	int other;
	res.reserve(6);
	for(int i=0;i<(int)_mesh_vertices[index].size();++i)
	{
		if((other=_edges[_mesh_vertices[index][i]].GetNeighbor(0))!=index)
			res.push_back(other);
		else
			res.push_back(_edges[_mesh_vertices[index][i]].GetNeighbor(1));
	}
	return res;
}

Tessellation* DiscTess::clone(void) const
{
	return new DiscTess(*this);
}

bool DiscTess::NearBoundary(int index) const
{
	int n=_mesh_vertices[index].size();
	int n1,n0;
	for(int i=0;i<n;++i)
	{
		n0=_edges[_mesh_vertices[index][i]].GetNeighbor(0);
		n1=_edges[_mesh_vertices[index][i]].GetNeighbor(1);
		if(n0<0||n1<0||n0>n||n1>n)
			return true;
	}
	return false;
}

void DiscTess::RemoveCells(vector<int> &ToRemovevector,vector<vector<int> > &VolIndex,
	vector<vector<double> > &Volratio)
{
	// Not implemented yet
	return;
}

void DiscTess::RefineCells(vector<int> const& ToRefine,vector<Vector2D> const&
	directions,double alpha)
{
	// Not implemented yet
	return;
}

void DiscTess::output(string filename)
{	
	int temp2=_edges.size();
	fstream myFile(filename.c_str(),ios::out | ios::binary);
	myFile.write ((char*)&temp2, sizeof (int));
	double temp;
	for(int i=0;i<temp2;i++)
	{
		temp=_edges[i].get_x(0);
		myFile.write ((char*)&temp, sizeof (double));
		temp=_edges[i].get_x(1);
		myFile.write ((char*)&temp, sizeof (double));
	}
	for(int i=0;i<temp2;i++)
	{
		temp=_edges[i].get_y(0);
		myFile.write ((char*)&temp, sizeof (double));
		temp=_edges[i].get_y(1);
		myFile.write ((char*)&temp, sizeof (double));
	}
	int j;
	for(int i=0;i<temp2;i++)
	{
		j=_edges[i].GetNeighbor(0);
		myFile.write ((char*)&j, sizeof (int));
		j=_edges[i].GetNeighbor(1);
		myFile.write ((char*)&j, sizeof (int));
	}
	temp2=cor.size();
	myFile.write ((char*)&temp2, sizeof (int));
	for(int i=0;i<temp2;i++)
	{
		temp=cor[i].x;
		myFile.write ((char*)&temp, sizeof (double));
		temp=cor[i].y;
		myFile.write ((char*)&temp, sizeof (double));
	}
	vector<int> indexes;
	int length;
	for(int i=0;i<temp2;i++)
	{
		indexes=GetCellEdges(i);
		length=indexes.size();
		myFile.write((char*)&length,sizeof(int));
		for(int j=0;j<(int)indexes.size();j++)
			myFile.write((char*)&indexes[j],sizeof(int));
	}
	myFile.close();
}

Vector2D DiscTess::CalcFaceVelocity(Vector2D wl, Vector2D wr,Vector2D rL,
	Vector2D rR,Vector2D f)const
{
	double r=abs(rL);
	Vector2D phi_hat_l(-rL.y/r,rL.x/r);
	double vphi_l=ScalarProd(phi_hat_l,wl);
	r=abs(rR);
	Vector2D phi_hat_r(-rR.y/r,rR.x/r);
	double vphi_r=ScalarProd(phi_hat_r,wr);
	return 0.5*(vphi_l+vphi_r)*Vector2D(-f.y,f.x)/abs(f);
}

vector<Vector2D> DiscTess::calc_edge_velocities(HydroBoundaryConditions const* hbc,
	vector<Vector2D> const& point_velocities,double time)const
{
	vector<Vector2D> facevelocity;
	facevelocity.resize(_edges.size());
	for(int i = 0; i < (int)_edges.size(); ++i)
	{
		if(hbc->IsBoundary(_edges[i],this))
		{
			// Boundary
			facevelocity[i] = hbc->CalcEdgeVelocity(this,point_velocities,_edges[i],
				time);
		}
		else
		{
			// Bulk
			facevelocity[i] = CalcFaceVelocity(
				point_velocities[_edges[i].GetNeighbor(0)],
				point_velocities[_edges[i].GetNeighbor(1)],
				GetMeshPoint(_edges[i].GetNeighbor(0)),
				GetMeshPoint(_edges[i].GetNeighbor(1)),
				0.5*(_edges[i].GetVertex(0)+_edges[i].GetVertex(1)));
		}
	}
	return facevelocity;

}