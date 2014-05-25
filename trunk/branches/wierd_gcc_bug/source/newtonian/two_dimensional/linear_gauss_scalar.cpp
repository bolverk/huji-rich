#include "linear_gauss_scalar.hpp"

namespace{
Vector2D calc_centroid(Edge const& edge)
{
	return 0.5*(edge.GetVertex(0)+
		edge.GetVertex(1));
}

double interp(double cell,Vector2D const& cell_cm,Vector2D const& slope,
	Vector2D const& target)
{
	return cell+ScalarProd((target-cell_cm),slope);
}
}

namespace
{

	vector<Edge> get_edge_list(Tessellation const* tess,
		vector<int> const& edge_indices)
	{
		vector<Edge> res(edge_indices.size());
		for(int i=0;i<(int)edge_indices.size();++i){
			res[i] = tess->GetEdge(edge_indices[i]);
		}
		return res;
	}

	vector<vector<double> > GetNeighborTracers(Tessellation const* tess,
		vector<Edge> const&	edges,int cell_index,
		vector<vector<double> > const& tracers,
		HydroBoundaryConditions const* hbc,double time)
	{
	  int n=(int)edges.size();
		//		int n2=tracers[0].size();
		int neigh0,neigh1;
		vector<vector<double> > res(n);
		for(int i=0;i<n;++i)
		{
			if(hbc->IsBoundary(edges[i],tess))
			{
				res[i]=hbc->GetBoundaryTracers(edges[i],tess,tracers,time);
			}
			else
			{
				neigh0=edges[i].GetNeighbor(0);
				neigh1=edges[i].GetNeighbor(1);
				if(neigh0==cell_index)
					res[i]=tracers[neigh1];
				else
					res[i]=tracers[neigh0];
			}
		}
		return res;
	}

	Vector2D GetReflectedPoint(Tessellation const* tess,int point,
		Edge const& edge)
	{
		Vector2D MeshPoint=tess->GetMeshPoint(point);
		Vector2D par=edge.GetVertex(1)-edge.GetVertex(0);
		if(abs(par.x)>abs(par.y))
		{
			// We are in the x direction
			if(MeshPoint.y>edge.get_y(0))
				MeshPoint.y-=MeshPoint.y-edge.get_y(0);
			else
				MeshPoint.y+=edge.get_y(0)-MeshPoint.y;
		}
		else
		{
			// We are in the y direction
			if(MeshPoint.x>edge.get_x(0))
				MeshPoint.x-=MeshPoint.x-edge.get_x(0);
			else
				MeshPoint.x+=edge.get_x(0)-MeshPoint.x;
		}
		return MeshPoint;
	}


	vector<Vector2D> GetNeighborMesh(Tessellation const* tess,vector<Edge> const&
		edges,int cell_index)
	{
	  int n=(int)edges.size();
		int neigh0,neigh1;
		vector<Vector2D> res(n);
		for(int i=0;i<n;++i)
		{
			neigh0=edges[i].GetNeighbor(0);
			neigh1=edges[i].GetNeighbor(1);
			if(neigh0==cell_index)
				if(neigh1>-1) // we are not near rigid wall
					res[i]=tess->GetMeshPoint(neigh1);
				else
					res[i]=GetReflectedPoint(tess,neigh0,edges[i]);
			else
				if(neigh0>-1)
					res[i]=tess->GetMeshPoint(neigh0);
				else
					res[i]=GetReflectedPoint(tess,neigh1,edges[i]);
		}
		return res;
	}


	vector<Vector2D> calc_naive_slope(vector<double> cell,Vector2D const& center,
		double cell_volume,vector<vector<double> > const& neighbors,
		vector<Vector2D> const& neighbor_centers,vector<Edge> const& edge_list)
	{
		vector<Vector2D> res(cell.size());
		cell_volume=1/cell_volume;
		for(int i=0;i<(int)edge_list.size();++i)
		{
			Vector2D c_ij = calc_centroid(edge_list[i])-
				0.5*(neighbor_centers[i]+center);
			Vector2D r_ij = center - neighbor_centers[i];
			double rij_1 = 1/abs(r_ij);
			for(int j=0;j<(int)cell.size();++j)
				res[j] += (edge_list[i].GetLength()*cell_volume)*
				((neighbors[i][j]-cell[j])*(c_ij*rij_1)-
				0.5*(cell[j]+neighbors[i][j])*(r_ij*rij_1));
		}
		return res;
	}

	vector<Vector2D> slope_limit(vector<double> cell,Vector2D const& center,
		vector<vector<double> > const& neighbors,
		vector<Edge> const& edge_list,vector<Vector2D> const& slope)
	{
		vector<Vector2D> res = slope;
		vector<vector<double> > centroid_vars(edge_list.size());
		for(int i=0;i<(int)edge_list.size();++i)
		{
			centroid_vars[i].resize(cell.size());
			for(unsigned int j=0;j<cell.size();++j)
				centroid_vars[i][j] = interp(cell[j],center,slope[j],
				calc_centroid(edge_list[i]));
		}
		for(unsigned int j=0;j<cell.size();++j)
		{
			double psi = 1;
			double min_val = min(neighbors[0][j],cell[j]);
			double max_val = max(neighbors[0][j],cell[j]);
			for(unsigned int i=1;i<edge_list.size();++i)
			{
				min_val = min(neighbors[i][j],min_val);
				max_val = max(neighbors[i][j],max_val);
			}
			for(int i=0;i<(int)edge_list.size();++i)
			{
				if(centroid_vars[i][j]>cell[j])
					psi = min(psi,(max_val-cell[j])/
					(centroid_vars[i][j]-cell[j]));
				else
					if(centroid_vars[i][j]<cell[j])
						psi = min(psi,(min_val-cell[j])/
						(centroid_vars[i][j]-cell[j]));
					else
						if(centroid_vars[i][j]==cell[j])
							psi = min(psi,1.);
						else
						{
							UniversalError eo("A strange error occurred in slope limit. Probably a nan in cell");
							throw eo;
						}
			}
			res[j].x *= psi;
			res[j].y *= psi;
		}
		return res;
	}

	vector<Vector2D> slope_limit2(vector<double> cell,Vector2D const& center,
		vector<vector<double> > const& neighbors,
		vector<Edge> const& edge_list,vector<Vector2D> const& slope,double theta)
	{
		vector<Vector2D> res = slope;
		vector<vector<double> > centroid_vars(edge_list.size());
		for(int i=0;i<(int)edge_list.size();++i)
		{
			centroid_vars[i].resize(cell.size());
			for(unsigned int j=0;j<cell.size();++j)
				centroid_vars[i][j] = interp(cell[j],center,slope[j],
				calc_centroid(edge_list[i]));
		}
		for(unsigned int j=0;j<cell.size();++j)
		{
			double psi = 1;
			double dphi;
			double maxdiff=0;
			for(int i=0;i<(int)edge_list.size();++i)
				maxdiff=max(maxdiff,abs(neighbors[i][j]-cell[j]));
			for(int i=0;i<(int)edge_list.size();++i)
			{
				dphi=centroid_vars[i][j]-cell[j];
				if(abs(dphi)<0.1*maxdiff&&(centroid_vars[i][j]*cell[j]>0))
					continue;
				if(dphi>0)
					psi = min(psi,theta*max((neighbors[i][j]-cell[j])/dphi,0.0));
				else 
					if(dphi<0)
						psi=min(psi,theta*max((neighbors[i][j]-cell[j])/dphi,0.0));
					else 
						if(dphi==0)
							psi = min(psi,1.);
						else
							throw UniversalError("A strange error occurred in slope limit. Probably a nan in cell");
			}
			res[j].x *= psi;
			res[j].y *= psi;
		}
		return res;
	}
}

namespace {
vector<Vector2D> calc_slope(Tessellation const* tess,
	vector<vector<double> > const& tracers,
	int cell_index, HydroBoundaryConditions const* bc,
	bool limit,double time,double theta)
{
	vector<int> edge_indices = tess->GetCellEdges(cell_index);
	vector<Edge> edge_list = get_edge_list(tess,edge_indices);
	vector<Vector2D> neighbor_mesh_list =GetNeighborMesh(tess,edge_list,cell_index);
	vector<vector<double> > neighbors_tracers = GetNeighborTracers(tess,edge_list,cell_index,
		tracers,bc,time);
	const vector<Vector2D> naive_slope = calc_naive_slope(tracers[cell_index],
		tess->GetMeshPoint(cell_index),tess->GetVolume(cell_index),
		neighbors_tracers,neighbor_mesh_list,edge_list);
	if(limit)
		return slope_limit2(tracers[cell_index],tess->GetCellCM(cell_index),
		neighbors_tracers,edge_list,naive_slope,theta);
	else
		return slope_limit(tracers[cell_index],tess->GetCellCM(cell_index),
		neighbors_tracers,edge_list,naive_slope);
}

  void calc_slopes(Tessellation const* /*tess*/,
		   vector<vector<double> > const& /*tracers*/,
		   HydroBoundaryConditions const* /*bc*/,
		   SpatialReconstruction const* /*sr*/,
		   vector<vector<Vector2D> >& /*slopes*/,
		   double /*time*/,
		   double /*theta*/,
		   HydroBoundaryConditions const* /*hbc*/)
{
  throw;
  /*
	for(int i=0;i<(int)tracers.size();++i)
		if(!hbc->IsGhostCell(i,tess))
			slopes[i] = calc_slope(tess,tracers,i,bc,sr->WasSlopeLimited(i),time,
			theta);		
  */
}
}

vector<double> LinearGaussScalar::Interpolate
(Tessellation const* tessellation,
 vector<vector<double> > const& tracers,double /*dt*/,
 Edge const& edge,int side,SInterpolationType interptype)const
{
	int cell_index = edge.GetNeighbor(side);
	Vector2D target = calc_centroid(edge);
	vector<double> res(tracers[0].size());
	if(interptype==SInBulk)
	{
		for(unsigned int j=0;j<tracers[0].size();++j)
			res[j]=interp(tracers[cell_index][j],
			tessellation->GetCellCM(cell_index),
			slopes_[cell_index][j], target);
		return res;
	}
	else
		if(interptype==SBoundary)
		{
			int other=edge.GetNeighbor((side+1)%2);
			for(unsigned int j=0;j<tracers[0].size();++j)
				res[j]=interp(tracers[tessellation->GetOriginalIndex(other)][j],
				tessellation->GetCellCM(other),
				slopes_[tessellation->GetOriginalIndex(other)][j], target);
			return res;
		}
		else
			throw UniversalError("Wrong interpolation type in linear_gauss_scalar");
}

LinearGaussScalar::LinearGaussScalar
(SpatialReconstruction const* sr,
 HydroBoundaryConditions const* hbc,double theta):
  sr_(sr),
  slopes_(vector<vector<Vector2D> >()),
  hbc_(hbc),theta_(theta){}

LinearGaussScalar::~LinearGaussScalar(){}

void LinearGaussScalar::Prepare(Tessellation const* tessellation,
				vector<vector<double> > const& tracers,
				double /*dt*/,double time)
{
	if(slopes_.size()!=tracers.size())
		slopes_.resize(tracers.size());
	calc_slopes(tessellation,tracers,hbc_,sr_,slopes_,time,theta_,hbc_);
}
