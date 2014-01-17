#include <limits>
#include "linear_gauss.hpp"

namespace {

	PrimitiveGradient2D operator*(Vector2D const& v,
		Primitive const& p)
	{
		return PrimitiveGradient2D(v.x*p,
			v.y*p);
	}

	PrimitiveGradient2D operator*(Primitive const& p,
		Vector2D const& v)
	{
		return PrimitiveGradient2D(v.x*p,
			v.y*p);
	}

	PrimitiveGradient2D operator*(double s,
		PrimitiveGradient2D const& pg)
	{
		return PrimitiveGradient2D(s*pg.x,s*pg.y);
	}

	PrimitiveGradient2D operator-(PrimitiveGradient2D const& pg1,
		PrimitiveGradient2D const& pg2)
	{
		return PrimitiveGradient2D(pg1.x-pg2.x,
			pg1.y-pg2.y);
	}

	PrimitiveGradient2D operator+(PrimitiveGradient2D const& pg1,
		PrimitiveGradient2D const& pg2)
	{
		return PrimitiveGradient2D(pg1.x+pg2.x,
			pg1.y+pg2.y);
	}

	Primitive operator*(Vector2D const& v,
		PrimitiveGradient2D const& pg)
	{
		return v.x*pg.x+v.y*pg.y;
	}
}

namespace {
	vector<int> GetOnlyRealNeighbors(Tessellation const* tess,int cell_index)
	{
		vector<int> Neigh=tess->GetNeighbors(cell_index);
		vector<int> res;
		res.reserve(Neigh.size());
		for(int i=0;i<(int)Neigh.size();++i)
			if(Neigh[i]>-1)
				res.push_back(Neigh[i]);
		return res;
	}

	vector<Primitive> GetRealNeighborPrimitive
		(Tessellation const* tess,
		vector<Primitive> const& cells,vector<int> const& Neigh)
	{
		vector<Primitive> res;
		for(int i=0;i<(int)Neigh.size();++i){
			if(Neigh[i]>-1)
				res.push_back(cells[tess->GetOriginalIndex(Neigh[i])]);
		}
		return res;
	}

	vector<Vector2D> GetRealNeighborMesh
		(Tessellation const* tess,
		vector<int> const& neigh)
	{
		vector<Vector2D> res;
		res.resize(neigh.size());
		for(int i=0;i<(int)neigh.size();++i){
			res[i]=tess->GetCellCM(neigh[i]);
		}
		return res;
	}

	Vector2D LinearGrad
		(vector<Vector2D> const& X,vector<double> const& val,
		Vector2D const& X0,double val0)
	{
		vector<double> dx,dy;
		dx.resize(X.size());
		dy.resize(X.size());
		for(int i=0;i<(int)X.size();++i){
			dx[i]=X[i].x-X0.x;
			dy[i]=X[i].y-X0.y;
		}
		double c1=0,c2=0,c3=0,c4=0,c5=0,r_inv;
		for(int i=0;i<(int)X.size();++i)
		{
			r_inv=1/(dx[i]*dx[i]+dy[i]*dy[i]);
			c1+=dy[i]*dy[i]*r_inv;
			c2+=dx[i]*dx[i]*r_inv;
			c3+=dx[i]*dy[i]*r_inv;
			c4+=dx[i]*(val[i]-val0)*r_inv;
			c5+=dy[i]*(val[i]-val0)*r_inv;
		}
		double b=(c5-c4*c3/c2)/(c1-c3*c3/c2);
		return Vector2D((c4-b*c3)/c2,b);
	}

	PrimitiveGradient2D LinearSlope
		(vector<Primitive> const& cells,
		Tessellation const* tess,int point)
	{
		vector<int> neigh=GetOnlyRealNeighbors(tess,point);
		vector<Primitive> neighbors=GetRealNeighborPrimitive(tess,cells,neigh);
		vector<Vector2D> neighbor_mesh=GetRealNeighborMesh(tess,neigh);
		Vector2D mesh_point=tess->GetCellCM(point);
		PrimitiveGradient2D res;
		vector<double> vals(neigh.size());
		Vector2D grad;
		for(int i=0;i<cells[0].GetVarNo();++i)
		{
			for(int j=0;j<(int)neigh.size();++j)
				vals[j]=neighbors[j][i];
			grad=LinearGrad(neighbor_mesh,vals,mesh_point,cells[point][i]);
			res.x[i]=grad.x;
			res.y[i]=grad.y;
		}
		return res;
	}
}

template<typename T> bool is_it_finite(T arg)
{
	return arg == arg && 
		arg != std::numeric_limits<T>::infinity() &&
		arg != -std::numeric_limits<T>::infinity();
}

namespace {
	vector<Edge> RemoveBoundaryEdge(vector<Edge> const& edges)
	{
		vector<Edge> res;
		res.reserve(edges.size());
		for(int i=0;i<(int)edges.size();++i)
		{
			if(edges[i].GetNeighbor(0)!=-1&&edges[i].GetNeighbor(1)!=-1)
				res.push_back(edges[i]);
		}
		return res;
	}

	bool BoundaryNeighbor(vector<int> const& neighbors)
	{
		for(int i=0;i<(int)neighbors.size();++i)
			if(neighbors[i]==-1)
				return true;
		return false;
	}

	bool IsBoundaryPoint(Tessellation const* tess,int point)
	{
		return BoundaryNeighbor(tess->GetNeighbors(point));
	}

	UniversalError non_finite_value
		(vector<Vector2D> const& points,
		vector<Primitive> const& neighbors,vector<Edge> const& edges)
	{
		UniversalError res("Non finite value in Linear Guass naive slope");
		for(int i=0;i<(int) points.size();++i)
		{
			res.AddEntry("Neighbor X cor",points[i].x);
			res.AddEntry("Neighbor Y cor",points[i].y);
			res.AddEntry("Neighbor pressure",neighbors[i].Pressure);
			res.AddEntry("Neighbor density",neighbors[i].Density);
			res.AddEntry("Neighbor velocity x",neighbors[i].Velocity.x);
			res.AddEntry("Neighbor velocity y",neighbors[i].Velocity.y);
		}
		for(int i=0;i<(int) edges.size();++i)
		{
			res.AddEntry("Edge X cor",edges[i].get_x(0));
			res.AddEntry("Edge Y cor",edges[i].get_y(0));
			res.AddEntry("Edge X cor",edges[i].get_x(1));
			res.AddEntry("Edge Y cor",edges[i].get_y(1));
			res.AddEntry("Edge length",edges[i].GetLength());
		}
		return res;
	}

	UniversalError non_finite_slope_limit
		(Primitive const& cell,
		vector<Primitive> const& neighbors)
	{
		UniversalError res("Non finite value in Linear Guass slope limit");
		res.AddEntry("Cell pressure",cell.Pressure);
		res.AddEntry("Cell density",cell.Density);
		res.AddEntry("Cell velocity x",cell.Velocity.x);
		res.AddEntry("Cell velocity y",cell.Velocity.y);
		for(int i=0;i<(int) neighbors.size();++i)
		{
			res.AddEntry("Interpolated pressure",neighbors[i].Pressure);
			res.AddEntry("Interpolated density",neighbors[i].Density);
			res.AddEntry("Interpolated velocity x",neighbors[i].Velocity.x);
			res.AddEntry("Interpolated velocity y",neighbors[i].Velocity.y);
		}
		return res;
	}
}

LinearGauss::LinearGauss
	(OuterBoundary const& obc,
	HydroBoundaryConditions const* hbc,
	bool slf, bool soitf,double delta_v,double theta,
	double delta_P,bool rigidflag):
slope_limited_(vector<bool>()),
	slopes_(vector<PrimitiveGradient2D>()),
	obc_(&obc),
	hbc_(hbc),slf_(slf), soitf_(soitf),shockratio_(delta_v),
	diffusecoeff_(theta),
	pressure_ratio_(delta_P),_rigidflag(rigidflag) {}

namespace {
	int GetNeighborIndex(int cell_index,
		Edge const& edge)
	{
		if(edge.GetNeighbor(0)==cell_index)
			return edge.GetNeighbor(1);
		else if(edge.GetNeighbor(1)==cell_index)
			return edge.GetNeighbor(0);
		else
			throw UniversalError("Error in GetNeighborIndex: edge does not bound cell");
	}

	Primitive interp(Primitive const& cell,
		Vector2D const& cell_cm,
		PrimitiveGradient2D const& slope,
		Vector2D const& target)
	{
		return cell+(target-cell_cm)*slope;
	}

	Vector2D GetReflectedPoint
		(Tessellation const* tess,int point,
		OuterBoundary const& /*obc*/,Edge const& edge)
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

	/*
	Primitive GetReflectedCell
	(Tessellation const* tess,int point,
	Edge const& edge,vector<Primitive> const& cells)
	{
	Primitive res=cells[point];
	Vector2D par=edge.GetVertex(1)-edge.GetVertex(0);
	if(abs(par.x)>abs(par.y))
	{
	// We are in the x direction
	res.Velocity.y*=-1;
	}
	else
	{
	// We are in the y direction
	res.Velocity.x*=-1;
	}
	return res;
	}
	*/

	vector<Vector2D> GetNeighborMesh
		(Tessellation const* tess,vector<Edge> const&
		edges,int cell_index,OuterBoundary const& obc)
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
					res[i]=GetReflectedPoint(tess,neigh0,obc,edges[i]);
			else
				if(neigh0>-1)
					res[i]=tess->GetMeshPoint(neigh0);
				else
					res[i]=GetReflectedPoint(tess,neigh1,obc,edges[i]);
		}
		return res;
	}

	vector<Primitive> GetNeighborPrimitive
		(Tessellation const* tess,vector<Edge> const&
		edges,int cell_index,OuterBoundary const& /*obc*/,
		vector<Primitive> const& cells,HydroBoundaryConditions const* hbc,
		double time)
	{
		int n=(int)edges.size();
		int neigh0,neigh1;
		vector<Primitive> res(n);
		for(int i=0;i<n;++i)
		{
			if(hbc->IsBoundary(edges[i],tess))
			{
				res[i]=hbc->GetBoundaryPrimitive(edges[i],tess,cells,time);
			}
			else
			{
				neigh0=edges[i].GetNeighbor(0);
				neigh1=edges[i].GetNeighbor(1);
				if(neigh0==cell_index)
					res[i]=cells[neigh1];
				else
					res[i]=cells[neigh0];
			}
		}
		return res;
	}

	PrimitiveGradient2D CalcNaiveSlope
		(Primitive const& cell,
		Vector2D const& center,
		double cell_volume,
		vector<Primitive> const& neighbors,
		vector<Vector2D> const& neighbor_centers,
		vector<Edge> const& edge_list)
	{
		PrimitiveGradient2D res;
		Primitive phi_i = cell;
		const double cell_volume_inv=1.0/cell_volume;
		int iter_num=(int)edge_list.size();
		for(int i=0;i<iter_num;++i)
		{
			Vector2D c_ij = CalcCentroid(edge_list[i])-
				0.5*(neighbor_centers[i]+center);
			Vector2D r_ij = center - neighbor_centers[i];
			double rij_1=1/abs(r_ij);
			Primitive phi_j = neighbors[i];
			res += (edge_list[i].GetLength()*cell_volume_inv)*
				((phi_j-phi_i)*(c_ij*rij_1)-
				0.5*(phi_i+phi_j)*(r_ij*rij_1));
		}
		for(int i=0;i<phi_i.GetVarNo();++i)
		{
			if(!is_it_finite(res.x[i])||!is_it_finite(res.y[i]))
			{
				//	    cout<<"Nan in calc naive slope"<<endl;
				throw non_finite_value(neighbor_centers,neighbors,edge_list);
			}
		}
		return res;
	}

	bool is_between(double xl, double xm, double xr)
	{
		return (xr-xm)*(xm-xl)<=0;		 
	}

	bool is_between(Primitive const& pl,
		Primitive const& pm,
		Primitive const& pr)
	{
		for(int i=0;i<pl.GetVarNo();++i){
			if(!is_between(pl[i],pm[i],pr[i]))
				return false;
		}
		return true;
	}

	double get_max_prop(vector<Primitive> const& vp,
		int prop_idx)
	{
		double res = vp[0][prop_idx];
		for(int i=1;i<(int)vp.size();++i){
			res = max(res,vp[i][prop_idx]);
		}
		return res;
	}

	double get_min_prop(vector<Primitive> const& vp,
		int prop_idx)
	{
		double res = vp[0][prop_idx];
		for(int i=1;i<(int)vp.size();++i){
			res = min(res,vp[i][prop_idx]);
		}
		return res;
	}

	PrimitiveGradient2D slope_limit(Primitive const& cell,
		Vector2D const& center,
		vector<Primitive> const& neighbors,
		vector<Edge> &edge_list,
		PrimitiveGradient2D const& slope)
	{
		PrimitiveGradient2D res = slope;
		vector<Primitive> centroid_vars(edge_list.size());
		for(int i=0;i<(int)edge_list.size();++i)
			centroid_vars[i] = interp(cell,
			center,
			slope,
			CalcCentroid(edge_list[i]));
		for(int i=0;i<cell.GetVarNo();++i){
			double psi = 1;
			double min_prop = min(get_min_prop(neighbors,i),cell[i]);
			double max_prop = max(get_max_prop(neighbors,i),cell[i]);
			double maxdiff=max(max_prop-cell[i],cell[i]-min_prop);
			for(int j=0;j<(int)edge_list.size();++j)
			{
				double dphi=centroid_vars[j][i]-cell[i];
				if(abs(dphi)<0.1*maxdiff&&(centroid_vars[j][i]*cell[i]>0))
					continue;
				if(centroid_vars[j][i]>cell[i])
					psi = min(psi,(max_prop-cell[i])/(centroid_vars[j][i]-cell[i]));
				else if(centroid_vars[j][i]<cell[i])
					psi = min(psi,(min_prop-cell[i])/(centroid_vars[j][i]-cell[i]));
				else if(centroid_vars[j][i]==cell[i])
					psi = min(psi,1.);
				else
				{
					/*	      cout<<edge_list[j].GetNeighbor(0)<<" "<<edge_list[j].GetNeighbor(1)<<endl;
					cout<<cell.Density<<" "<<cell.Pressure<<endl;
					cout<<"A strange error occurred in slope_limit. Probably a nan in cell ";*/
					throw non_finite_slope_limit(cell,centroid_vars);
				}
			}
			res.x[i] *= psi;
			res.y[i] *= psi;
		}
		return res;
	}

	PrimitiveGradient2D slope_limit2(Primitive const& cell,
		Vector2D const& center,
		vector<Primitive> const& neighbors,
		vector<Edge> &edge_list,
		PrimitiveGradient2D const& slope,double diffusecoeff)
	{
		PrimitiveGradient2D res = slope;
		vector<Primitive> centroid_vars(edge_list.size());
		for(int i=0;i<(int)edge_list.size();++i)
			centroid_vars[i] = interp(cell,
			center,
			slope,
			CalcCentroid(edge_list[i]));
		double dphi;
		for(int i=0;i<cell.GetVarNo();++i){
			double psi = 1;
			double maxdiff=0;
			for(int j=0;j<(int)edge_list.size();++j)
				maxdiff=max(maxdiff,abs(neighbors[j][i]-cell[i]));
			for(int j=0;j<(int)edge_list.size();++j)
			{
				dphi=centroid_vars[j][i]-cell[i];
				if(abs(dphi)<0.1*maxdiff&&(centroid_vars[j][i]*cell[i]>0))
					continue;
				if(dphi>0)
					psi = min(psi,diffusecoeff*max((neighbors[j][i]-cell[i])/dphi,0.0));
				else if(dphi<0)
					psi = min(psi,diffusecoeff*max((neighbors[j][i]-cell[i])/dphi,0.0));
				else if(dphi==0)
					psi = min(psi,1.0);
				else
				{
					/*	      cout<<edge_list[j].GetNeighbor(0)<<" "<<edge_list[j].GetNeighbor(1)<<endl;
					cout<<cell.Density<<" "<<cell.Pressure<<endl;
					cout<<"A strange error occurred in slope_limit2. Probably a nan in cell ";*/
					throw non_finite_slope_limit(cell,centroid_vars);
				}
			}
			res.x[i] *= psi;
			res.y[i] *= psi;
		}
		return res;
	}

	vector<int> GetNeighborIndices(int cell_index,
		vector<Edge> const& edge_list)
	{
		vector<int> res(edge_list.size(),0);
		for(int i=0;i<(int)edge_list.size();++i){
			res[i] = GetNeighborIndex(cell_index,edge_list[i]);
		}
		return res;
	}

	vector<Primitive> GetNeighborList(vector<Primitive> const& cells, 
		vector<int> const& cell_indices)
	{
		vector<Primitive> res(cell_indices.size());
		for(int i=0;i<(int)cell_indices.size();++i){
			res[i] = cells[cell_indices[i]];
		}
		return res;
	}

	vector<Vector2D> GetNeighborMeshList(Tessellation const* tess,
		vector<int> const& neighbor_indices)
	{
		vector<Vector2D> res(neighbor_indices.size());
		for(int i=0;i<(int)neighbor_indices.size();++i){
			res[i] = tess->GetMeshPoint(neighbor_indices[i]);
		}
		return res;
	}

	vector<Edge> GetEdgeList(Tessellation const* tess,
		vector<int> const& edge_indices)
	{
		vector<Edge> res(edge_indices.size());
		for(int i=0;i<(int)edge_indices.size();++i){
			res[i] = tess->GetEdge(edge_indices[i]);
		}
		return res;
	}

	/*  double ShockedRatio(Tessellation const* tess,int index,
	vector<Primitive> const& cells)
	{
	vector<int> edges=tess->GetCellEdges(index);
	double DivV=0;
	Edge temp;
	int n0,n1;
	Vector2D n,dv;
	for(int i=0;i<(int)edges.size();++i)
	{
	temp=tess->GetEdge(edges[i]);
	n0=tess->GetOriginalIndex(temp.GetNeighbor(0));
	n1=tess->GetOriginalIndex(temp.GetNeighbor(1));
	if(n0<0)
	{
	n=tess->GetMeshPoint(n1)-temp.GetVertex(0);
	Vector2D p = Parallel(temp);
	n=n-p*ScalarProd(n,p)/pow(abs(p),2);
	n=n/sqrt(n.x*n.x+n.y*n.y);
	dv=cells[n1].Velocity-Reflect(cells[n1].Velocity,temp.GetVertex(1)-
	temp.GetVertex(0));
	}
	else
	if(n1<0)
	{
	n=temp.GetVertex(0)-tess->GetMeshPoint(n0);
	Vector2D p = Parallel(temp);
	n=n-p*ScalarProd(n,p)/pow(abs(p),2);
	n=n/sqrt(n.x*n.x+n.y*n.y);
	dv=Reflect(cells[n0].Velocity,temp.GetVertex(1)-
	temp.GetVertex(0))-cells[n0].Velocity;
	}
	else
	{
	n=(tess->GetMeshPoint(n1)-tess->GetMeshPoint(n0));
	n=n/sqrt(n.x*n.x+n.y*n.y);
	dv=cells[n1].Velocity-cells[n0].Velocity;
	}
	DivV+=ScalarProd(dv,n)*temp.GetLength();
	}
	return DivV/(tess->GetWidth(index)*cells[index].SoundSpeed);
	}
	*/

	double PressureRatio(Primitive cell,vector<Primitive> const& neigh)
	{
		int n=(int)neigh.size();
		double res=1;
		double p=cell.Pressure;
		for(int i=0;i<n;++i)
		{
			if(p>neigh[i].Pressure)
				res=min(res,neigh[i].Pressure/p);
			else
				res=min(res,p/neigh[i].Pressure);
		}
		return res;
	}

	PrimitiveGradient2D CalcSlope
		(Tessellation const* tess,
		vector<Primitive> const& cells,
		int cell_index,
		bool slf,vector<bool> const& mask,
		OuterBoundary const& obc,HydroBoundaryConditions const* hbc,
		double shockratio,double diffusecoeff,double pressure_ratio,
		bool &limited,double time,bool rigidflag)
	{
		vector<int> edge_indices = tess->GetCellEdges(cell_index);
		vector<Edge> edge_list = GetEdgeList(tess,edge_indices);
		if(!mask[cell_index])
			return PrimitiveGradient2D();
		vector<Vector2D> neighbor_mesh_list  = 
			GetNeighborMesh(tess,edge_list,cell_index, obc);
		vector<Primitive> neighbor_list = GetNeighborPrimitive
			(tess,edge_list,cell_index,
			obc,cells,hbc,time);

		PrimitiveGradient2D naive_slope;
		if(IsBoundaryPoint(tess,cell_index)&&rigidflag)
		{
			naive_slope=LinearSlope(cells,tess,cell_index);
			edge_list=RemoveBoundaryEdge(edge_list);
		}
		else
			naive_slope =CalcNaiveSlope
			(cells[cell_index],
			tess->GetMeshPoint(cell_index),
			tess->GetVolume(cell_index),neighbor_list,
			neighbor_mesh_list,edge_list);

		limited=false;
		if(slf)
		{
			if(((naive_slope.x.Velocity.x+naive_slope.y.Velocity.y)*
				tess->GetWidth(cell_index)>(-shockratio*cells[cell_index].
				SoundSpeed)&&(PressureRatio(cells[cell_index],neighbor_list)>
				pressure_ratio)))
			{
				return slope_limit(cells[cell_index],tess->GetCellCM(cell_index),
					neighbor_list,edge_list,naive_slope);
			}
			else
			{
				limited=true;
				return slope_limit2(cells[cell_index],tess->GetCellCM(cell_index),
					neighbor_list,edge_list,naive_slope,diffusecoeff);
			}
		}
		else
			return naive_slope;
	}

	/*  void DisplayDebugInfo(Tessellation const* tess,
	vector<Primitive> const& cells,
	Edge const& edge,
	int side, PrimitiveGradient2D const& slope)
	{
	using namespace std;

	int cell_index = edge.GetNeighbor(side);
	Vector2D cell_cm = tess->GetCellCM(cell_index);
	if(cell_cm.x>0.3&&
	cell_cm.x<0.7&&
	cell_cm.y>0.3&&
	cell_cm.y<0.7){

	cout << endl;
	cout << "Cell index = " << cell_index << endl;
	cout << "cell cm = ( " << cell_cm.x << " , "
	<< cell_cm.y << " ) " << endl;
	cout << "cell density = " << cells[cell_index].Density << endl;
	cout << "cell volume = " << tess->GetVolume(cell_index) << endl;

	vector<int> edge_indices = tess->GetCellEdges(cell_index);
	vector<Edge> edge_list = GetEdgeList(tess,edge_indices);
	cout << "start of edge vertices: " << endl;
	for(int i=0;i<(int)edge_list.size();++i){
	cout << "( " << edge_list[i].GetVertex(0).x << " , "
	<< edge_list[i].GetVertex(0).y << " ) ; ( "
	<< edge_list[i].GetVertex(1).x << " , "
	<< edge_list[i].GetVertex(1).y << " ) " << endl;
	}
	cout << "end of edge vertices" << endl;

	vector<int> neighbor_indices = 
	GetNeighborIndices(cell_index, edge_list);
	vector<Primitive> neighbot_list = 
	GetNeighborList(cells,neighbor_indices);
	vector<Vector2D> neighbor_cm_list = 
	GetNeighborMeshList(tess,neighbor_indices);
	cout <<"neighbor cm and density:" << endl;
	for(int i=0;i<(int)neighbor_indices.size();++i){
	cout << "( " << neighbor_cm_list[i].x << " , "
	<< neighbor_cm_list[i].y << " ) "
	<< neighbot_list[i].Density << endl;
	}
	cout <<"end of neighbor data" << endl;

	cout << "interpolated data at centroids:" << endl;
	for(int i=0;i<(int)edge_list.size();++i){
	Vector2D centroid = CalcCentroid(edge_list[i]);
	Primitive temp = interp(cells[cell_index],
	cell_cm,
	slope,
	CalcCentroid(edge_list[i]));
	cout << "( " << centroid.x << " , "
	<< centroid.y << " ) "
	<< temp.Density << endl;
	}    
	cout << "end interpolated data" << endl;

	throw;
	}    
	}*/
}

void LinearGauss::Prepare
	(Tessellation const* tessellation,
	vector<Primitive> const& cells,
	double /*dt*/,vector<bool> const& mask,double time)
{
	if(tessellation->GetPointNo()!=(int)slopes_.size()){
		slopes_.resize(tessellation->GetPointNo());
		slope_limited_.resize(slopes_.size());
	}
	bool limited = false;
	for(int i=0;i<tessellation->GetPointNo();++i)
	{
		if(!hbc_->IsGhostCell(i,tessellation))
		{
			slopes_[i] = CalcSlope(tessellation,cells,i,slf_,mask,*obc_,
				hbc_,shockratio_,diffusecoeff_,pressure_ratio_,limited,time,
				_rigidflag);
			slope_limited_[i]=limited;
		}
	}
}

namespace {
	Primitive CalcDtFlux
		(Primitive const&cell,
		PrimitiveGradient2D const&grad,double dt)
	{
		Primitive res;
		res.Density-=0.5*dt*(cell.Velocity.x*grad.x.Density+
			cell.Velocity.y*grad.y.Density);
		res.Velocity.x-=0.5*dt*(cell.Velocity.x*
			grad.x.Velocity.x+cell.Velocity.y*
			grad.y.Velocity.x+grad.x.Pressure/cell.Density);
		res.Velocity.y-=0.5*dt*(cell.Velocity.x*
			grad.x.Velocity.y+cell.Velocity.y*
			grad.y.Velocity.y+grad.y.Pressure/cell.Density);
		res.Pressure-=0.5*dt*(cell.Pressure*(1+cell.Pressure/
			(cell.Energy*cell.Density))*(grad.x.Velocity.x+
			grad.y.Velocity.y)+grad.x.Pressure*cell.Velocity.x
			+grad.y.Pressure*cell.Velocity.y);
		return res;
	}
}

Primitive LinearGauss::Interpolate(Tessellation const* tess,
	vector<Primitive> const& cells,double dt,Edge const& edge,int side,
	InterpolationType interptype) const
{
	int cell_index = edge.GetNeighbor(side);
	Vector2D target = CalcCentroid(edge);
	Primitive res;
	if(interptype==InBulk)
	{
		res = interp(cells[cell_index],tess->GetCellCM(cell_index),
			slopes_[cell_index], target);
		if(soitf_)
			res+=CalcDtFlux(cells[cell_index],slopes_[cell_index],dt);
		return res;
	}
	else
		if(interptype==Boundary)
		{
			int other=edge.GetNeighbor((side+1)%2);	
			Primitive res2 = interp(cells[tess->GetOriginalIndex(other)],
				tess->GetMeshPoint(other), slopes_[tess->GetOriginalIndex(other)],
				target);
			if(soitf_)
			{
				res2+=CalcDtFlux(cells[cell_index],slopes_[cell_index],dt);
			}
			return res2;
		}
		else
			throw UniversalError("Wrong interpolation type");
}

PrimitiveGradient2D LinearGauss::GetSlope(int cell_index)const
{
	return slopes_[cell_index];
}

bool LinearGauss::WasSlopeLimited(int index)const
{
	return slope_limited_[index];
}

LinearGauss::~LinearGauss(void) {}
