#include <limits>
#include <boost/foreach.hpp>
#include "linear_gauss_consistent.hpp"
#include "../../common/hydrodynamics.hpp"

ReducedPrimitiveGradient2D::ReducedPrimitiveGradient2D(void):
density(),
	pressure(),
	xvelocity(),
	yvelocity() {}

ReducedPrimitiveGradient2D::ReducedPrimitiveGradient2D
	(Vector2D const& d,
	Vector2D const& p,
	Vector2D const& vx,
	Vector2D const& vy):
density(d),
	pressure(p),
	xvelocity(vx),
	yvelocity(vy) {}

ReducedPrimitiveGradient2D& ReducedPrimitiveGradient2D::operator+=
	(ReducedPrimitiveGradient2D const& source)
{
	density += source.density;
	pressure += source.pressure;
	xvelocity += source.xvelocity;
	yvelocity += source.yvelocity;
	return *this;
}

ReducedPrimitiveGradient2D& ReducedPrimitiveGradient2D::operator*=
	(double s)
{
	density *= s;
	pressure *= s;
	xvelocity *= s;
	yvelocity *= s;
	return *this;
}

class ReducedPrimitive
{
public:
	ReducedPrimitive(void):
	  density(0),
		  pressure(0),
		  xvelocity(0),
		  yvelocity(0) {}

	  ReducedPrimitive(Primitive const& p):
	  density(p.Density),
		  pressure(p.Pressure),
		  xvelocity(p.Velocity.x),
		  yvelocity(p.Velocity.y) {}

	  ReducedPrimitive(double d,
		  double p,
		  double vx,
		  double vy):
	  density(d),
		  pressure(p),
		  xvelocity(vx),
		  yvelocity(vy) {}

	  double density;
	  double pressure;
	  double xvelocity;
	  double yvelocity;
};

namespace {

	ReducedPrimitive operator*(Vector2D const& v,
		ReducedPrimitiveGradient2D const& rpg)
	{
		return ReducedPrimitive(ScalarProd(v,rpg.density),
			ScalarProd(v,rpg.pressure),
			ScalarProd(v,rpg.xvelocity),
			ScalarProd(v,rpg.yvelocity));
	}

	ReducedPrimitiveGradient2D operator*(double s,
		ReducedPrimitiveGradient2D const& rpg)
	{
		return ReducedPrimitiveGradient2D(s*rpg.density,
			s*rpg.pressure,
			s*rpg.xvelocity,
			s*rpg.yvelocity);
	}

	ReducedPrimitiveGradient2D operator-(ReducedPrimitiveGradient2D const& rpg1,
		ReducedPrimitiveGradient2D const& rpg2)
	{
		return ReducedPrimitiveGradient2D(rpg1.density-rpg2.density,
			rpg1.pressure-rpg2.pressure,
			rpg1.xvelocity-rpg2.xvelocity,
			rpg1.yvelocity-rpg2.yvelocity);
	}

	ReducedPrimitive operator*(double s,
		ReducedPrimitive const& rp)
	{
		return ReducedPrimitive(s*rp.density,
			s*rp.pressure,
			s*rp.xvelocity,
			s*rp.yvelocity);
	}

	ReducedPrimitive operator+(ReducedPrimitive const& rp1,
		ReducedPrimitive const& rp2)
	{
		return ReducedPrimitive(rp1.density+rp2.density,
			rp1.pressure+rp2.pressure,
			rp1.xvelocity+rp2.xvelocity,
			rp1.yvelocity+rp2.yvelocity);
	}

	ReducedPrimitive operator-(ReducedPrimitive const& rp1,
		ReducedPrimitive const& rp2)
	{
		return ReducedPrimitive(rp1.density-rp2.density,
			rp1.pressure-rp2.pressure,
			rp1.xvelocity-rp2.xvelocity,
			rp1.yvelocity-rp2.yvelocity);
	}

	ReducedPrimitiveGradient2D operator*(ReducedPrimitive const& rp,
		Vector2D const& v)
	{
		return ReducedPrimitiveGradient2D(v*rp.density,
			v*rp.pressure,
			v*rp.xvelocity,
			v*rp.yvelocity);
	}

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

LinearGaussConsistent::LinearGaussConsistent
	(EquationOfState const& eos,
	OuterBoundary const& obc,
	HydroBoundaryConditions const* hbc,
	bool slf, bool soitf,double delta_v,double theta,
	double delta_P,bool rigidflag):
eos_(eos),
	slope_limited_(vector<bool>()),
	rslopes_(),
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

	ReducedPrimitive interp(Primitive const& cell,
		Vector2D const& cell_cm,
		ReducedPrimitiveGradient2D const& slope,
		Vector2D const& target)
	{
		return ReducedPrimitive(cell)+(target-cell_cm)*slope;
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

	ReducedPrimitiveGradient2D calc_naive_slope
		(Primitive const& cell,
		Vector2D const& center,
		double cell_volume,
		vector<Primitive> const& neighbors,
		vector<Vector2D> const& neighbor_centers,
		vector<Edge> const& edge_list)
	{
		ReducedPrimitiveGradient2D res;
		const ReducedPrimitive phi_i(cell);
		for(int i=0;i<(int)edge_list.size();++i)
		{
			const Vector2D c_ij = CalcCentroid(edge_list[i])-
				0.5*(neighbor_centers[i]+center);
			const Vector2D r_ij = center - neighbor_centers[i];
			const ReducedPrimitive phi_j(neighbors[i]);
			res += (edge_list[i].GetLength()/cell_volume)*
				((phi_j-phi_i)*(c_ij/abs(r_ij))-
				0.5*(phi_i+phi_j)*(r_ij/abs(r_ij)));
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

	class PropertyGetter
	{
	public:

		virtual double get_property(ReducedPrimitive const& rp) const = 0;

		virtual double get_property(Primitive const& p) const = 0;

		virtual Vector2D& get_property(ReducedPrimitiveGradient2D& rpg) const = 0;

		virtual ~PropertyGetter(void) {}
	};

	class DensityGetter: public PropertyGetter
	{
	public:
		double get_property(ReducedPrimitive const& rp) const
		{
			return rp.density;
		}

		double get_property(Primitive const& p) const
		{
			return p.Density;
		}

		Vector2D& get_property(ReducedPrimitiveGradient2D& rpg) const
		{
			return rpg.density;
		}
	};

	class PressureGetter: public PropertyGetter
	{
	public:
		double get_property(ReducedPrimitive const& rp) const
		{
			return rp.pressure;
		}

		double get_property(Primitive const& p) const
		{
			return p.Pressure;
		}

		Vector2D& get_property(ReducedPrimitiveGradient2D& rpg) const
		{
			return rpg.pressure;
		}
	};

	class XVelocityGetter: public PropertyGetter
	{
	public:

		double get_property(ReducedPrimitive const& rp) const
		{
			return rp.xvelocity;
		}

		double get_property(Primitive const& p) const
		{
			return p.Velocity.x;
		}

		Vector2D& get_property(ReducedPrimitiveGradient2D& rpg) const
		{
			return rpg.xvelocity;
		}
	};

	class YVelocityGetter: public PropertyGetter
	{
	public:

		double get_property(ReducedPrimitive const& rp) const
		{
			return rp.yvelocity;
		}

		double get_property(Primitive const& p) const
		{
			return p.Velocity.y;
		}

		Vector2D& get_property(ReducedPrimitiveGradient2D& rpg) const
		{
			return rpg.yvelocity;
		}
	};

	
	ReducedPrimitiveGradient2D slope_limit
		(Primitive const& cell,
		Vector2D const& center,
		vector<Primitive> const& neighbors,
		vector<Edge> const& edge_list,
		ReducedPrimitiveGradient2D const& slope)
	{
		ReducedPrimitiveGradient2D res = slope;
		vector<ReducedPrimitive> centroid_vars(edge_list.size());
		for(int i=0;i<(int)edge_list.size();++i)
		{
			centroid_vars[i] = interp(cell,
				center,
				slope,
				CalcCentroid(edge_list[i]));
		}
		DensityGetter density_getter;
		PressureGetter pressure_getter;
		XVelocityGetter xvelocity_getter;
		YVelocityGetter yvelocity_getter;
		boost::array<PropertyGetter const*,4> hydro_vars;
		hydro_vars[0]=&density_getter;
		hydro_vars[1]=&pressure_getter;
		hydro_vars[2]=&xvelocity_getter;
		hydro_vars[3]=&yvelocity_getter;
		BOOST_FOREACH(PropertyGetter const* rpg, hydro_vars)
		{
			PropertyGetter const& pg = *rpg;
			double psi = 1;
			boost::container::static_vector<double,17> neighbor_vals;
			BOOST_FOREACH(Primitive const& pr,neighbors)
			{
				neighbor_vals.push_back(pg.get_property(pr));
			}
			const double my_val = pg.get_property(cell);
			//const double min_val = min(vmin2(neighbor_vals),my_val);
			//const double max_val = max(vmax2(neighbor_vals),my_val);
			const double min_val = min(*std::min_element(neighbor_vals.begin(),
				neighbor_vals.end()),my_val);
			const double max_val = max(*std::max_element(neighbor_vals.begin(),
				neighbor_vals.end()),my_val);
			const double max_diff = max(max_val-my_val,
				my_val-min_val);
			for(int i=0;i<(int)edge_list.size();++i)
			{
				const double centroid_val = pg.get_property(centroid_vars[i]);
				const double dphi = centroid_val - my_val;
				if(abs(dphi)<0.1*max_diff&&
					centroid_val*my_val>0)
					continue;
				else if(centroid_val>my_val)
					psi = min(psi,(max_val-my_val)/(centroid_val-my_val));
				else if(centroid_val<my_val)
					psi = min(psi,(min_val-my_val)/(centroid_val-my_val));
				else if(centroid_val==my_val)
					psi = min(psi,1.0);
				else
					throw "Something bad happened in LinearGaussConsistent::Prepare";
			}
			pg.get_property(res) *= psi;
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
					throw non_finite_slope_limit(cell,centroid_vars);
				}
			}
			res.x[i] *= psi;
			res.y[i] *= psi;
		}
		return res;
	}

	ReducedPrimitiveGradient2D shocked_slope_limit
		(Primitive const& cell,
		Vector2D const& center,
		vector<Primitive> const& neighbors,
		vector<Edge> const& edge_list,
		ReducedPrimitiveGradient2D const& slope,
		double diffusecoeff)
	{
		ReducedPrimitiveGradient2D res = slope;
		vector<ReducedPrimitive> centroid_vars(edge_list.size());
		for(int i=0;i<(int)edge_list.size();++i)
		{
			centroid_vars[i] = interp(cell,
				center,
				slope,
				CalcCentroid(edge_list[i]));
		}
		DensityGetter density_getter;
		PressureGetter pressure_getter;
		XVelocityGetter xvelocity_getter;
		YVelocityGetter yvelocity_getter;
		boost::array<PropertyGetter const*,4> hydro_vars;
		hydro_vars[0]=&density_getter;
		hydro_vars[1]=&pressure_getter;
		hydro_vars[2]=&xvelocity_getter;
		hydro_vars[3]=&yvelocity_getter;
		BOOST_FOREACH(PropertyGetter const* rpg, hydro_vars)
		{
			PropertyGetter const& pg = *rpg;
			double psi = 1;
			boost::container::static_vector<double,17> neighbor_vals;
			BOOST_FOREACH(Primitive const& pr,neighbors)
			{
				neighbor_vals.push_back(pg.get_property(pr));
			}
			const double my_val = pg.get_property(cell);
			double maxdiff = 0;
			BOOST_FOREACH(double nv, neighbor_vals)
			{
				maxdiff = max(maxdiff,abs(nv-my_val));
			}
			for(int i=0;i<(int)edge_list.size();++i)
			{
				const double centroid_val = pg.get_property(centroid_vars[i]);
				const double dphi = centroid_val - my_val;
				if(abs(dphi)<0.1*maxdiff&&centroid_val*my_val>0)
					continue;
				else if(dphi>0)
					psi = min(psi,diffusecoeff*max((neighbor_vals[i]-my_val)/dphi,0.0));
				else if(dphi<0)
					psi = min(psi,diffusecoeff*max((neighbor_vals[i]-my_val)/dphi,0.0));
				else if(dphi==0)
					psi = min(psi,1.0);
				else
					throw "Something bad happened in LinearGaussConsistent::Prepare";
			}
			pg.get_property(res) *= psi;
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

	bool is_shock(ReducedPrimitiveGradient2D const& naive_slope,
		double cell_width,
		double shock_ratio,
		Primitive const& cell,
		vector<Primitive> const& neighbor_list,
		double pressure_ratio)
	{
		const bool cond1 = 
			(naive_slope.xvelocity.x+naive_slope.yvelocity.y)*
			cell_width>(-shock_ratio)*cell.SoundSpeed;
		const bool cond2 = PressureRatio(cell,neighbor_list)>pressure_ratio;
		return cond1&&cond2;
	}

	ReducedPrimitiveGradient2D calc_slope
		(Tessellation const* tess,
		vector<Primitive> const& cells,
		int cell_index,
		bool slf,vector<bool> const& mask,
		OuterBoundary const& obc,HydroBoundaryConditions const* hbc,
		double shockratio,double diffusecoeff,double pressure_ratio,
		bool &limited,double time,bool /*rigidflag*/)
	{
		vector<int> edge_indices = tess->GetCellEdges(cell_index);
		vector<Edge> edge_list = GetEdgeList(tess,edge_indices);
		if(!mask[cell_index])
			return ReducedPrimitiveGradient2D();
		vector<Vector2D> neighbor_mesh_list = 
			GetNeighborMesh(tess,edge_list,cell_index, obc);
		vector<Primitive> neighbor_list = GetNeighborPrimitive
			(tess, edge_list, cell_index,
			obc, cells, hbc, time);

		const ReducedPrimitiveGradient2D naive_slope = calc_naive_slope
			(cells[cell_index],
			tess->GetMeshPoint(cell_index),
			tess->GetVolume(cell_index),
			neighbor_list,
			neighbor_mesh_list,
			edge_list);

		if(slf)
		{
			if(is_shock
				(naive_slope,
				tess->GetWidth(cell_index),
				shockratio,
				cells[cell_index],
				neighbor_list,
				pressure_ratio))
			{
				limited = false;
				return slope_limit(cells[cell_index],
					tess->GetCellCM(cell_index),
					neighbor_list,
					edge_list,
					naive_slope);
			}
			else
			{
				limited = true;
				return shocked_slope_limit(cells[cell_index],
					tess->GetCellCM(cell_index),
					neighbor_list,
					edge_list,
					naive_slope,
					diffusecoeff);
			}
		}
		else
		{
			return naive_slope;
		}
	}
}


void LinearGaussConsistent::Prepare
	(Tessellation const* tessellation,
	vector<Primitive> const& cells,
	double /*dt*/,vector<bool> const& mask,double time)
{
	if(tessellation->GetPointNo()!=(int)rslopes_.size()){
		rslopes_.resize(tessellation->GetPointNo());
		slope_limited_.resize(rslopes_.size());
	}
	for(int i=0;i<tessellation->GetPointNo();++i)
	{
		if(!hbc_->IsGhostCell(i,tessellation))
		{
			bool limited = false;
			rslopes_[i] = calc_slope(tessellation,
				cells,
				i,
				slf_,
				mask,*obc_,
				hbc_,shockratio_,
				diffusecoeff_,pressure_ratio_,
				limited,time,
				_rigidflag);
			slope_limited_[i]=limited;
		}
	}
}

namespace {

	Primitive CalcDtFlux
		(Primitive const&cell,
		ReducedPrimitiveGradient2D const& grad,
		double dt)
	{
		Primitive res;
		res.Density-=0.5*dt*(cell.Velocity.x*grad.density.x+
			cell.Velocity.y*grad.density.y);
		res.Velocity.x-=0.5*dt*(cell.Velocity.x*
			grad.xvelocity.x+cell.Velocity.y*
			grad.yvelocity.x+grad.pressure.x/cell.Density);
		res.Velocity.y-=0.5*dt*(cell.Velocity.x*
			grad.xvelocity.y+cell.Velocity.y*
			grad.yvelocity.y+grad.pressure.y/cell.Density);
		res.Pressure-=0.5*dt*(cell.Pressure*(1+cell.Pressure/
			(cell.Energy*cell.Density))*(grad.xvelocity.x+
			grad.yvelocity.y)+grad.pressure.x*cell.Velocity.x
			+grad.pressure.y*cell.Velocity.y);
		return res;
	}

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

Primitive LinearGaussConsistent::Interpolate(Tessellation const* tess,
	vector<Primitive> const& cells,double dt,Edge const& edge,int side,
	InterpolationType interptype) const
{
	int cell_index = edge.GetNeighbor(side);
	Vector2D target = CalcCentroid(edge);
	if(interptype==InBulk)
	{
		const ReducedPrimitive temp = interp
			(cells[cell_index],
			tess->GetCellCM(cell_index),
			rslopes_[cell_index],
			target);
		Primitive res = CalcPrimitive(temp.density,
			temp.pressure,
			Vector2D(temp.xvelocity,
			temp.yvelocity),
			eos_);
		if(soitf_){
			res+=CalcDtFlux(cells[cell_index],rslopes_[cell_index],dt);
			res = CalcPrimitive(res.Density,
				res.Pressure,
				Vector2D(res.Velocity.x,
				res.Velocity.y),
				eos_);
		}
		return res;
	}
	else
		if(interptype==Boundary)
		{
			const int other=edge.GetNeighbor((side+1)%2);
			const ReducedPrimitive temp = interp
				(cells[tess->GetOriginalIndex(other)],
				tess->GetMeshPoint(other),
				rslopes_[tess->GetOriginalIndex(other)],
				target);
			Primitive res = CalcPrimitive(temp.density,
				temp.pressure,
				Vector2D(temp.xvelocity,
				temp.yvelocity),
				eos_);
			/*
			Primitive res2 = interp(cells[tess->GetOriginalIndex(other)],
			tess->GetMeshPoint(other), slopes_[tess->GetOriginalIndex(other)],
			target);
			*/
			if(soitf_)
			{
				res+=CalcDtFlux(cells[cell_index],rslopes_[cell_index],dt);
				res = CalcPrimitive(res.Density,
					res.Pressure,
					res.Velocity,
					eos_);
			}
			return res;
		}
		else
			throw UniversalError("Wrong interpolation type");
}

bool LinearGaussConsistent::WasSlopeLimited(int index)const
{
	return slope_limited_[index];
}

LinearGaussConsistent::~LinearGaussConsistent(void) {}
