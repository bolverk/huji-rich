#include <cmath>
#include <cstddef>
#include "round_cells.hpp"

using std::max;

RoundCells::RoundCells(PointMotion& pm,HydroBoundaryConditions const& hbc,
	double chi,
	double eta,bool coldflows,
	int innerNum,
	OuterBoundary const* outer):
  pm_(pm), hbc_(hbc),
  chi_(chi),
	eta_(eta),
	inner_(innerNum),
	outer_(outer),
	coldflows_(coldflows),
	lastdt_(0),
	evencall_(false),
	external_dt_(-1)
{}

namespace {
	void FixRefinedCells(vector<Vector2D> &vel,Tessellation const& tess,
		HydroBoundaryConditions const& hbc,vector<CustomEvolution*> const& cevolve)
	{
	  for(size_t i=0;i<static_cast<size_t>(tess.GetPointNo());++i)
		{
			if(cevolve[i]!=0)
				continue;
			const vector<int> edge_index(tess.GetCellEdges(static_cast<int>(i)));
			double R=tess.GetWidth(static_cast<int>(i));
			Vector2D r=tess.GetMeshPoint(static_cast<int>(i));
			for(size_t j=0;j<edge_index.size();++j)
			{
				if(DistanceToEdge(r,tess.GetEdge(edge_index[j]))<0.2*R)
				{
				  if(static_cast<int>(i)==tess.GetEdge(edge_index[j]).neighbors.first||
						hbc.IsGhostCell(tess.GetEdge(edge_index[j]).neighbors.second,
						tess))
					{
						if(tess.GetEdge(edge_index[j]).neighbors.second==-1)
							continue;
						int other=tess.GetEdge(edge_index[j]).neighbors.second;
#ifndef RICH_MPI
						other=tess.GetOriginalIndex(other);
#endif
						Vector2D p = Parallel(tess.GetEdge(edge_index[j]));
						p=p/abs(p);
						Vector2D n = r - tess.GetEdge(edge_index[j]).vertices.first;
						n-=ScalarProd(n,p)*p;
						n=n/abs(n);
						double v_avg=0.5*(ScalarProd(vel[i],p)+ScalarProd(vel[static_cast<size_t>(other)],p));
						vel[i]=ScalarProd(vel[i],n)*n+v_avg*p;
						if(other<tess.GetPointNo()&&(cevolve[static_cast<size_t>(other)]==0))
							vel[static_cast<size_t>(other)]=ScalarProd(vel[static_cast<size_t>(other)],n)*n+v_avg*p;
						continue;
					}
				}
			}
		}
	}

	Vector2D calc_dw(int index,
		Tessellation const& tess,
		double c,
		double /*v*/,
		double eta,
		double chi)
	{
		const double R = tess.GetWidth(index);
		const Vector2D& s = tess.GetCellCM(index);
		const Vector2D r = tess.GetMeshPoint(index);
		const double d = abs(s-r);

		if(d<=(0.9*eta*R))
			return Vector2D();
		else
			return chi*c*(s-r)/d;
	}
}

Vector2D RoundCells::CalcVelocity(int index,Tessellation const& tessellation,
	vector<Primitive> const& primitives,double time)
{
	const int nreal=tessellation.GetPointNo();
	if(index<inner_)
		return Vector2D(0,0);
	else
	{
		const Vector2D res = pm_.CalcVelocity
			(index, tessellation, primitives,time);
		if(coldflows_)
			return res;
		Vector2D dw(0,0);
		if(index<nreal)
		  dw = calc_dw(index,tessellation,primitives[static_cast<size_t>(index)].SoundSpeed,
			       abs(primitives[static_cast<size_t>(index)].Velocity),eta_,chi_);
		return res + dw;
	}
}

namespace {
	double numeric_velocity_scale(Tessellation const& tess,int index,double dt,
		vector<Primitive> const& cells)
	{
		const Vector2D& s = tess.GetCellCM(index);
		const Vector2D r = tess.GetMeshPoint(index);
		return max(min(abs(s-r)/dt,3*abs(cells[static_cast<size_t>(index)].Velocity)),cells[static_cast<size_t>(index)].SoundSpeed);
	}
}

vector<Vector2D> RoundCells::calcAllVelocities
	(Tessellation const& tess,
	 vector<Primitive> const& cells,double time,vector<CustomEvolution*> &cevolve, const vector<vector<double> >& /*tracers*/)
{
	vector<Vector2D> res;
	const int n=tess.GetPointNo();
	res.reserve(static_cast<size_t>(n));
	for(size_t i=0;i<static_cast<size_t>(n);++i)
	{
		if(cevolve[i]==0)
		  res.push_back(CalcVelocity(static_cast<int>(i),tess,cells,time));
		else
		  res.push_back(cevolve[i]->CalcVelocity(static_cast<int>(i),tess,cells,time));
	}
#ifdef RICH_MPI
	SendRecvVelocity(res,tess.GetDuplicatedPoints(),tess.GetDuplicatedProcs(),
		tess.GetGhostIndeces(),tess.GetTotalPointNumber());
#endif
	FixRefinedCells(res,tess,hbc_,cevolve);
	return res;
}

namespace {
	void LimitNeighborVelocity(vector<Vector2D> &vel,Tessellation const& tess,
		int index,double factor)
	{
		vector<int> neigh=tess.GetNeighbors(index);
		Vector2D r=tess.GetMeshPoint(index);
		double R=tess.GetWidth(index);
		for(size_t i=0;i<neigh.size();++i)
		{
			if(neigh[i]>-1)
			{
				if(r.distance(tess.GetMeshPoint(neigh[i]))<0.1*R)
				{
				  vel[static_cast<size_t>(neigh[i])]=vel[static_cast<size_t>(neigh[i])]*factor;
					return;
				}
			}
		}
	}
}

void RoundCells::CorrectPointsOverShoot(vector<Vector2D> &v,double dt,
	Tessellation const& tess) const
{
	// check that we don't go outside grid
	int n=tess.GetPointNo();
	const double inv_dt=1.0/dt;
	for(size_t i=static_cast<size_t>(inner_);i<static_cast<size_t>(n);++i)
	{
	  Vector2D point(tess.GetMeshPoint(static_cast<int>(i)));
		if((v[i].x*dt*2+point.x)>outer_->GetGridBoundary(Right))
		{
			double factor=0.9*(outer_->GetGridBoundary(Right)-
				point.x)*inv_dt/abs(v[i]);
			v[i]=v[i]*factor;
			LimitNeighborVelocity(v,tess,static_cast<int>(i),factor);
		}
		if((v[i].x*dt*2+point.x)<outer_->GetGridBoundary(Left))
		{
			double factor=0.9*(point.x-
				outer_->GetGridBoundary(Left))*inv_dt/abs(v[i]);
			v[i]=v[i]*factor;
			LimitNeighborVelocity(v,tess,static_cast<int>(i),factor);
		}
		if((v[i].y*dt*2+point.y)>outer_->GetGridBoundary(Up))
		{
			double factor=0.9*(outer_->GetGridBoundary(Up)-point.y)*
				inv_dt/abs(v[i]);
			v[i]=v[i]*factor;
			LimitNeighborVelocity(v,tess,static_cast<int>(i),factor);
		}
		if((v[i].y*dt*2+point.y)<outer_->GetGridBoundary(Down))
		{
			double factor=0.9*(point.y-outer_->GetGridBoundary(Down))*
				inv_dt/abs(v[i]);
			v[i]=v[i]*factor;
			LimitNeighborVelocity(v,tess,static_cast<int>(i),factor);
		}
	}
	return;
}

void RoundCells::SetExternalTimeStep(double dt)
{
	external_dt_=dt;
}

void RoundCells::ApplyFix(Tessellation const& tess, vector<Primitive> const& cells, double time,
	vector<CustomEvolution*> &cevolve, const vector<vector<double> >& tracers, double dt, vector < Vector2D >
	& velocities)
{
	const int n = tess.GetPointNo();
	if (!evencall_)
	{
		lastdt_ = dt;
		evencall_ = true;
	}
	else
	{
		evencall_ = false;
		if (dt<lastdt_)
			dt = lastdt_;
	}
	if (outer_ != 0 && !coldflows_)
	{
		CorrectPointsOverShoot(velocities, dt, tess);
		return;
	}
	if (coldflows_)
	{
		if (external_dt_>0)
			dt = min(dt, external_dt_);
		for (size_t i = 0; i<static_cast<size_t>(n); ++i)
		{
			if (cevolve[i] != 0)
				continue;
			if (i<static_cast<size_t>(inner_))
			{
				velocities[i].Set(0, 0);
			}
			else
			{
				velocities[i] = pm_.CalcVelocity(static_cast<int>(i), tess, cells, time);
				const double nvs = numeric_velocity_scale(tess, static_cast<int>(i), dt, cells);
				const Vector2D dw = calc_dw(static_cast<int>(i), tess, nvs, nvs, eta_, chi_);
				velocities[i] = velocities[i] + dw;
			}
		}
#ifdef RICH_MPI
		SendRecvVelocity(velocities, tess.GetDuplicatedPoints(), tess.GetDuplicatedProcs(),
			tess.GetGhostIndeces(), tess.GetTotalPointNumber());
#endif
		FixRefinedCells(velocities, tess, hbc_, cevolve);
		if (outer_ != 0)
			CorrectPointsOverShoot(velocities, dt, tess);
	}
}
