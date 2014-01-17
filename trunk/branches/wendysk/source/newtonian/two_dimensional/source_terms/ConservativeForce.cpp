#include "ConservativeForce.hpp"


namespace {
Vector2D MassFlux(Tessellation const* tess,int point,
	vector<Conserved> const& fluxes,HydroBoundaryConditions const*hbc)
{
	Vector2D dm;
	vector<int> edge_index=tess->GetCellEdges(point);
	Vector2D center=tess->GetMeshPoint(point);
	int n=(int)edge_index.size();
	Edge edge;
	for(int i=0;i<n;++i)
	{
		edge=tess->GetEdge(edge_index[i]);
		if(point==edge.GetNeighbor(0))
			if(edge.GetNeighbor(1)>-1)
			{
				dm-=edge.GetLength()*fluxes[edge_index[i]].Mass*
					(center-tess->GetMeshPoint(edge.GetNeighbor(1)));
			}
			else
			{
				Vector2D r=hbc->Normal(edge,tess);
				dm-=2*edge.GetLength()*fluxes[edge_index[i]].Mass*r;
			}
		else 
			if(point==edge.GetNeighbor(1))
				if(edge.GetNeighbor(0)>-1)
				{
					dm+=edge.GetLength()*fluxes[edge_index[i]].Mass*
						(center-tess->GetMeshPoint(edge.GetNeighbor(0)));
				}
				else
				{
					Vector2D r=hbc->Normal(edge,tess);
					dm+=2*edge.GetLength()*fluxes[edge_index[i]].Mass*r;
				}
			else
			  throw UniversalError("Error in ConservativeForce MassFlux: Cell and edge are not mutual neighbors");
	}
	return dm;
}
}

ConservativeForce::ConservativeForce(Acceleration *acc,bool DtCalc):acc_(acc),
	DtCalc_(DtCalc),dt_(-1),last_time_(-1){}

ConservativeForce::~ConservativeForce(void){}

Conserved ConservativeForce::Calculate(Tessellation const* tess,
		vector<Primitive> const& cells,int point,vector<Conserved> const& fluxes,
		vector<Vector2D> const& point_velocity,HydroBoundaryConditions const*hbc,
		vector<vector<double> > const &/*tracer_extensive*/,vector<double> &/*dtracer*/,
		double time,double dt)
{
	Conserved res;
	Vector2D acc=acc_->Calculate(tess,cells,point,fluxes,point_velocity,hbc,time,
		dt);
	double volume=tess->GetVolume(point);
	res.Momentum=volume*cells[point].Density*acc;
	res.Energy=cells[point].Density*volume*ScalarProd(point_velocity[point],acc)+
		0.5*ScalarProd(MassFlux(tess,point,fluxes,hbc),acc);
	if(DtCalc_)
	{
		if(time>last_time_)
			dt_=sqrt(volume/abs(acc));
		else
			dt_=max(dt_,sqrt(volume/abs(acc)));
	}
	return res;
}

double ConservativeForce::GetTimeStep(void) const
{
	return dt_;
}

Acceleration::~Acceleration(void) {}
