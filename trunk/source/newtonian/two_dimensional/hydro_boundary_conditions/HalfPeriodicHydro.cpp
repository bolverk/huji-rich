#include "HalfPeriodicHydro.hpp"

Primitive HalfPeriodicHydro::GetBoundaryPrimitive(Edge const& edge,
	Tessellation const& Data,vector<Primitive> const& cells,double time)const
{
	if(obc_.AreWeReflective(edge))
		return rigid.GetBoundaryPrimitive(edge,Data,cells,time);
	else
		return periodic.GetBoundaryPrimitive(edge,Data,cells,time);
}

HalfPeriodicHydro::HalfPeriodicHydro(HalfPeriodicBox const& obc,RiemannSolver &rs)
	:obc_(obc),rigid(rs),periodic(rs){}

HalfPeriodicHydro::~HalfPeriodicHydro(){}

Conserved HalfPeriodicHydro::CalcFlux(Tessellation const& tessellation,
	vector<Primitive> const& cells,Vector2D const& edge_velocity,
	Edge const& edge,SpatialReconstruction const& interp,double dt,double time) const
{
	if(obc_.AreWeReflective(edge))
		return rigid.CalcFlux(tessellation,cells,edge_velocity,edge,interp,
		dt,time);
	else
		return periodic.CalcFlux(tessellation,cells,edge_velocity,edge,
		interp,dt,time);
}

Vector2D HalfPeriodicHydro::CalcEdgeVelocity(Tessellation const& tessellation,
	vector<Vector2D> const& point_velocities,
	Edge const& edge,double time) const
{
	if(obc_.AreWeReflective(edge))
		return rigid.CalcEdgeVelocity(tessellation,point_velocities,edge,time);
	else
		return periodic.CalcEdgeVelocity(tessellation,point_velocities,edge,time);
}

bool HalfPeriodicHydro::IsBoundary(Edge const& edge,Tessellation const& tessellation)const
{
	if((edge.neighbors.first<0)||(edge.neighbors.first>=tessellation.GetPointNo()))
		return true;
	if((edge.neighbors.second<0)||(edge.neighbors.second>=tessellation.GetPointNo()))
		return true;
	return false;
}

bool HalfPeriodicHydro::IsGhostCell(int i,Tessellation const& Data) const
{
	if(i>Data.GetPointNo()||i<0)
		return true;
	else
		return false;
}

vector<double> HalfPeriodicHydro::GetBoundaryTracers(Edge const& edge,Tessellation const& Data,
	vector<vector<double> > const& tracers,double time)const
{
	if(obc_.AreWeReflective(edge))
		return rigid.GetBoundaryTracers(edge,Data,tracers,time);
	else
		return periodic.GetBoundaryTracers(edge,Data,tracers,time);
}

vector<double> HalfPeriodicHydro::CalcTracerFlux
(Tessellation const& tessellation,
 vector<Primitive> const& cells,
 vector<vector<double> > const& tracers,double dm,
 Edge const& edge,int index,double dt,
 double time,SpatialReconstruction const& interp,
 Vector2D const& edge_velocity) const
{
	if(obc_.AreWeReflective(edge))
		return rigid.CalcTracerFlux(tessellation,cells,tracers,dm,edge,
		index,dt,time,interp,edge_velocity);
	else
		return periodic.CalcTracerFlux(tessellation,cells,tracers,dm,edge,
		index,dt,time,interp,edge_velocity);
}
