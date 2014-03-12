#include "CustomOuter.hpp"
#include "../../../misc/universal_error.hpp"

vector<double> CustomOuter::GetBoundaryTracers(Edge const& edge,Tessellation const& Data,
	vector<vector<double> > const& tracers,double time)const
{
	Vector2D n = Normal(edge,Data);
	double angle;
	if(IsGhostCell(edge.GetNeighbor(0),Data))
		angle=atan2(-n.y,-n.x);
	else
		angle=atan2(n.y,n.x);
	if(angle>-0.25*M_PI&&angle<0.25*M_PI)
		return _right.GetBoundaryTracers(edge,Data,tracers,time);
	if(angle<0.75*M_PI&&angle>0.25*M_PI)
		return _up.GetBoundaryTracers(edge,Data,tracers,time);
	if(abs(angle)>0.75*M_PI)
		return _left.GetBoundaryTracers(edge,Data,tracers,time);
	if(angle>-0.75*M_PI&&angle-0.25*M_PI)
		return _down.GetBoundaryTracers(edge,Data,tracers,time);
	throw UniversalError("error in custom Outer conditions GetTracer");
}

Primitive CustomOuter::GetBoundaryPrimitive(Edge const& edge,
	Tessellation const& tessellation,vector<Primitive> const& cells,
	double time)const
{
	Vector2D n = Normal(edge, tessellation);
	double angle;
	if(IsGhostCell(edge.GetNeighbor(0),tessellation))
		angle=atan2(-n.y,-n.x);
	else
		angle=atan2(n.y,n.x);
	if(angle>=-0.25*M_PI&&angle<0.25*M_PI)
		return _right.GetBoundaryPrimitive(edge,tessellation,cells,time);
	if(angle<0.75*M_PI&&angle>=0.25*M_PI)
		return _up.GetBoundaryPrimitive(edge,tessellation,cells,time);
	if(abs(angle)>0.75*M_PI)
		return _left.GetBoundaryPrimitive(edge,tessellation,cells,time);
	if(angle>=-0.75*M_PI&&angle<-0.25*M_PI)
		return _down.GetBoundaryPrimitive(edge,tessellation,cells,time);
	throw UniversalError("error in custom Outer conditions GetPrimitive");
}

CustomOuter::CustomOuter(HydroBoundaryConditions& Left,HydroBoundaryConditions& Right,
	HydroBoundaryConditions& Down,HydroBoundaryConditions& Up)
	:_left(Left),_right(Right),_down(Down),_up(Up){}

Conserved CustomOuter::CalcFlux
	(Tessellation const& tessellation,
	vector<Primitive> const& cells,
	Vector2D const& edge_velocity,
	Edge const& edge,SpatialReconstruction const& interp,double dt,
	double time) const
{
	Vector2D n = Normal(edge, tessellation);
	double angle;
	if(IsGhostCell(edge.GetNeighbor(0),tessellation))
		angle=atan2(-n.y,-n.x);
	else
		angle=atan2(n.y,n.x);
	if(angle>-0.25*M_PI&&angle<0.25*M_PI)
	  return _right.CalcFlux(tessellation,cells,edge_velocity,edge,
				  interp,dt,time);
	if(angle<0.75*M_PI&&angle>0.25*M_PI)
	  return _up.CalcFlux(tessellation,cells,edge_velocity,edge,
			       interp,dt,time);
	if(abs(angle)>0.75*M_PI)
	  return _left.CalcFlux(tessellation,cells,edge_velocity,edge,
				 interp,dt,time);
	if(angle>-0.75*M_PI&&angle-0.25*M_PI)
	  return _down.CalcFlux(tessellation,cells,edge_velocity,edge,
				 interp,dt,time);
	throw UniversalError("error in custom Outer conditions");
}

Vector2D CustomOuter::CalcEdgeVelocity
	(Tessellation const& /*tessellation*/,
	vector<Vector2D> const& /*point_velocities*/,
	Edge const& /*edge*/, double /*time*/) const
{
	Vector2D res(0,0);
	return res;
}

bool CustomOuter::IsBoundary(Edge const& edge,Tessellation const& tessellation)const
{
	if((edge.GetNeighbor(0)<0)||(edge.GetNeighbor(0)>=tessellation.GetPointNo()))
		return true;
	if((edge.GetNeighbor(1)<0)||(edge.GetNeighbor(1)>=tessellation.GetPointNo()))
		return true;
	return false;
}

bool CustomOuter::IsGhostCell(int i,Tessellation const& Data) const
{
	if(i>Data.GetPointNo()||i<0)
		return true;
	else
		return false;
}

vector<double> CustomOuter::CalcTracerFlux
(Tessellation const& tessellation,vector<Primitive> const& cells,
 vector<vector<double> > const& tracers,double dm,
 Edge const& edge,int index,double dt,
 double time,SpatialReconstruction const& interp,
 Vector2D const& edge_velocity) const
{
	Vector2D n = Normal(edge, tessellation);
	double angle;
	if(IsGhostCell(edge.GetNeighbor(0),tessellation))
		angle=atan2(-n.y,-n.x);
	else
		angle=atan2(n.y,n.x);
	if(angle>-0.25*M_PI&&angle<0.25*M_PI)
		return _right.CalcTracerFlux(tessellation,cells,tracers,dm,edge,index,
		dt,time,interp,edge_velocity);
	if(angle<0.75*M_PI&&angle>0.25*M_PI)
		return _up.CalcTracerFlux(tessellation,cells,tracers,dm,edge,index,
		dt,time,interp,edge_velocity);
	if(abs(angle)>0.75*M_PI)
		return _left.CalcTracerFlux(tessellation,cells,tracers,dm,edge,index,
		dt,time,interp,edge_velocity);
	if(angle>-0.75*M_PI&&angle-0.25*M_PI)
		return _down.CalcTracerFlux(tessellation,cells,tracers,dm,edge,index,
		dt,time,interp,edge_velocity);
	throw UniversalError("error in custom Outer conditions");
}
