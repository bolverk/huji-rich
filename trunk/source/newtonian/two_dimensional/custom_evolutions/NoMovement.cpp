#include "NoMovment.hpp"

bool NoMovement::ShouldForceTracerReset(void)const
{
	return false;
}

NoMovement::NoMovement(void){}

NoMovement::~NoMovement(void){}

Conserved NoMovement::CalcFlux(Tessellation const& tess,
		vector<Primitive> const& cells,	double dt,
		SpatialReconstruction& interpolation,Edge const& edge,
		Vector2D const& face_velocity,
		RiemannSolver const& rs,int index,
		HydroBoundaryConditions const& boundaryconditions,double time,
		vector<vector<double> > const& tracers)
{
	const Vector2D normal_dir = 
			tess.GetMeshPoint(edge.GetNeighbor(1))-
			tess.GetMeshPoint(edge.GetNeighbor(0));

		const Vector2D paral_dir = 
			edge.vertices.second - edge.vertices.first;

		const Primitive left = interpolation.Interpolate
			(tess,cells,dt,edge,0,InBulk,face_velocity);

		const Primitive right = interpolation.Interpolate
			(tess,cells,dt,edge,1,InBulk,face_velocity);

		return FluxInBulk(normal_dir, paral_dir,
			left, right,
			face_velocity, rs);
}

Primitive NoMovement::UpdatePrimitive(vector<Conserved> const& conservedintensive,
		EquationOfState const& eos,vector<Primitive>& /*cells*/,int index,
		Tessellation const& /*tess*/,double /*time*/,
		vector<vector<double> > const& /*tracers*/)
{
	Conserved res = conservedintensive[index];
	return Conserved2Primitive(res, eos);
}

vector<double> NoMovement::UpdateTracer(int index,vector<vector<double> >
		const& tracers,vector<vector<double> > const& tracerchange,
		vector<Primitive> const& cells,Tessellation const& /*tess*/,
		double /*time*/)
{
	int dim=int(tracers[0].size());
	vector<double> res(tracers[index]);
	for(int j=0;j<dim;++j)
		res[j]+=tracerchange[index][j];
	return res;
}

vector<double> NoMovement::CalcTracerFlux(Tessellation const& tess,
		vector<Primitive> const& cells,vector<vector<double> > const& old_trace,
		double dm,Edge const& edge,int cell_index,double dt,double time,
		SpatialReconstruction const& interp,Vector2D const& vface)
{
	const vector<double> temp2=interp.interpolateTracers(tess,cells,old_trace,dt,
		edge,dm<0,InBulk,vface);
	int n=(int)temp2.size();
	vector<double> res(n);
	int sign=(edge.GetNeighbor(1)==cell_index)? 1 : -1;
	for(int i=0;i<n;++i)
		res[i]=sign*dm*dt*edge.GetLength()*temp2[i];
	return res;
}

bool NoMovement::TimeStepRelevant(void)const
{
	return true;
}
