#include "NoMovment.hpp"

bool NoMovement::DensityFloorRelevant(void) const
{
	return true;
}


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
			       RiemannSolver const& rs,int /*index*/,
			       HydroBoundaryConditions const& /*boundaryconditions*/,double /*time*/,
			       vector<vector<double> > const& /*tracers*/)
{
	const Vector2D normal_dir =
			tess.GetMeshPoint(edge.neighbors.second)-
			tess.GetMeshPoint(edge.neighbors.first);

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
  Conserved res = conservedintensive[static_cast<size_t>(index)];
	return Conserved2Primitive(res, eos);
}

vector<double> NoMovement::UpdateTracer(int index,vector<vector<double> >
		const& tracers,vector<vector<double> > const& tracerchange,
					vector<Primitive> const& /*cells*/,Tessellation const& /*tess*/,
		double /*time*/)
{
	int dim=int(tracers[0].size());
	vector<double> res(tracers[static_cast<size_t>(index)]);
	for(int j=0;j<dim;++j)
	  res[static_cast<size_t>(j)]+=tracerchange[static_cast<size_t>(index)][static_cast<size_t>(j)];
	return res;
}

vector<double> NoMovement::CalcTracerFlux(Tessellation const& tess,
		vector<Primitive> const& cells,vector<vector<double> > const& old_trace,
					  double dm,Edge const& edge,int /*cell_index*/,double dt,double /*time*/,
		SpatialReconstruction const& interp,Vector2D const& vface)
{
	const vector<double> temp2=interp.interpolateTracers(tess,cells,old_trace,dt,
		edge,dm<0,InBulk,vface);
	int n=static_cast<int>(temp2.size());
	vector<double> res(static_cast<size_t>(n));
	for(int i=0;i<n;++i)
	  res[static_cast<size_t>(i)]=dm*dt*edge.GetLength()*temp2[static_cast<size_t>(i)];
	return res;
}

bool NoMovement::TimeStepRelevant(void)const
{
	return true;
}

bool NoMovement::isRelevantToInterpolation(void)const
{
	return true;
}
