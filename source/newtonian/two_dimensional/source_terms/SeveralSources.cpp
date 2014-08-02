#include "SeveralSources.hpp"


SeveralSources::SeveralSources(vector<SourceTerm*> forces):
sources_(vector<SourceTerm*> ())
{
	for(int i=0;i<(int)forces.size();++i)
		sources_.push_back(forces[i]);
}


SeveralSources::~SeveralSources(void)
{
}


Conserved SeveralSources::Calculate(Tessellation const& tess,
	vector<Primitive> const& cells,int point,
	vector<Conserved> const& fluxes,
	vector<Vector2D> const& point_velocity,
	HydroBoundaryConditions const& hbc,
	vector<vector<double> > const &tracer_extensive,vector<double> &dtracer,
	vector<double> const& lengthes,
	double time,double dt)
{
	Conserved res;
	vector<double> dtracer_local;
	if(!dtracer.empty())
		dtracer_local.assign(dtracer.size(),0);
	for(int i=0;i<(int)sources_.size();++i)
	{
		res+=sources_[i]->Calculate(tess,cells,point,fluxes,
		point_velocity,hbc,tracer_extensive,dtracer_local,lengthes,time,dt);
		if(!dtracer.empty())
			dtracer=dtracer+dtracer_local;
	}
	return res;
}
