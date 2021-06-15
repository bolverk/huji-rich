#include "SourceTerm3D.hpp"

SourceTerm3D::~SourceTerm3D(void){}

double SourceTerm3D::SuggestInverseTimeStep(void)const
{
	return 0;
}

void ZeroForce3D::operator()(const Tessellation3D& /*tess*/, const vector<ComputationalCell3D>& /*cells*/,
		const vector<Conserved3D>& /*fluxes*/, const vector<Vector3D>& /*point_velocities*/, const double /*t*/, 
		double /*dt*/, vector<Conserved3D> &/*extensives*/) const {} 
