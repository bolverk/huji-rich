#include "AlphaDisc.hpp"

AlphaDisc::AlphaDisc(double alpha, double height_ratio, SpatialReconstruction &grad) : Viscosity(0, grad),alpha_(alpha),
height_ratio_(height_ratio), grads_(grad){}


double AlphaDisc::GetNu(Tessellation const& tess, PhysicalGeometry const& /*pg*/, vector<Primitive> const& cells, int point)const
{
	return abs(tess.GetMeshPoint(point)*height_ratio_*cells[point].SoundSpeed*alpha_);
}