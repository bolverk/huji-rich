#include "ConservativeForce3D.hpp"

namespace
{
	Vector3D MassFlux(Tessellation3D const& tess, size_t point,vector<Conserved3D> const& fluxes)
	{
		Vector3D dm;
		Vector3D center = tess.GetMeshPoint(point);
		face_vec const& faces = tess.GetCellFaces(point);
		size_t Nfaces = faces.size();

		for (size_t i = 0; i < Nfaces; ++i)
		{
			if (point == tess.GetFaceNeighbors(faces[i]).first)
				dm -= tess.GetArea(faces[i])*fluxes[faces[i]].mass*(center -
					tess.GetMeshPoint(tess.GetFaceNeighbors(faces[i]).second));
			else
				dm += tess.GetArea(faces[i])*fluxes[faces[i]].mass*(center -
					tess.GetMeshPoint(tess.GetFaceNeighbors(faces[i]).first));
		}
		return dm;
	}
}

ConservativeForce3D::ConservativeForce3D(const Acceleration3D& acc,bool mass_flux) : acc_(acc),mass_flux_(mass_flux),dt_(0) {}

ConservativeForce3D::~ConservativeForce3D(void) {}

void ConservativeForce3D::operator()(const Tessellation3D& tess,const vector<ComputationalCell3D>& cells,
	const vector<Conserved3D>& fluxes,const vector<Vector3D>& point_velocities,const double t,double dt,
	TracerStickerNames const& tracerstickernames,vector<Conserved3D> & extensives) const
{
	size_t N = tess.GetPointNo();
	vector<Vector3D> acc;
	acc_(tess, cells, fluxes, t, tracerstickernames, acc);
	dt_ = 0;
	for (size_t i = 0; i < N; ++i)
	{
		dt_ = std::max(dt_, std::sqrt(fastabs(acc[i])/ tess.GetWidth(i)));
		double volume = tess.GetVolume(i);
		double Ek = 0.5 * ScalarProd(extensives[i].momentum, extensives[i].momentum) / extensives[i].mass;
		extensives[i].momentum += volume*cells[i].density*acc[i]*dt;
		if (mass_flux_ && (fastabs(acc[i])*tess.GetWidth(i)*cells[i].density)<(0.5*cells[i].pressure))
		{
			double part0 = volume*cells[i].density*ScalarProd(point_velocities[i], acc[i]);
			double part1 = 0.5*ScalarProd(MassFlux(tess, i, fluxes), acc[i]);
			extensives[i].energy += (part0+part1)*dt;
		}
		else
		{
			double Eknew = 0.5 * ScalarProd(extensives[i].momentum, extensives[i].momentum) / extensives[i].mass;
			//extensives[i].energy += volume * cells[i].density * ScalarProd(acc[i], cells[i].velocity) * dt;
			extensives[i].energy += Eknew - Ek;
		}
	}
}

Acceleration3D::~Acceleration3D(void) {}

double ConservativeForce3D::SuggestInverseTimeStep(void)const
{
	return dt_;
}