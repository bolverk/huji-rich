#include "ConservativeForce.hpp"
#include "../../../misc/lazy_list.hpp"

using std::min;
using std::sqrt;

namespace
{
	Vector2D MassFlux(Tessellation const& tess, int point,
		vector<Extensive> const& fluxes)
	{
		Vector2D dm;
		vector<int> edge_index = tess.GetCellEdges(point);
		//Vector2D center = tess.GetCellCM(point);
		Vector2D center = tess.GetMeshPoint(point);
		int n = static_cast<int>(edge_index.size());
		Edge edge;
		for (int i = 0; i < n; ++i)
		{
			edge = tess.GetEdge(edge_index[static_cast<size_t>(i)]);
			if (point == edge.neighbors.first)
			{
				dm -= edge.GetLength() * fluxes[static_cast<size_t>(edge_index[static_cast<size_t>(i)])].mass *
					//(center - 0.5*(edge.vertices.first + edge.vertices.second));
					(center - tess.GetMeshPoint(edge.neighbors.second));
			}
			else
				if (point == edge.neighbors.second)
				{
					dm += edge.GetLength() * fluxes[static_cast<size_t>(edge_index[static_cast<size_t>(i)])].mass *
						//(center - 0.5*(edge.vertices.first + edge.vertices.second));
						(center - tess.GetMeshPoint(edge.neighbors.first));
				}
				else
					throw UniversalError("Error in ConservativeForce MassFlux: Cell and edge are not mutual neighbors");
		}
		return dm;
	}
}

ConservativeForce::ConservativeForce(const Acceleration& acc, bool mass_flux) :
	acc_(acc), mass_flux_(mass_flux) {}

ConservativeForce::~ConservativeForce(void) {}

vector<Extensive> ConservativeForce::operator()
(const Tessellation& tess,
	const PhysicalGeometry& /*pg*/,
	const CacheData& cd,
	const vector<ComputationalCell>& cells,
	const vector<Extensive>& fluxes,
	const vector<Vector2D>& point_velocities,
	const double t,
	TracerStickerNames const& tracerstickernames) const
{
	vector<Extensive> res(static_cast<size_t>(tess.GetPointNo()));
	for (size_t i = 0; i < res.size(); ++i)
	{
		const Vector2D acc = acc_(tess, cells, fluxes, t, static_cast<int>(i),tracerstickernames);
		const double volume = cd.volumes[i];
		res[i].mass = 0;
		res[i].momentum = volume*cells[i].density*acc;
		if (!mass_flux_)
			res[i].energy = volume*cells[i].density*ScalarProd(acc, cells[i].velocity);
		else
		{
			const Vector2D mass_flux = MassFlux(tess, static_cast<int>(i), fluxes);
			res[i].energy = volume*cells[i].density*ScalarProd(point_velocities[i], acc) + 
				0.5*ScalarProd(mass_flux, acc);
		}
//		res[i].tracers.resize(tracerstickernames.tracer_names.size(),0);
	}
	return res;
}

Acceleration::~Acceleration(void) {}
