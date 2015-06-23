#include "Viscosity.hpp"

Viscosity::Viscosity(double nu, SpatialReconstruction &grad) : nu_(nu), grads_(grad){}

Conserved Viscosity::Calculate(Tessellation const& tess, const PhysicalGeometry& pg, vector<Primitive> const& cells,
	int point, vector<Conserved> const& /*fluxes*/, vector<Vector2D> const& /*point_velocity*/, HydroBoundaryConditions const& hbc,
	vector<vector<double> > const& /*tracers*/, vector<double>& /*dtracer*/, vector<double> const& /*lengthes*/, double /*t*/, double /*dt*/)
{
	vector<ReducedPrimitiveGradient2D> const& grads = grads_.GetGradients();
	Conserved res;
	vector<int> const& edge_indeces = tess.GetCellEdges(point);
	ReducedPrimitiveGradient2D const& grad = grads[point];
	Primitive const& cell = cells[point];
	Vector2D const& cm = tess.GetCellCM(point);
	const double nu = GetNu(tess, pg, cells, point);
	for (size_t i = 0; i < edge_indeces.size(); ++i)
	{
		Edge const& edge = tess.GetEdge(edge_indeces[i]);
		if (edge.neighbors.first < 0 || edge.neighbors.second < 0) // Fix here for rigid walls
			continue;
		const Vector2D edge_center = 0.5*(edge.vertices.first + edge.vertices.second);
		const int other_org = tess.GetOriginalIndex(edge.neighbors.first) == point ? edge.neighbors.second :
			edge.neighbors.first;
		const int other = tess.GetOriginalIndex(other_org);
		Vector2D const& cm_other = tess.GetCellCM(other) + tess.GetMeshPoint(other_org) - tess.GetMeshPoint(other);
		ReducedPrimitiveGradient2D const& grad_other = grads[other];
		Primitive const& cell_other = cells[other];
		const double nu_other = GetNu(tess, pg, cells, other);
		double eta = 0.5*nu*(cell.Density + ScalarProd(grad.density, edge_center - cm)) +
			0.5*nu_other*(cell_other.Density + ScalarProd(grad_other.density, edge_center - cm_other));
		const double strain_xx = 0.5*eta*(4 * grad.xvelocity.x / 3.0 - 2 * grad.yvelocity.y / 3.0 +
			4 * grad_other.xvelocity.x / 3.0 - 2 * grad_other.yvelocity.y / 3.0);
		const double strain_xy = 0.5*eta*(grad.xvelocity.y + grad.yvelocity.x + grad_other.xvelocity.y + grad_other.yvelocity.x);
		const double strain_yy = 0.5*eta*(4 * grad.yvelocity.y / 3.0 - 2 * grad.xvelocity.x / 3.0 +
			4 * grad_other.yvelocity.y / 3.0 - 2 * grad_other.xvelocity.x / 3.0);
		const Vector2D norm = hbc.Normal(edge, tess);
		const Vector2D momentum_flux = Vector2D(strain_xx*norm.x + strain_xy*norm.y, strain_yy*norm.y + strain_xy*norm.x)*(edge.GetLength() / abs(norm));
		if (other_org == edge.neighbors.first)
		{
			res.Momentum -= momentum_flux;
			res.Energy -= ScalarProd(momentum_flux, 0.5*(cell.Velocity + cell_other.Velocity));
		}
		else
		{
			res.Momentum += momentum_flux;
			res.Energy += ScalarProd(momentum_flux, 0.5*(cell.Velocity + cell_other.Velocity));
		}

	}
	return res;
}

double Viscosity::GetNu(Tessellation const& /*tess*/, PhysicalGeometry const& /*pg*/, vector<Primitive> const& /*cells*/, int /*point*/)const
{
	return nu_;
}

double Viscosity::GetTimeStep(void)const
{
	return 1;
}