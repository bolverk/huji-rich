#include <limits>
#include <boost/foreach.hpp>
#include "LinearGaussImproved.hpp"
#include "../../common/hydrodynamics.hpp"

using std::min;
using std::max;
using std::abs;

namespace {
	class ReducedPrimitive
	{
	public:
		ReducedPrimitive(void) :
			density(0), pressure(0), xvelocity(0), yvelocity(0), tracers(vector<double>()) {}

		ReducedPrimitive(Primitive const& p, vector<double> const& trace
			= vector<double>()) :
			density(p.Density), pressure(p.Pressure), xvelocity(p.Velocity.x),
			yvelocity(p.Velocity.y), tracers(trace) {}

		ReducedPrimitive(double d, double p, double vx, double vy,
			vector<double> const& trace = vector<double>()) :
			density(d), pressure(p), xvelocity(vx), yvelocity(vy), tracers(trace) {}

		double density;
		double pressure;
		double xvelocity;
		double yvelocity;
		vector<double> tracers;
	};
}

namespace
{
	ReducedPrimitive AddRp(ReducedPrimitive const& rp1,
		ReducedPrimitive const& rp2)
	{
		return ReducedPrimitive(rp1.density + rp2.density, rp1.pressure + rp2.pressure,
			rp1.xvelocity + rp2.xvelocity, rp1.yvelocity + rp2.yvelocity,
			rp1.tracers + rp2.tracers);
	}

  /*
	ReducedPrimitive MinusRp(ReducedPrimitive const& rp1,
		ReducedPrimitive const& rp2)
	{
		return ReducedPrimitive(rp1.density - rp2.density, rp1.pressure - rp2.pressure,
			rp1.xvelocity - rp2.xvelocity, rp1.yvelocity - rp2.yvelocity,
			rp1.tracers - rp2.tracers);
	}
  */

	vector<double> ScalarProd(Vector2D const&v, vector<Vector2D> const& vec)
	{
		vector<double> res(vec.size());
		for (size_t i = 0; i<vec.size(); ++i)
			res[i] = ScalarProd(v, vec[i]);
		return res;
	}

	ReducedPrimitive operator*(Vector2D const& v,
		ReducedPrimitiveGradient2D const& rpg)
	{
		return ReducedPrimitive(ScalarProd(v, rpg.density),
			ScalarProd(v, rpg.pressure),
			ScalarProd(v, rpg.xvelocity),
			ScalarProd(v, rpg.yvelocity), ScalarProd(v, rpg.tracers));
	}

  /*
	ReducedPrimitiveGradient2D operator-(ReducedPrimitiveGradient2D const& rpg1,
		ReducedPrimitiveGradient2D const& rpg2)
	{
		return ReducedPrimitiveGradient2D(rpg1.density - rpg2.density, rpg1.pressure
			- rpg2.pressure, rpg1.xvelocity - rpg2.xvelocity, rpg1.yvelocity
			- rpg2.yvelocity, rpg1.tracers - rpg2.tracers);
	}
  */

	vector<Vector2D> operator*(Vector2D const&v, vector<double> const& vec)
	{
		if (vec.empty())
			return vector<Vector2D>();
		vector<Vector2D> res(vec.size());
		for (size_t i = 0; i<vec.size(); ++i)
			res[i] = v*vec[i];
		return res;
	}

	ReducedPrimitiveGradient2D operator*(ReducedPrimitive const& rp,
		Vector2D const& v)
	{
		return ReducedPrimitiveGradient2D(v*rp.density, v*rp.pressure, v*rp.xvelocity,
			v*rp.yvelocity, v*rp.tracers);
	}
}

LinearGaussImproved::LinearGaussImproved
(EquationOfState const& eos,
OuterBoundary const& obc,
HydroBoundaryConditions const& hbc,
bool slf, bool soitf, double delta_v, double theta,
double delta_P, bool rigidflag) :
eos_(eos), rslopes_(), obc_(obc), hbc_(hbc), slf_(slf), soitf_(soitf),
shockratio_(delta_v), diffusecoeff_(theta), pressure_ratio_(delta_P),
_rigidflag(rigidflag) {}

namespace
{
	ReducedPrimitive interp_all(Primitive const& cell, vector<double> const&
		cell_tracer, Vector2D const& cell_cm, ReducedPrimitiveGradient2D const& slope,
		Vector2D const& target)
	{
		return AddRp(ReducedPrimitive(cell, cell_tracer), (target - cell_cm)*slope);
	}

	Primitive interp_primitive(Primitive const& cell, Vector2D const& cell_cm,
		ReducedPrimitiveGradient2D const& slope, Vector2D const& target)
	{
		Primitive res(cell);
		res.Pressure += ScalarProd(target - cell_cm, slope.pressure);
		res.Density += ScalarProd(target - cell_cm, slope.density);
		res.Velocity.x += ScalarProd(target - cell_cm, slope.xvelocity);
		res.Velocity.y += ScalarProd(target - cell_cm, slope.yvelocity);
		return res;
	}

	vector<double> operator*(Vector2D const&v, vector<Vector2D> const&vec)
	{
		vector<double> res(vec.size());
		for (size_t i = 0; i<vec.size(); ++i)
			res[i] = ScalarProd(v, vec[i]);
		return res;
	}

	vector<double> interp_tracer(vector<double> const& cell,
		Vector2D const& cell_cm, ReducedPrimitiveGradient2D const& slope,
		Vector2D const& target)
	{
		return cell + (target - cell_cm)*slope.tracers;
	}

	Vector2D GetReflectedPoint(Tessellation const& tess, int point,
		OuterBoundary const& /*obc*/, Edge const& edge)
	{
		Vector2D par = edge.vertices.second - edge.vertices.first;
		par = par / abs(par);
		Vector2D norm = Vector2D(-par.y, par.x);
		Vector2D tofix = tess.GetMeshPoint(point) - edge.vertices.first;
		tofix -= 2 * ScalarProd(norm, tofix)*norm - edge.vertices.first;
		return tofix;
	}

	vector<Vector2D> GetNeighborMesh
		(Tessellation const& tess, vector<Edge> const&
		edges, int cell_index, OuterBoundary const& obc)
	{
		vector<Vector2D> res(edges.size());
		for (size_t i = 0; i<edges.size(); ++i)
		{
			const int neigh0 = edges[i].neighbors.first;
			const int neigh1 = edges[i].neighbors.second;
			if (neigh0 == cell_index)
				if (neigh1>-1) // we are not near rigid wall
					res[i] = tess.GetMeshPoint(neigh1);
				else
					res[i] = GetReflectedPoint(tess, neigh0, obc, edges[i]);
			else
				if (neigh0>-1)
					res[i] = tess.GetMeshPoint(neigh0);
				else
					res[i] = GetReflectedPoint(tess, neigh1, obc, edges[i]);
		}
		return res;
	}

	vector<Vector2D> GetNeighborCM
		(Tessellation const& tess, vector<Edge> const&
		edges, int cell_index, OuterBoundary const& obc)
	{
		vector<Vector2D> res(edges.size());
		for (size_t i = 0; i<edges.size(); ++i)
		{
			const int neigh0 = edges[i].neighbors.first;
			const int neigh1 = edges[i].neighbors.second;
			if (neigh0 == cell_index)
				if (neigh1>-1) // we are not near rigid wall
					res[i] = tess.GetCellCM(neigh1);
				else
					res[i] = tess.GetCellCM(neigh0) + GetReflectedPoint(tess, neigh0, obc, edges[i]) - tess.GetMeshPoint(cell_index);
			else
				if (neigh0>-1)
					res[i] = tess.GetCellCM(neigh0);
				else
					res[i] = tess.GetCellCM(neigh1) + GetReflectedPoint(tess, neigh1, obc, edges[i]) - tess.GetMeshPoint(cell_index);
		}
		return res;
	}

	vector<Primitive> GetNeighborPrimitive
		(Tessellation const& tess, vector<Edge> const&
		edges, int cell_index, OuterBoundary const& /*obc*/,
		vector<Primitive> const& cells, HydroBoundaryConditions const& hbc,
		double time, vector<bool> const& isrelevant)
	{
		int neigh0, neigh1;
		vector<Primitive> res(edges.size());
		for (size_t i = 0; i<edges.size(); ++i)
		{
			if (hbc.IsBoundary(edges[i], tess))
			{
				res[i] = hbc.GetBoundaryPrimitive(edges[i], tess, cells, time);
			}
			else
			{
				neigh0 = edges[i].neighbors.first;
				neigh1 = edges[i].neighbors.second;
				if (neigh0 == cell_index)
					if (isrelevant[static_cast<size_t>(neigh1)])
						res[i] = cells[static_cast<size_t>(neigh1)];
					else
						res[i] = cells[static_cast<size_t>(neigh0)];
				else
					if (isrelevant[static_cast<size_t>(neigh0)])
						res[i] = cells[static_cast<size_t>(neigh0)];
					else
						res[i] = cells[static_cast<size_t>(neigh1)];
			}
		}
		return res;
	}

	vector<vector<double> > GetNeighborTracers(Tessellation const& tess,
		vector<Edge> const&	edges, int cell_index,
		vector<vector<double> > const& tracers,
		HydroBoundaryConditions const& hbc, double time)
	{
		int neigh0, neigh1;
		vector<vector<double> > res(edges.size());
		for (size_t i = 0; i<edges.size(); ++i)
		{
			if (hbc.IsBoundary(edges[i], tess))
			{
				res[i] = hbc.GetBoundaryTracers(edges[i], tess, tracers, time);
			}
			else
			{
				neigh0 = edges[i].neighbors.first;
				neigh1 = edges[i].neighbors.second;
				if (neigh0 == cell_index)
					res[i] = tracers[static_cast<size_t>(neigh1)];
				else
					res[i] = tracers[static_cast<size_t>(neigh0)];
			}
		}
		return res;
	}

	ReducedPrimitiveGradient2D calc_naive_slope(Primitive const& cell,
		Vector2D const& center, Vector2D const& cell_cm, double cell_volume, vector<Primitive> const& neighbors,
		vector<Vector2D> const& neighbor_centers, vector<Vector2D> const& neigh_cm, vector<Edge> const& edge_list,
		vector<double> const& cell_tracer = vector<double>(), vector<vector<double> > const& neighbor_tracers
		= vector<vector<double> >())
	{
		ReducedPrimitiveGradient2D res;
		const ReducedPrimitive phi_i(cell, cell_tracer);
		int n = static_cast<int>(edge_list.size());
		if (n>20)
		{
			UniversalError eo("Cell has too many neighbors");
			eo.AddEntry("Cell x cor", center.x);
			eo.AddEntry("Cell y cor", center.y);
			throw eo;
		}
		// Create the matrix to invert and the vector to compare
		vector<double> m(4, 0);
		ReducedPrimitiveGradient2D vec_compare;
		for (size_t i = 0; i < edge_list.size(); ++i)
		{
			const Vector2D c_ij = CalcCentroid(edge_list[i]) -
				0.5*(neigh_cm[i] + cell_cm);
			Vector2D r_ij = neighbor_centers[i] - center;
			const double e_length = edge_list[i].GetLength();
			r_ij = r_ij*e_length / abs(r_ij);
			m[0] -= c_ij.x*r_ij.x;
			m[1] -= c_ij.y*r_ij.x;
			m[2] -= c_ij.x*r_ij.y;
			m[3] -= c_ij.y*r_ij.y;

			ReducedPrimitive phi_j;
			if (cell_tracer.empty())
				phi_j = ReducedPrimitive(neighbors[i]);
			else
				phi_j = ReducedPrimitive(neighbors[i], neighbor_tracers[i]);
			r_ij *= 0.5;
			vec_compare += AddRp(phi_i, phi_j)*r_ij;
		}
		m[0] += cell_volume;
		m[3] += cell_volume;
		// Find the det
		const double det = m[0] * m[3] - m[1] * m[2];
		// Check none singular
		if (abs(det) < 1e-10*cell_volume*cell_volume)
		{
			UniversalError eo("Singular matrix");
			eo.AddEntry("Cell x cor", center.x);
			eo.AddEntry("Cell y cor", center.y);
			eo.AddEntry("Cell volume", cell_volume);
			eo.AddEntry("Det was", det);
			throw eo;
		}
		// Invert the matrix
		vector<double> m_inv(4);
		const double det_inv = 1.0 / det;
		m_inv[0] = m[3] * det_inv;
		m_inv[1] = -m[1] * det_inv;
		m_inv[2] = -m[2] * det_inv;
		m_inv[3] = m[0] * det_inv;
		// Calculate the gradient
		res.density.x = m_inv[0] * vec_compare.density.x;
		res.density.x += m_inv[1] * vec_compare.density.y;
		res.density.y = m_inv[2] * vec_compare.density.x;
		res.density.y += m_inv[3] * vec_compare.density.y;

		res.pressure.x = m_inv[0] * vec_compare.pressure.x;
		res.pressure.x += m_inv[1] * vec_compare.pressure.y;
		res.pressure.y = m_inv[2] * vec_compare.pressure.x;
		res.pressure.y += m_inv[3] * vec_compare.pressure.y;

		res.xvelocity.x = m_inv[0] * vec_compare.xvelocity.x;
		res.xvelocity.x += m_inv[1] * vec_compare.xvelocity.y;
		res.xvelocity.y = m_inv[2] * vec_compare.xvelocity.x;
		res.xvelocity.y += m_inv[3] * vec_compare.xvelocity.y;

		res.yvelocity.x = m_inv[0] * vec_compare.yvelocity.x;
		res.yvelocity.x += m_inv[1] * vec_compare.yvelocity.y;
		res.yvelocity.y = m_inv[2] * vec_compare.yvelocity.x;
		res.yvelocity.y += m_inv[3] * vec_compare.yvelocity.y;

		if (!cell_tracer.empty())
		{
			for (size_t j = 0; j < cell_tracer.size(); ++j)
			{
				res.tracers[j].x = m_inv[0] * vec_compare.tracers[j].x;
				res.tracers[j].x += m_inv[1] * vec_compare.tracers[j].y;
				res.tracers[j].y = m_inv[2] * vec_compare.tracers[j].x;
				res.tracers[j].y += m_inv[3] * vec_compare.tracers[j].y;
			}
		}
		return res;
	}


	class PropertyGetter
	{
	public:

		virtual double get_property(ReducedPrimitive const& rp) const = 0;

		virtual double get_property(Primitive const& p) const = 0;

		virtual Vector2D& get_property(ReducedPrimitiveGradient2D& rpg) const = 0;

		virtual ~PropertyGetter(void) {}
	};

	class DensityGetter : public PropertyGetter
	{
	public:
		double get_property(ReducedPrimitive const& rp) const
		{
			return rp.density;
		}

		double get_property(Primitive const& p) const
		{
			return p.Density;
		}

		Vector2D& get_property(ReducedPrimitiveGradient2D& rpg) const
		{
			return rpg.density;
		}
	};

	class PressureGetter : public PropertyGetter
	{
	public:
		double get_property(ReducedPrimitive const& rp) const
		{
			return rp.pressure;
		}

		double get_property(Primitive const& p) const
		{
			return p.Pressure;
		}

		Vector2D& get_property(ReducedPrimitiveGradient2D& rpg) const
		{
			return rpg.pressure;
		}
	};

	class XVelocityGetter : public PropertyGetter
	{
	public:

		double get_property(ReducedPrimitive const& rp) const
		{
			return rp.xvelocity;
		}

		double get_property(Primitive const& p) const
		{
			return p.Velocity.x;
		}

		Vector2D& get_property(ReducedPrimitiveGradient2D& rpg) const
		{
			return rpg.xvelocity;
		}
	};

	class YVelocityGetter : public PropertyGetter
	{
	public:

		double get_property(ReducedPrimitive const& rp) const
		{
			return rp.yvelocity;
		}

		double get_property(Primitive const& p) const
		{
			return p.Velocity.y;
		}

		Vector2D& get_property(ReducedPrimitiveGradient2D& rpg) const
		{
			return rpg.yvelocity;
		}
	};

	ReducedPrimitiveGradient2D slope_limit(Primitive const& cell,
		Vector2D const& center, vector<Primitive> const& neighbors,
		vector<Edge> const& edge_list, ReducedPrimitiveGradient2D const& slope,
		vector<vector<double> > const& neighbor_tracers = vector<vector<double> >(),
		vector<double> const& cell_trace = vector<double>())
	{
		ReducedPrimitiveGradient2D res = slope;
		boost::array<ReducedPrimitive, 20> centroid_vars;
		int n = static_cast<int>(edge_list.size());
		for (size_t i = 0; i<edge_list.size(); ++i)
		{
			centroid_vars[i] = interp_all(cell, cell_trace, center, slope,
				CalcCentroid(edge_list[i]));
		}
		DensityGetter density_getter;
		PressureGetter pressure_getter;
		XVelocityGetter xvelocity_getter;
		YVelocityGetter yvelocity_getter;
		boost::array<PropertyGetter const*, 4> hydro_vars;
		hydro_vars[0] = &density_getter;
		hydro_vars[1] = &pressure_getter;
		hydro_vars[2] = &xvelocity_getter;
		hydro_vars[3] = &yvelocity_getter;
		BOOST_FOREACH(PropertyGetter const* rpg, hydro_vars)
		{
			PropertyGetter const& pg = *rpg;
			double psi = 1;
			boost::container::static_vector<double, 20> neighbor_vals;
			BOOST_FOREACH(Primitive const& pr, neighbors)
			{
				neighbor_vals.push_back(pg.get_property(pr));
			}
			const double my_val = pg.get_property(cell);
			const double min_val = min(*std::min_element(neighbor_vals.begin(),
				neighbor_vals.end()), my_val);
			const double max_val = max(*std::max_element(neighbor_vals.begin(),
				neighbor_vals.end()), my_val);
			const double max_diff = max(max_val - my_val,
				my_val - min_val);
			for (size_t i = 0; i<static_cast<size_t>(n); ++i)
			{
				const double centroid_val = pg.get_property(centroid_vars[i]);
				const double dphi = centroid_val - my_val;
				if (abs(dphi)<0.1*max_diff&&centroid_val*my_val>0)
					continue;
				else if (centroid_val>my_val)
					psi = min(psi, (max_val - my_val) / (centroid_val - my_val));
				else if (centroid_val<my_val)
					psi = min(psi, (min_val - my_val) / (centroid_val - my_val));
				else if (fabs(centroid_val - my_val)<1e-9)
					psi = min(psi, 1.0);
				else
					throw "Something bad happened in LinearGaussConsistent::Prepare";
			}
			pg.get_property(res) *= psi;
		}
		// Deal with tracers
		if (!res.tracers.empty())
		{
			int tracer_number = static_cast<int>(res.tracers.size());
			for (size_t i = 0; i<static_cast<size_t>(tracer_number); ++i)
			{
				double psi = 1;
				boost::container::static_vector<double, 20> neighbor_vals;
				for (size_t j = 0; j<static_cast<size_t>(n); ++j)
					neighbor_vals.push_back(neighbor_tracers[j][i]);
				const double my_val = cell_trace[i];
				const double min_val = min(*std::min_element(neighbor_vals.begin(),
					neighbor_vals.end()), my_val);
				const double max_val = max(*std::max_element(neighbor_vals.begin(),
					neighbor_vals.end()), my_val);
				const double max_diff = max(max_val - my_val,
					my_val - min_val);
				for (size_t k = 0; k<static_cast<size_t>(n); ++k)
				{
					const double centroid_val = centroid_vars[k].tracers[i];
					const double dphi = centroid_val - my_val;
					if (abs(dphi)<0.1*max_diff&&centroid_val*my_val>0)
						continue;
					else
						if (centroid_val>my_val)
							psi = min(psi, (max_val - my_val) / (centroid_val - my_val));
						else
							if (centroid_val<my_val)
								psi = min(psi, (min_val - my_val) /
								(centroid_val - my_val));
							else
								if (fabs(centroid_val - my_val)<1e-9)
									psi = min(psi, 1.0);
								else
								{
									UniversalError eo("Something bad happened in LinearGaussConsistent::Prepare");
									throw eo;
								}
				}
				res.tracers[i] *= psi;
			}
		}
		return res;
	}

	ReducedPrimitiveGradient2D shocked_slope_limit(Primitive const& cell,
		Vector2D const& center, vector<Primitive> const& neighbors,
		vector<Edge> const& edge_list, ReducedPrimitiveGradient2D const& slope,
		double diffusecoeff,
		vector<vector<double> > const& neighbor_tracers = vector<vector<double> >(),
		vector<double> const& cell_trace = vector<double>())
	{
		ReducedPrimitiveGradient2D res = slope;
		boost::array<ReducedPrimitive, 20> centroid_vars;
		int n = static_cast<int>(edge_list.size());
		for (size_t i = 0; i<static_cast<size_t>(n); ++i)
		{
			centroid_vars[i] = interp_all(cell, cell_trace, center, slope,
				CalcCentroid(edge_list[i]));
		}
		DensityGetter density_getter;
		PressureGetter pressure_getter;
		XVelocityGetter xvelocity_getter;
		YVelocityGetter yvelocity_getter;
		boost::array<PropertyGetter const*, 4> hydro_vars;
		hydro_vars[0] = &density_getter;
		hydro_vars[1] = &pressure_getter;
		hydro_vars[2] = &xvelocity_getter;
		hydro_vars[3] = &yvelocity_getter;
		BOOST_FOREACH(PropertyGetter const* rpg, hydro_vars)
		{
			PropertyGetter const& pg = *rpg;
			double psi = 1;
			boost::container::static_vector<double, 20> neighbor_vals;
			BOOST_FOREACH(Primitive const& pr, neighbors)
			{
				neighbor_vals.push_back(pg.get_property(pr));
			}
			const double my_val = pg.get_property(cell);
			double maxdiff = 0;
			BOOST_FOREACH(double nv, neighbor_vals)
			{
				maxdiff = max(maxdiff, static_cast<double>(abs(nv - my_val)));
			}
			for (size_t i = 0; i<edge_list.size(); ++i)
			{
				const double centroid_val = pg.get_property(centroid_vars[i]);
				const double dphi = centroid_val - my_val;
				if (abs(dphi)<0.1*maxdiff&&centroid_val*my_val>0)
					continue;
				else if (dphi>0)
					psi = min(psi, diffusecoeff*max((neighbor_vals[i] - my_val) / dphi, 0.0));
				else if (dphi<0)
					psi = min(psi, diffusecoeff*max((neighbor_vals[i] - my_val) / dphi, 0.0));
				else if (fabs(dphi)<1e-9)
					psi = min(psi, 1.0);
				else
					throw "Something bad happened in LinearGaussConsistent::Prepare";
			}
			pg.get_property(res) *= psi;
		}

		// Deal with tracers
		if (!res.tracers.empty())
		{
			int tracer_number = static_cast<int>(res.tracers.size());
			for (size_t i = 0; i<static_cast<size_t>(tracer_number); ++i)
			{
				double psi = 1;
				boost::container::static_vector<double, 20> neighbor_vals;
				for (size_t j = 0; j<static_cast<size_t>(n); ++j)
					neighbor_vals.push_back(neighbor_tracers[j][i]);
				//const double my_val = cell_trace[i];
				const double my_val2 = cell_trace[i];
				double maxdiff = 0;
				BOOST_FOREACH(double nv, neighbor_vals)
				{
					maxdiff = max(maxdiff, static_cast<double>(abs(nv - my_val2)));
				}
				for (size_t j = 0; j<static_cast<size_t>(n); ++j)
				{
					const double centroid_val = centroid_vars[j].tracers[i];
					const double dphi = centroid_val - my_val2;
					if (abs(dphi)<0.1*maxdiff&&centroid_val*my_val2>0)
						continue;
					else
						if (dphi>0)
							psi = min(psi, diffusecoeff*
							max((neighbor_vals[i] - my_val2) / dphi, 0.0));
						else
							if (dphi<0)
								psi = min(psi, diffusecoeff*
								max((neighbor_vals[i] - my_val2) / dphi, 0.0));
							else
								if (fabs(dphi)<1e-9)
									psi = min(psi, 1.0);
								else
									throw "Something bad happened in LinearGaussConsistent::Prepare";
				}
				res.tracers[i] *= psi;
			}
		}
		return res;
	}

	vector<Edge> GetEdgeList(Tessellation const& tess,
		vector<int> const& edge_indices)
	{
		vector<Edge> res(edge_indices.size());
		for (size_t i = 0; i<edge_indices.size(); ++i){
			res[i] = tess.GetEdge(edge_indices[i]);
		}
		return res;
	}

	double PressureRatio(Primitive cell, vector<Primitive> const& neigh)
	{
		int n = static_cast<int>(neigh.size());
		double res = 1;
		double p = cell.Pressure;
		for (size_t i = 0; i<static_cast<size_t>(n); ++i)
		{
			if (p>neigh[i].Pressure)
				res = min(res, neigh[i].Pressure / p);
			else
				res = min(res, p / neigh[i].Pressure);
		}
		return res;
	}

	bool is_shock(ReducedPrimitiveGradient2D const& naive_slope,
		double cell_width,
		double shock_ratio,
		Primitive const& cell,
		vector<Primitive> const& neighbor_list,
		double pressure_ratio)
	{
		const bool cond1 =
			(naive_slope.xvelocity.x + naive_slope.yvelocity.y)*
			cell_width<(-shock_ratio)*cell.SoundSpeed;
		const bool cond2 = PressureRatio(cell, neighbor_list)<pressure_ratio;
		return cond1 || cond2;
	}

	ReducedPrimitiveGradient2D calc_slope(Tessellation const& tess,
		vector<Primitive> const& cells, vector<vector<double> > const& tracers,
		vector<bool> const& isrelevant, int cell_index, bool slf,
		OuterBoundary const& obc, HydroBoundaryConditions const& hbc,
		double shockratio, double diffusecoeff, double pressure_ratio, double time,
		bool /*rigidflag*/)
	{
		vector<int> edge_indices = tess.GetCellEdges(cell_index);
		vector<Edge> edge_list = GetEdgeList(tess, edge_indices);
		vector<Vector2D> neighbor_mesh_list = GetNeighborMesh(tess, edge_list,
			cell_index, obc);
		vector<Vector2D> neighbor_cm_list = GetNeighborCM(tess, edge_list,
			cell_index, obc);
		vector<Primitive> neighbor_list = GetNeighborPrimitive(tess, edge_list, cell_index,
			obc, cells, hbc, time, isrelevant);
		vector<vector<double> > neighbor_tracers;
		if (!tracers.empty())
			neighbor_tracers = GetNeighborTracers(tess, edge_list, cell_index, tracers,
			hbc, time);
		ReducedPrimitiveGradient2D naive_slope, s_compare;
		if (!tracers.empty())
		{
			naive_slope = calc_naive_slope
				(cells[static_cast<size_t>(cell_index)], tess.GetMeshPoint(cell_index), tess.GetCellCM(cell_index),
				tess.GetVolume(cell_index), neighbor_list,
				neighbor_mesh_list, neighbor_cm_list, edge_list, tracers[static_cast<size_t>(cell_index)], neighbor_tracers);
		}
		else
		{
			naive_slope = calc_naive_slope
				(cells[static_cast<size_t>(cell_index)], tess.GetMeshPoint(cell_index), tess.GetCellCM(cell_index),
				tess.GetVolume(cell_index), neighbor_list,
				neighbor_mesh_list, neighbor_cm_list, edge_list);
		}
		if (slf)
		{
			if (!is_shock(naive_slope, tess.GetWidth(cell_index), shockratio,
				cells[static_cast<size_t>(cell_index)], neighbor_list, pressure_ratio))
			{
				if (!tracers.empty())
					return slope_limit(cells[static_cast<size_t>(cell_index)], tess.GetCellCM(cell_index),
					neighbor_list, edge_list, naive_slope, neighbor_tracers,
					tracers[static_cast<size_t>(cell_index)]);
				else
					return slope_limit(cells[static_cast<size_t>(cell_index)], tess.GetCellCM(cell_index),
					neighbor_list, edge_list, naive_slope);
			}
			else
			{
				if (!tracers.empty())
					return shocked_slope_limit(cells[static_cast<size_t>(cell_index)],
					tess.GetCellCM(cell_index), neighbor_list, edge_list,
					naive_slope, diffusecoeff, neighbor_tracers, tracers[static_cast<size_t>(cell_index)]);
				else
					return shocked_slope_limit(cells[static_cast<size_t>(cell_index)],
					tess.GetCellCM(cell_index), neighbor_list, edge_list,
					naive_slope, diffusecoeff);
			}
		}
		else
		{
			return naive_slope;
		}
	}
}

void LinearGaussImproved::Prepare(Tessellation const& tessellation,
	vector<Primitive> const& cells, vector<vector<double> > const& tracers,
	vector<bool> const& isrelevant, double /*dt*/, double time)
{
	if (tessellation.GetPointNo() != static_cast<int>(rslopes_.size()))
	{
		rslopes_.resize(static_cast<size_t>(tessellation.GetPointNo()));
	}
	for (size_t i = 0; i<static_cast<size_t>(tessellation.GetPointNo()); ++i)
	{
		if (!hbc_.IsGhostCell(static_cast<int>(i), tessellation))
			rslopes_[i] = calc_slope(tessellation, cells, tracers, isrelevant, static_cast<int>(i), slf_, obc_, hbc_,
			shockratio_, diffusecoeff_, pressure_ratio_, time, _rigidflag);
	}
}

namespace {
	Primitive CalcDtFlux(Primitive const&cell,
		ReducedPrimitiveGradient2D const& grad, double dt, Vector2D const& vface)
	{
		Primitive res;
		res.Density -= 0.5*dt*((cell.Velocity.x - vface.x)*grad.density.x +
			(cell.Velocity.y - vface.y)*grad.density.y +
			cell.Density*(grad.xvelocity.x + grad.yvelocity.y));
		res.Velocity.x -= 0.5*dt*((cell.Velocity.x - vface.x)*
			grad.xvelocity.x + (cell.Velocity.y - vface.y)*
			grad.yvelocity.x + grad.pressure.x / cell.Density);
		res.Velocity.y -= 0.5*dt*((cell.Velocity.x - vface.x)*
			grad.xvelocity.y + (cell.Velocity.y - vface.y)*
			grad.yvelocity.y + grad.pressure.y / cell.Density);
		res.Pressure -= 0.5*dt*(cell.Pressure*(1 + cell.Pressure /
			(cell.Energy*cell.Density))*(grad.xvelocity.x +
			grad.yvelocity.y) + grad.pressure.x*(cell.Velocity.x - vface.x)
			+ grad.pressure.y*(cell.Velocity.y - vface.y));
		return res;
	}

	vector<double> CalcDtFlux(vector<double> const& tracer, Primitive const& cell,
		ReducedPrimitiveGradient2D const& grad, double dt, Vector2D const& vface)
	{
		vector<double> res(tracer.size());
		//		int n=static_cast<int>(tracer.size());
		for (size_t i = 0; i<tracer.size(); ++i)
			res[i] = -0.5*dt*(grad.tracers[i].x*(cell.Velocity.x - vface.x) +
			grad.tracers[i].y*(cell.Velocity.y - vface.y));
		return res;
	}
}

vector<double> LinearGaussImproved::interpolateTracers
(Tessellation const& tess, vector<Primitive> const& cells,
vector<vector<double> > const& tracers, double dt, Edge const& edge,
int side, InterpolationType interptype, Vector2D const& vface) const
{
	int cell_index = pair_member(edge.neighbors, side);
	Vector2D target = CalcCentroid(edge);
	if (interptype == InBulk)
	{
		vector<double> res = interp_tracer(tracers[static_cast<size_t>(cell_index)],
			tess.GetCellCM(cell_index), rslopes_[static_cast<size_t>(cell_index)], target);
		if (soitf_)
		{
			vector<double> temp = CalcDtFlux(tracers[static_cast<size_t>(cell_index)], cells[static_cast<size_t>(cell_index)],
				rslopes_[static_cast<size_t>(cell_index)], dt, vface);
			int n = static_cast<int>(temp.size());
			for (size_t i = 0; i<static_cast<size_t>(n); ++i)
				res[i] += temp[i];
		}
		return res;
	}
	else
		if (interptype == Boundary)
		{
			int other = pair_member(edge.neighbors, (side + 1) % 2);
			vector<double> res = interp_tracer(tracers[static_cast<size_t>(other)],
				tess.GetMeshPoint(other), rslopes_[static_cast<size_t>(other)],
				target);
			if (soitf_)
			{
				vector<double> temp = CalcDtFlux(tracers[static_cast<size_t>(other)], cells[static_cast<size_t>(other)], rslopes_[static_cast<size_t>(other)],
					dt, vface);
				int n = static_cast<int>(temp.size());
				for (size_t i = 0; i<static_cast<size_t>(n); ++i)
					res[i] += temp[i];
			}
			return res;
		}
		else
			throw UniversalError("Wrong interpolation type in linear_gauss_scalar");
}

namespace
{
	UniversalError InterpolationError(UniversalError &eo, int cell_index, Edge const& edge,
		Primitive const& cell, ReducedPrimitiveGradient2D const& slope, Vector2D const& dr)
	{
		eo.AddEntry("Error in interpolation of cell", cell_index);
		eo.AddEntry("Error in interpolation of edge, n0 ", edge.neighbors.first);
		eo.AddEntry("Error in interpolation of edge, n1 ", edge.neighbors.second);
		eo.AddEntry("Density", cell.Density);
		eo.AddEntry("Pressure", cell.Pressure);
		eo.AddEntry("Density slope x", slope.density.x);
		eo.AddEntry("Density slope y", slope.density.y);
		eo.AddEntry("Pressure slope x", slope.pressure.x);
		eo.AddEntry("Pressure slope y", slope.pressure.y);
		eo.AddEntry("Interpolation dx", dr.x);
		eo.AddEntry("Interpolation dy", dr.y);
		throw eo;
	}
}

Primitive LinearGaussImproved::Interpolate(Tessellation const& tess,
	vector<Primitive> const& cells, double dt, Edge const& edge, int side,
	InterpolationType interptype, Vector2D const& vface) const
{
	int cell_index = pair_member(edge.neighbors, side);
	Vector2D target = CalcCentroid(edge);
	if (interptype == InBulk)
	{
		const Primitive temp = interp_primitive(cells[static_cast<size_t>(cell_index)],
			tess.GetCellCM(cell_index), rslopes_[static_cast<size_t>(cell_index)], target);
		try
		{
			Primitive res = CalcPrimitive(temp.Density, temp.Pressure,
				temp.Velocity, eos_);
			if (soitf_)
			{
				res += CalcDtFlux(cells[static_cast<size_t>(cell_index)], rslopes_[static_cast<size_t>(cell_index)], dt, vface);
				res = CalcPrimitive(res.Density, res.Pressure, Vector2D(res.Velocity.x,
					res.Velocity.y), eos_);
			}
			return res;
		}
		catch (UniversalError &eo)
		{
			throw InterpolationError(eo, cell_index, edge, cells[static_cast<size_t>(cell_index)],
				rslopes_[static_cast<size_t>(cell_index)], target - tess.GetCellCM(cell_index));
		}
	}
	else
		if (interptype == Boundary)
		{
			const int other = pair_member(edge.neighbors, (side + 1) % 2);
			const Primitive temp = interp_primitive(cells[static_cast<size_t>(other)], tess.GetMeshPoint(other),
				rslopes_[static_cast<size_t>(other)], target);
			try
			{
				Primitive res = CalcPrimitive(temp.Density, temp.Pressure, temp.Velocity, eos_);
				if (soitf_)
				{
					res += CalcDtFlux(cells[static_cast<size_t>(other)], rslopes_[static_cast<size_t>(other)], dt, vface);
					res = CalcPrimitive(res.Density, res.Pressure, res.Velocity, eos_);
				}
				return res;
			}
			catch (UniversalError &eo)
			{
				throw InterpolationError(eo, other, edge, cells[static_cast<size_t>(other)],
					rslopes_[static_cast<size_t>(other)], target - tess.GetMeshPoint(other));
			}
		}
		else
			throw UniversalError("Wrong interpolation type");
}

LinearGaussImproved::~LinearGaussImproved(void) {}

vector<ReducedPrimitiveGradient2D>& LinearGaussImproved::GetGradients(void)
{
	return rslopes_;
}
